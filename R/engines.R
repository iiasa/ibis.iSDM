#' @include bdproto-engine.R utils-inla.R bdproto-distributionmodel.R
NULL

#' Use INLA as engine
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param optional_mesh A directly supplied [`INLA`] mesh (Default: NULL)
#' @param max_distance Range of maximum distance between two nodes is between 50 and 5000 meter
#' @param offset Offset for INLA mesh
#' @param ... Other variables
#' @name engine_inla
NULL
#' @import INLA
#' @rdname engine_inla
#' @export
engine_inla <- function(x, optional_mesh = NULL, max_distance = c(10,1000), offset = c(1,1), verbose = TRUE,...) {
  # TODO:
  # Find a better way to pass on parameters...
  # assert that arguments are valid
  assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                          inherits(x$background,'Raster'),
                          inherits(optional_mesh,'inla.mesh') || is.null(optional_mesh),
                          is.vector(max_distance),
                          assertthat::is.flag(verbose),
                          requireNamespace("INLA", quietly = TRUE))

  # Set the projection mesh
  if(inherits(optional_mesh,'inla.mesh')) {
    mesh <- optional_mesh
  } else {
    # Background points
    dat <- raster::rasterToPoints(x$background)[,c('x','y')]

    # Make a boundary from the background
    bounds <- raster::boundaries(background, type = 'outer', asNA = TRUE)

    bdry <- INLA::as.inla.mesh.segment(rasterToPolygons(bounds),join = TRUE)
    bdry$loc <- inla.mesh.map(bdry$loc)

    # Prepare the mesh
    suppressWarnings(
      mesh <- INLA::inla.mesh.2d(
        loc = dat[,c("x", "y")], # Define initial triangulation points
        max.edge = max_distance, # maximum distance between two nodes
        # FIXME: This still does not work correctly
#        boundary = bdry, # boundary
        offset = offset, # Offset of outer boundaries
        crs = inla.CRS(projargs = sp::proj4string(x$background))
      )
    )
    rm(dat)
  }

  # Print a message in case there is already an engine object
  if(!is.Waiver(x$engine)) message('Replacing currently selected engine.')

  # Set engine in distribution object
  x$set_engine(
    bdproto(
      "INLA-Engine",
      Engine,
      name = "<INLA>",
      data = list(
        'mesh' = mesh
      ),
      # parameters = parameters(
      #   numeric_parameter("gap", gap, lower_limit = 0),
      #   integer_parameter("time_limit", time_limit, lower_limit = -1L,
      #                     upper_limit = as.integer(.Machine$integer.max)),
      #   integer_parameter("presolve", presolve, lower_limit = -1L,
      #                     upper_limit = 2L),
      #   integer_parameter("threads", threads, lower_limit = 1L,
      #                     upper_limit = parallel::detectCores(TRUE)),
      #   binary_parameter("first_feasible", first_feasible),
      #   binary_parameter("numeric_focus", numeric_focus),
      #   binary_parameter("verbose", verbose)
      #   ),
      # Spatial latent function
      calc_latent_spatial = function(self, alpha = 2,...){
        # Define Matern SPDE model and save
        self$data$latentspatial <- inla.spde2.matern(self$data$mesh,alpha = alpha,...)
        # Make index for spatial field
        self$data$s.index <- inla.spde.make.index(name = "spatial.field",
                                                  n.spde = self$data$latentspatial$n.spde)
        assertthat::assert_that(
          inherits(self$data$latentspatial,'inla.spde'),
          length(self$data$s.index$spatial.field) == self$data$mesh$n
        )
        invisible()
      },
      # Get latent spatial equation bit
      get_equation_latent_spatial = function(self,spatial_object = 'spde'){
        if( 'latentspatial' %notin% names(self$data) ) self$get_equation_latent_spatial()
        return(
          paste0('f(spatial.field, model = ',spatial_object,')')
        )
      },
      # Setup computation function
      setup = function(self, model, ...){
        # TODO: Some assert calls
        # TODO: Potentially add a setup log later for this function

        # For poipo
        # Create INLA projection matrix from mesh to observed data from the points
        mat_proj <- INLA::inla.spde.make.A(
          mesh = self$get_data('mesh'),
          loc = as.matrix(model$data$poipo_values[,c('x','y')])
        )

        # Create INLA stack
        # The three main inla.stack() arguments are a vector list with the data (data),
        # a list of projector matrices (each related to one block effect,
        # A) and the list of effects (effects).

        # Response for inla stack
        resp <- model$data$poipo_response
        ll_resp <- list()
        ll_resp[[ resp ]] <- model$data$poipo[,resp]
        # A column for the offset
        ll_resp[[ 'e' ]] <- rep(0, nrow(model$data$poipo_values) )

        # Effects matrix
        ll_effects <- list()
        pred_names <- model$predictors_names
        # Note, order adding this is important apparently...
        ll_effects[['intercept']] <- list(intercept = rep(1, self$get_data('mesh')$n ))
        ll_effects[['predictors']] <- model[['data']][['poipo_values']][,pred_names] # Get only the covariates to use
        # Add latent
        if(x$get_latent()=='<Spatial>'){
          if('latentspatial' %notin% names(self$data)) self$calc_latent_spatial()
          # Get spatial object
          spde <- self$get_data('latentspatial')
          iset <- self$get_data('s.index')
          ll_effects[['spatial.field']] <- iset#list(Bnodes = 1:spde$n.spde)
          # Define projection matrix
          # FIXME: Check that the below formulation is correct!
          A = list(mat_proj,1,mat_proj)
        } else {
          spde <- NULL # Set SPDE to NULL
          A = list(mat_proj, 1)
        }

        # Create the stack for the fitted model
        stk_resp <-
          INLA::inla.stack(
            data =  ll_resp,             # Response
            A    =  A,                   # Predictor projection matrix. 1 is included to make a list
            effects = ll_effects,        # Effects matrix
            tag = paste0('obs_','poipo') # Description tag
          )
        # Save this stack in data
        self$set_data('stk_resp', stk_resp)

        # Create and join prediction stack
        stk_pred <- inla_make_prediction_stack(stk_resp = stk_resp,
                                               cov = model$predictors,
                                               mesh = self$get_data('mesh'),
                                               type = 'poipo',
                                               spde = spde)
        # Save this stack in data
        self$set_data('stk_pred', stk_pred)
        # And also the joined stack
        self$set_data('stk_full', INLA::inla.stack(stk_resp, stk_pred) )

        invisible()
      },
      # Main INLA training function ----
      train = function(self, model, varsel = FALSE) {
        # TODO: Implement variable selection
        # Check that all inputs are there
        assertthat::assert_that(
          is.list(model),length(model)>1,
          'stk_resp' %in% names(self$data), inherits(self$get_data('stk_resp'),'inla.data.stack'),
          'stk_pred' %in% names(self$data), inherits(self$get_data('stk_pred'),'inla.data.stack'),
          'stk_full' %in% names(self$data), inherits(self$get_data('stk_full'),'inla.data.stack')
        )
        # Get the datasets
        if(x$get_latent()=='<Spatial>'){
            stack_data_resp <- INLA::inla.stack.data(self$get_data('stk_resp'),spde = self$get_data('latentspatial'))
            stack_data_full <- INLA::inla.stack.data(self$get_data('stk_full'),spde = self$get_data('latentspatial'))
          } else {
            stack_data_resp <- INLA::inla.stack.data(self$get_data('stk_resp'))
            stack_data_full <- INLA::inla.stack.data(self$get_data('stk_full'))
          }

        # Train the model on the response
        fit_resp <- INLA::inla(formula = model$equation$poipo, # The specified formula
                          data  = stack_data_resp,  # The data stack
                          quantiles = c(0.05, 0.5, 0.95),
                          E = INLA::inla.stack.data(self$get_data('stk_resp'))$e, # Exposure (Eta) for Poisson model
                          family= 'poisson',   # Family the data comes from
                          control.family = list(link = "log"), # Control options
                          control.predictor=list(A = INLA::inla.stack.A(self$get_data('stk_resp')),link = NULL, compute = TRUE),  # Compute for marginals of the predictors
                          control.compute = list(cpo = TRUE, waic = TRUE, config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                          verbose = FALSE, # To see the log of the model runs
                          control.inla(int.strategy = "eb", # Empirical bayes for integration
                                       strategy = 'simplified.laplace', huge = TRUE), # To make it run faster...
                          num.threads = parallel::detectCores()-1
        )

        # Predict on full
        fit_pred <- INLA::inla(formula = model$equation$poipo, # The specified formula
                               data  = stack_data_full,  # The data stack
                               quantiles = c(0.05, 0.5, 0.95),
                               E = INLA::inla.stack.data(self$get_data('stk_full'))$e,
                               family= 'poisson',   # Family the data comes from
                               control.family = list(link = "log"), # Control options
                               control.predictor = list(A = INLA::inla.stack.A(self$get_data('stk_full')), link = NULL, compute = TRUE),  # Compute for marginals of the predictors
                               control.compute = list(cpo = FALSE, waic = FALSE, config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                               verbose = FALSE, # To see the log of the model runs
                               control.inla(int.strategy = "eb", # Empirical bayes for integration
                                            strategy = 'simplified.laplace'), # To make it run faster...
                               num.threads = parallel::detectCores() - 1
        )

        # Create a spatial prediction
        index.pred <- INLA::inla.stack.index(self$get_data('stk_full'), 'pred_poipo')$data
        post <- fit_pred$summary.linear.predictor[index.pred, ]
        assertthat::assert_that(nrow(post)>0)
        # TODO: Wrap Backtransform into a distribution dependent function
        post[,c('mean','0.05quant','0.5quant','0.95quant','mode')] <-
          exp(
            post[,c('mean','0.05quant','0.5quant','0.95quant','mode')]
            )
        post <- subset(post, select = c('mean','sd','0.05quant','0.5quant','0.95quant','mode','kld') )
        names(post) <- make.names(names(post))
        # Fill output rasters
        prediction <- fill_rasters(post = post,background = model$background)
        prediction <- raster::mask(prediction, model$background) # Mask with background

        # Definition of INLA Model object ----
        out <- bdproto(
              "INLA-Model",
              DistributionModel,
              name = model$runname,
              id = model$id,
              fits = list(
                "fit_best" = fit_resp,
                "fit_pred" = fit_pred,
                "mesh"     = self$get_data('mesh'),
                "prediction" = prediction
                )
        )
        return(out)
      }
      ))
}

# TODO:
# Other engines, e.g. STAN or simple estimators such as Maxent?
