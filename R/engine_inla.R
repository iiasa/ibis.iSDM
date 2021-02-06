#' @include bdproto-engine.R utils-inla.R bdproto-distributionmodel.R
NULL

#' Use INLA as engine
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param optional_mesh A directly supplied [`INLA`] mesh (Default: NULL)
#' @param max_distance Range of maximum distance between two nodes
#' @param offset Offset for INLA mesh
#' @param ... Other variables
#' @name engine_inla
NULL
#' @import INLA
#' @rdname engine_inla
#' @export
engine_inla <- function(x, optional_mesh = NULL,
                        max.edge = c(1,5),
                        offset = c(1,1),
                        cutoff = 0.1,
                        ...) {
  # TODO:
  # Find a better way to pass on parameters such as those related to the mesh size...
  # assert that arguments are valid
  assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                          inherits(x$background,'sf'),
                          inherits(optional_mesh,'inla.mesh') || is.null(optional_mesh),
                          is.vector(max.edge),
                          is.vector(offset) || is.numeric(offset),
                          is.numeric(cutoff),
                          requireNamespace("INLA", quietly = TRUE))

  # Convert the study region
  region.poly <- as(sf::st_geometry(x$background), "Spatial")

  # Set the projection mesh
  if(inherits(optional_mesh,'inla.mesh')) {
    # Load a provided on
    mesh <- optional_mesh
  } else {
    # Create a new one
    # Convert to boundary object for later
    suppressWarnings(
      bdry <- INLA::inla.sp2segment(
        sp = region.poly,
        join = TRUE,
        crs = INLA::inla.CRS(projargs = sp::proj4string(region.poly))
      )
    )
    bdry$loc <- INLA::inla.mesh.map(bdry$loc)

    # --- #
    # FIXME: Move this to the man description above
    # Create the mesh
    # A good mesh needs to have triangles as regular as possible in size and shape: equilateral
    suppressWarnings(
      mesh <- INLA::inla.mesh.2d(
        # Boundary object
        boundary = bdry,
        # The largest allowed triangle edge length, must be in the same scale units as the coordinates
        # Lower bounds affect the density of triangles
        max.edge = max.edge,
        # The automatic extension distance.
        # If positive: same scale units.
        # If negative, interpreted as a factor relative to the approximate data diameter;
        #   i.e., a value of -0.10 will add a 10% of the data diameter as outer extension.
        offset = offset,
        # The minimum allowed distance between points,
        # it means that points at a closer distance than the supplied value are replaced by a single vertex.
        # it is critical when there are some points very close to each other,
        #   either for point locations or in the domain boundary.
        cutoff = cutoff,
        # Define the CRS
        crs = bdry$crs
      )
    )
  }

  # Calculate area in km²
  ar <- suppressWarnings(
    mesh_area(mesh = mesh,region.poly = region.poly, variant = 'gpc')
  )

  # Print a message in case there is already an engine object
  if(!is.Waiver(x$engine)) message('Replacing currently selected engine.')

  # Set engine in distribution object
  x$set_engine(
    bdproto(
      "INLA-Engine",
      Engine,
      name = "<INLA>",
      data = list(
        'mesh' = mesh,
        'mesh.area' = ar
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
      calc_latent_spatial = function(self,type = 'pc', alpha = 2,
                                     prior.range = c(5,0.001),
                                     prior.sigma = c(0.5,0.05),
                                     ...){
        if(type=='pc'){
          # Define PC Matern SPDE model and save
          self$data$latentspatial <- INLA::inla.spde2.pcmatern(
            self$data$mesh,
            alpha = alpha,
            # P(Range < 1°) = 0.001  and P(sigma > 0.5) = 0.05
            prior.range = prior.range,
            prior.sigma = prior.sigma)
        } else {
          self$data$latentspatial <- INLA::inla.spde2.matern(mesh = self$data$mesh, alpha = alpha)
        }
        # Make index for spatial field
        self$data$s.index <- INLA::inla.spde.make.index(name = "spatial.field",
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

        # Calculate the offset
        e <- rep(0, nrow(model$data$poipo_values) )
        # Calculate e as relative area from the points
        # sfmesh <- mesh_as_sf( self$get_data('mesh'))
        # # Set those not intersecting with the background to 0
        # ind <- suppressMessages(
        #   sf::st_join(sfmesh, x$background, join = st_within)[,3] %>% st_drop_geometry()
        # )
        # sfmesh[which( is.na(ind[,1]) ),'relarea'] <- NA
        # e <- suppressMessages(
        #   point_in_polygon(
        #     poly = sfmesh,
        #     points = model$data$poipo
        #   )$relarea
        # )

        # Create INLA stack
        # The three main inla.stack() arguments are a vector list with the data (data),
        # a list of projector matrices (each related to one block effect,
        # A) and the list of effects (effects).

        # Response for inla stack
        resp <- model$data$poipo_response
        ll_resp <- list()
        ll_resp[[ resp ]] <- cbind( model$data$poipo[,resp] )
        # A column for the offset
        ll_resp[[ 'e' ]] <- e

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
          ll_effects[['spatial.field']] <- list(spatial.field = 1:spde$n.spde) #iset
          # Define projection matrix, has to match ll_effects
          # FIXME: Check that the below formulation is correct!
          A = list(mat_proj,1,mat_proj)
        } else {
          spde <- NULL # Set SPDE to NULL
          A = list(mat_proj, 1)
        }

        # Create the stack for the fitted model
        stk_obs <-
          INLA::inla.stack(
            data =  ll_resp,             # Response
            A    =  A,                   # Predictor projection matrix.
            effects = ll_effects,        # Effects matrix
            tag = paste0('obs_','poipo') # Description tag
          )
        # Save this stack in data
        self$set_data('stk_obs', stk_obs)

        # Create and join prediction stack
        stk_pred <- inla_make_prediction_stack(stk_resp = stk_obs,
                                               cov = model$predictors,
                                               mesh = self$get_data('mesh'),
                                               type = 'poipo',
                                               spde = spde)
        # Save this stack in data
        self$set_data('stk_pred', stk_pred)
        # And also the joined stack
        self$set_data('stk_full', INLA::inla.stack(stk_obs, stk_pred) )

        invisible()
      },
      # Main INLA training function ----
      train = function(self, model, varsel = TRUE) {
        # Check that all inputs are there
        assertthat::assert_that(
          is.list(model),length(model)>1,
          'stk_obs' %in% names(self$data), inherits(self$get_data('stk_obs'),'inla.data.stack'),
          'stk_pred' %in% names(self$data), inherits(self$get_data('stk_pred'),'inla.data.stack'),
          'stk_full' %in% names(self$data), inherits(self$get_data('stk_full'),'inla.data.stack')
        )
        # Get the datasets
        if(x$get_latent()=='<Spatial>'){
            stack_data_resp <- INLA::inla.stack.data(self$get_data('stk_obs'),spde = self$get_data('latentspatial'))
            stack_data_full <- INLA::inla.stack.data(self$get_data('stk_full'),spde = self$get_data('latentspatial'))
          } else {
            stack_data_resp <- INLA::inla.stack.data(self$get_data('stk_obs'))
            stack_data_full <- INLA::inla.stack.data(self$get_data('stk_full'))
          }
        # ----------- #
        # Provided or default formula
        master_form <- model$equation$poipo

        # Perform variable selection
        if(varsel){
          message('Performing variable selection...')
          #
          if(x$get_latent()=='<Spatial>') speq <- self$get_equation_latent_spatial() else speq <- NULL
          # Get formula list
          lf <- formula_combinations(varnames = model$predictors_names,
                               response = all.vars(master_form)[1],
                               InterceptOnly = FALSE,
                               spde_term = speq,
                               type = 'forward'
                               )

          # Set progress bar
          pb <- progress::progress_bar$new(total = length(lf))

          results <- data.frame(stringsAsFactors = FALSE)
          # Now loop through
          for(k in 1:length(lf)){
            pb$tick()

            # Train the model on the response
            fit <- INLA::inla(formula = as.formula(lf[k]), # The specified formula
                                   data  = stack_data_resp,  # The data stack
                                   E = INLA::inla.stack.data(self$get_data('stk_obs'))$e, # Exposure (Eta) for Poisson model
                                   family= 'poisson',   # Family the data comes from
                                   control.family = list(link = "log"), # Control options
                                   control.predictor=list(A = INLA::inla.stack.A(self$get_data('stk_obs')),
                                                          link = NULL, compute = TRUE),  # Compute for marginals of the predictors
                                   control.compute = list(cpo = TRUE,dic = TRUE, waic = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                                   control.inla(int.strategy = "eb", # Empirical bayes for integration
                                                strategy = 'simplified.laplace', huge = TRUE), # To make it run faster...
                                   num.threads = parallel::detectCores() - 1
            )
            # Add results
            results <- rbind(results,
                             data.frame(form = lf[k],
                                        converged = fit$ok,
                                        waic = fit$waic$waic,
                                        dic = fit$dic$dic,
                                        # conditional predictive ordinate values
                                        cpo = sum(log(fit$cpo$cpo)) * -2,
                                        mean.deviance = fit$dic$mean.deviance ) )
            rm(fit)
          }
          rm(pb)
          # Determine best model and set of variables
          results <- results[order(results$dic,decreasing = FALSE),]
          results$diff_dic <- c(0, diff(results$dic) )
          # Set new master form
          master_form <- as.formula(
            as.character(results$form[results$diff_dic==0])
          )
        }

        # Train the model on the response
        fit_resp <- INLA::inla(formula = master_form, # The specified formula
                          data  = stack_data_resp,  # The data stack
                          quantiles = c(0.05, 0.5, 0.95),
                          E = INLA::inla.stack.data(self$get_data('stk_obs'))$e, # Exposure (Eta) for Poisson model
                          family= 'poisson',   # Family the data comes from
                          control.family = list(link = "log"), # Control options
                          control.predictor=list(A = INLA::inla.stack.A(self$get_data('stk_obs')),link = NULL, compute = TRUE),  # Compute for marginals of the predictors
                          control.compute = list(cpo = TRUE, waic = TRUE, config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                          verbose = FALSE, # To see the log of the model runs
                          control.inla(int.strategy = "eb", # Empirical bayes for integration
                                       strategy = 'simplified.laplace', huge = TRUE), # To make it run faster...
                          num.threads = parallel::detectCores()-1
        )

        # Predict on full
        fit_pred <- INLA::inla(formula = master_form, # The specified formula
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

        # Create an empty raster from the predictor coordinates
        ra <- rasterFromXYZ(xyz = model$predictors[,c('x','y')],
                            crs = raster::projection(x$background)
                            )
        # Fill output rasters
        prediction <- fill_rasters(post = post,background = ra)
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
