#' @include utils.R bdproto-biodiversitydistribution.R utils-spatial.R
NULL

#' Train the model from a given engine
#'
#' Train a [distribution()] model with the specified engine.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object).
#' @param runname A [`character`] name of the trained run
#' @param ... further arguments passed on.
#'
#' @details
#'
#' @return A distribution prediction object
#' @name train
#' @exportMethod train
#' @aliases train, train-method
#' @export
NULL

#' @name train
#' @rdname train
#' @exportMethod train
#' @export
methods::setGeneric(
  "train",
  signature = methods::signature("x", "runname"),
  function(x, runname,...) standardGeneric("train"))

#' @name train
#' @rdname train
#' @usage \S4method{train}{BiodiversityDistribution}(x)
methods::setMethod(
  "train",
  methods::signature(x = "BiodiversityDistribution", runname = "character"),
  function(x, runname, ...) {
    # Make load checks
    assertthat::assert_that(
      inherits(x, "BiodiversityDistribution"),
      is.character(runname)
    )
    # Now make checks on completeness of the object
    assertthat::assert_that(!is.Waiver(x$engine),
                            msg = 'No engine set for training the distribution model.')
    assertthat::assert_that( x$show_biodiversity_length() > 0,
                             msg = 'No biodiversity data specified.')

    # --- #
    # Set model object for fitting
    model <- list()

    # Set model name
    model[['runname']] <- 'test' #runname

    # Specify a unique id for the run
    model[['id']] <- new_id()

    # Get biodiversity data
    model[['data']] <- list()
    types <- names( x$biodiversity$get_types() )
    for(ty in types) model[['data']][[ty]] <- x$biodiversity$get_data(ty)

    # Get covariates
    if(is.Waiver(x$get_predictor_names())) {
      # Dummy covariate of background raster
      dummy <- x$background; names(dummy) <- 'dummy'
      model[['predictors']] <- dummy
      } else { model[['predictors']] <- x$predictors$get_data(df = TRUE, na.rm = FALSE) }

    # assign default priors
    if(is.Waiver( x$priors )){
      # TODO: Define prior objects. Also look at PC priors https://www.tandfonline.com/doi/abs/10.1080/01621459.2017.1415907?journalCode=uasa20
      message('TODO: Define prior objects')
      model[['priors']] <- NULL
      #x$set_prior( add_default_priors() )
    }
    # Get latent variables
    if(!is.Waiver(x$latentfactors)){
      # Calculate latent spatial factor (saved in engine data)
      if(x$get_latent()=="<Spatial>") x$engine$calc_latent_spatial()
    }

    # Format formulas
    model[['equation']] <- list()
    types <- names( x$biodiversity$get_types() )
    for(ty in types) {
      # Default equation found
      if(x$biodiversity$get_equations()[[ty]]=='<Default>'){
        # Construct formula with all variables
        f <- formula(
                  paste( x$biodiversity$get_columns_occ()[[ty]], '~ ',
                         0, #ifelse(x$show_biodiversity_length()==1,1,0),
                         ' +',
                         paste( x$get_predictor_names(), collapse = ' + ' )  )
                )
        if(x$get_latent()=="<Spatial>"){
          # Update with spatial term
          f <- update.formula(f, paste0(" ~ . + ",x$engine$get_equation_latent_spatial() ) )
        }
      } else{
        stop('TBD')
        # FIXME: Also make checks for correct formula, e.g. if variable is contained within object
        }
      model[['equation']][[ty]] <- f
      rm(f)
    }

    # Engine specific preparations

    # Sample nearest predictor values
    types <- names( x$biodiversity$get_types() )
    if('poipo' %in% types){
      model[['data']][['poipo_values']] <-
        get_ngbvalue(
          coords = x$biodiversity$get_coordinates('poipo'),
          env = x$predictors$get_data(df = TRUE, na.rm = FALSE),
          field_space = c('x','y'),
          longlat = raster::isLonLat(x$background)
        )
    }

    # Create INLA projection matrix from mesh to observed data from the points
    mat_proj <- INLA::inla.spde.make.A(
                        mesh = x$engine$get_data('mesh'),
                        loc = as.matrix(model[['data']][['poipo_values']][,c('x','y')])
                                  )
    #

    # Create INLA stack
    # The three main inla.stack() arguments are a vector list with the data (data),
    # a list of projector matrices (each related to one block effect,
    # A) and the list of effects (effects).

    # Response for inla stack
    resp <- x$biodiversity$get_columns_occ()[['poipo']]
    ll_resp <- list()
    ll_resp[[ resp ]] <- cbind( x$biodiversity$get_data('poipo')[,resp])
    # A column for the offset
    ll_resp[[ 'e' ]] <- rep(0, nrow(model$data$poipo_values) )

    # Effects matrix
    ll_effects <- list()
    # Note, order adding this is important apparently...
    ll_effects[['predictors']] <- model[['data']][['poipo_values']][,x$get_predictor_names()] # Get only the covariates to use
    ll_effects[['intercept']] <- list(intercept = c(1:x$engine$get_data('mesh')$n) )
    if(x$get_latent()=='<Spatial>'){
      spde <- x$engine$data$latentspatial
      ll_effects[['spatial.field']] <- list(Bnodes = 1:spde$n.spde)
      # Define projection matrix
      A = list(mat_proj, 1, 1)
    } else {
      spde <- NULL
      A = list(1, mat_proj)
    }

    # Create the stack for the fitted model
    stk_resp <-
      INLA::inla.stack(
        data =  ll_resp,             # Response
        A    =  A,   # Predictor projection matrix. 1 is included to make a list
        effects = ll_effects,        # Effects matrix
        tag = paste0('obs_','poipo') # Description tag
      )

    # Create and join prediction stack
    stk_pred <- inla_make_prediction_stack(stk_resp = stk_resp,
                                           cov = model$predictors,
                                           mesh = x$engine$get_data('mesh'),
                                           type = 'poipo',
                                           spde = spde)
    stk_full <- INLA::inla.stack(stk_resp, stk_pred)

    # Train the object
    fit <- INLA::inla(formula = model$equation$poipo, # The specified formula
                      data  = inla.stack.data(stk_full),  # The data stack
                      quantiles = c(0.05, 0.5, 0.95),
                      E = inla.stack.data(stk_full)$e,
                      family= 'poisson',   # Family the data comes from
                      control.family = list(link = "log"), # Control options
                      control.predictor=list(A = inla.stack.A(stk_full), compute = FALSE),  # Compute for marginals of the predictors
                      control.compute = list(cpo = TRUE, waic = TRUE, config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                      verbose = FALSE, # To see the log of the model runs
#                      control.inla(strategy = 'simplified.laplace', huge = TRUE), # To make it run faster...
                      num.threads = parallel::detectCores()-1
                      )


    index.pred <- INLA::inla.stack.index(stk_full, 'pred_poipo')$data
    post <- fit$summary.linear.predictor[index.pred, ]
    # TODO: Wrap Backtransform into a distribution dependent function
    post[,c('mean','0.05quant','0.5quant','0.95quant','mode')] <- exp(post[,c('mean','0.05quant','0.5quant','0.95quant','mode')])
    post <- subset(post, select = c('mean','sd','0.05quant','0.5quant','0.95quant','mode') )
    names(post) <- make.names(names(post))

    # Fill output rasters
    out <- fill_rasters(post = post,background = background)
    out <- raster::mask(out, background) # Mask with background
    plot(out$mean, main = 'Posterior prediction (mean lambda) using INLA')
    plot(as(virtual_points,'Spatial'),add =TRUE)

    #out <- x$engine$train()

    # TODO: Implement spatial crossvalidation, e.g. https://arxiv.org/pdf/2004.02324.pdf

    # return output object
    #return(out)
  }
)
