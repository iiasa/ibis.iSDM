#' @include bdproto-engine.R bdproto-distributionmodel.R
NULL
#' Use STAN as engine
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param chains A positive integer specifying the number of Markov chains. The default is 4.
#' @param iter A positive integer specifying the number of iterations for each chain (including warmup). The default is 2000.
#' @param warmup positive integer specifying the number of warmup (aka burnin) iterations per chain. If step-size adaptation is on (which it is by default), this also controls the number of iterations for which adaptation is run (and hence these warmup samples should not be used for inference). The number of warmup iterations should be smaller than iter and the default is iter/2.
#' @param cores If set to NULL take values from specified ibis option getOption('ibis.nthread')
#' @param init Initial values for parameters. Default: 'random'. Can also be specified as list (see: rstan::stan)
#' @param algorithm Mode used to sample from the posterior. See [`rstanarm::stan_glm`]
#' @param control See rstan::stan for more details on specifying the controls
#' @param ... Other variables
#' @seealso rstan
#' @name engine_stan
NULL
#' @rdname engine_stan
#' @export
engine_stan <- function(x,
                        chains = 4,
                        iter = 2000,
                        warmup = floor(iter/2),
                        init = "random",
                        cores = getOption("ibis.nthread"),
                        algorithm = 'sampling',
                        control = NULL,
                        ...) {
  # Check whether INLA package is available
  check_package('rstan')
  if(!isNamespaceLoaded("rstan")) { attachNamespace("rstan");requireNamespace('rstan') }
  # assert that arguments are valid
  assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                          inherits(x$background,'sf')
  )
  # Use rstanarm for now (latter impossible to specify spatial latent effects yet, but see https://github.com/stan-dev/rstanarm/issues/207)
  # https://github.com/ConnorDonegan/Stan-IAR
  # rstanarm object to improve speed through approximation
  check_package('rstanarm')
  if(!isNamespaceLoaded("rstanarm")) { attachNamespace("rstanarm");requireNamespace('rstanarm') }
  # assert that arguments are valid
  assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                          inherits(x$background,'sf')
  )

  # Other checks of parameters
  assertthat::assert_that(
    is.numeric(chains), is.numeric(iter), is.numeric(warmup),
    is.null(cores) || is.numeric(cores),
    is.character(init) || is.list(init),
    is.null(control) || is.list(control),
    is.character(algorithm),
    msg = 'Input parameters wrongly specified!'
  )
  # Match algorithm
  algorithm <- match.arg(algorithm, c("sampling", "optimizing", "meanfield", "fullrank"),several.ok = FALSE)

  if(is.null(cores)) cores <- getOption('ibis.nthread')

  # Create a background raster
  if(is.Waiver(x$predictors)){
    # Create from background
    template <- raster::raster(
      ext = raster::extent(x$background),
      crs = raster::projection(x$background),
      res = c(diff( (sf::st_bbox(x$background)[c(1,3)]) ) / 100, # Simplified assumption for resolution
              diff( (sf::st_bbox(x$background)[c(1,3)]) ) / 100
      )
    )
  } else {
    # If predictor existing, use them
    template <- emptyraster(x$predictors$get_data() )
  }
  # Burn in the background
  template <- raster::rasterize(x$background, template, field = 0)

  # Set engine in distribution object
  x$set_engine(
    bdproto(
      "STAN-Engine",
      Engine,
      name = "<STAN>",
      data = list(
        'template' = template
      ),
      # Stan options
      stan_param = list(
        chains = chains, iter = iter,
        warmup = warmup, init = init,
        cores = cores, algorithm = algorithm,
        control = control
      ),
      # Function to respecify the control parameters
      set_control = function(self,
                             chains = 4,
                             iter = 2000,
                             warmup = 500,
                             init = "random",
                             cores = NULL,
                             control = NULL){

        # Overwrite existing
        self$stan_param <- list(
          chains = chains, iter = iter,
          warmup = warmup, init = init,
          cores = cores, algorithm = algorithm,
          control = control
        )

      },
      # Spatial latent effect
      get_equation_latent_spatial = function(self){
        return( NULL )
      },
      # Setup a model
      setup = function(self, model, settings = NULL, ...){
        # Simple security checks
        assertthat::assert_that(
          assertthat::has_name(model, 'background'),
          assertthat::has_name(model, 'biodiversity'),
          inherits(settings,'Settings') || is.null(settings),
          nrow(model$predictors) == ncell(self$get_data('template')),
          length(model$biodiversity) == 1 # Only works with single likelihood. To be processed separately
        )
        # Detect and format the family
        fam <- model$biodiversity[[1]]$family
        self$data[['family']] <- fam

        # Check that all stan parameters are appropriately set
        assertthat::assert_that(
          is.numeric(self$stan_param$chains),
          is.numeric(self$stan_param$iter),
          is.numeric(self$stan_param$warmup)
        )
        # Set cores
        options(mc.cores = self$stan_param$cores)

        # Add pseudo-absence points if necessary
        # Include nearest predictor values for each
        if('poipo' == model$biodiversity[[1]]$type) {

          # Get background layer
          bg <- x$engine$get_data('template')
          assertthat::assert_that(!is.na(cellStats(bg,min)))

          # Sample pseudo absences
          abs <- create_pseudoabsence(
            env = model$predictors,
            presence = model$biodiversity[[1]]$observations,
            bias = settings$get('bias_variable'),
            template = bg,
            npoints = ifelse(ncell(bg)<10000,ncell(bg),10000), # FIXME: Ideally query this from settings
            replace = TRUE
          )
          abs$intercept <- 1 # Add dummy intercept
          # Combine absence and presence and save
          abs_observations <- abs[,c('x','y')]; abs_observations[['observed']] <- 0
          # Furthermore rasterize observed presences
          pres <- raster::rasterize(model$biodiversity[[1]]$predictors[,c('x','y')], bg,
                                    fun = 'count', background = 0)
          # If family is not poisson, assume factor distribution
          # FIXME: Ideally this is better organized through family
          if(model$biodiversity[[1]]$family != 'poisson') pres[] <- ifelse(pres[]==1,1,0)
          obs <- cbind( data.frame(observed = raster::extract(pres, model$biodiversity[[1]]$observations[,c('x','y')])),
                        model$biodiversity[[1]]$observations[,c('x','y')] )
          model$biodiversity[[1]]$observations <- rbind(obs, abs_observations)

          # Format out
          df <- rbind(model$biodiversity[[1]]$predictors,
                      abs[,c('x','y','intercept', model$biodiversity[[1]]$predictors_names)]) %>%
            subset(., complete.cases(.) )

          # Preprocessing security checks
          assertthat::assert_that( all( model$biodiversity[[1]]$observations[['observed']] >= 0 ),
                                   any(!is.na(rbind(obs, abs_observations)[['observed']] )),
                                   nrow(df) == nrow(model$biodiversity[[1]]$observations)
          )
          # Add offset if existent
          if(!is.Waiver(x$offset)) df[[x$get_offset()]] <- raster::extract(x$offset, df[,c('x','y')])

          # Define expectation as very small vector following Renner et al.
          w <- ppm_weights(df = df,
                           pa = model$biodiversity[[1]]$observations[['observed']],
                           bg = bg,
                           weight = 1e-6
          )
          df$w <- w # Also add as column

          model$biodiversity[[1]]$predictors <- df
          model$biodiversity[[1]]$expect <- w
        }

        # Messenger
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Engine setup.')

        # Return modified model object
        return(model)
      },
      train = function(self, model, settings, ...){
        # Messenger
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Starting fitting...')

        # Define an algorithm for MCMC sampling
        # an be "sampling" for MCMC (the default), "optimizing" for optimization,
        # "meanfield" for variational inference with independent normal distributions,
        # or "fullrank" for variational inference with a multivariate normal distribution.
        settings$set('algorithm', self$stan_param$algorithm)

        # Set a model seed for reproducibility
        settings$set('seed', 31337)

        # --- #
        # Get output raster
        prediction <- self$get_data('template')

        # Collect data
        equation <- model$biodiversity[[1]]$equation
        fam <- self$get_data('family')
        data <- cbind(model$biodiversity[[1]]$predictors,
                      data.frame(observed = model$biodiversity[[1]]$observations[,'observed']) )

        # Full data for prediction
        full <- model$predictors
        full <- subset(full, select = c('x','y',model$biodiversity[[1]]$predictors_names))
        full$cellid <- rownames(full) # Add rownames
        full[is.na(full)] <- 0 # For Prediction set all NA to 0. Mask later
        # full <- subset(full, complete.cases(full))

        # Make parameter list
        params <- list()

        # Weights
        if(model$biodiversity[[1]]$family=='poisson'){
          params$w <- model$biodiversity[[1]]$expect
        }

        # Priors
        params$prior <- rstanarm::hs() # Horseshoe priors
        params$prior_intercept <- rstanarm::normal(0, 2.5)
        params$prior_aux <- rstanarm::exponential(1) # Half-cauchy prior

        params$adapt_delta <- 0.95
        params$seed <- settings$get('seed')

        # Offsets
        if('offset' %in% names(model$biodiversity[[1]]) ){
          # Add offset to full prediction and load vector
          params$off <- model$biodiversity[[1]]$offset[, names(model$offset)[3] ]
          # Also add to full
        }

        # Model estimation
        params$algorithm <- settings$get('algorithm')
        if(settings$get('algorithm')=="sampling"){
           params$chains = self$stan_param$chains
           params$iter = self$stan_param$iter
           params$cores = self$stan_param$cores
           params$warmup = self$stan_param$warmup

        }
        # ---- #
        # Fitting
        fit_stan <- do.call( rstanarm::stan_glm, c(
          list(formula = equation, data = data, family = fam), params )
        )

        # Prediction
        if(!settings$get('inference_only')){
          # Messager
          if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Starting prediction...')

          # Set target variables to bias_value for prediction if specified
          if(!is.Waiver(settings$get('bias_variable'))){
            for(i in 1:length(settings$get('bias_variable'))){
              if(settings$get('bias_variable')[i] %notin% names(full)) next()
              full[[settings$get('bias_variable')[i]]] <- settings$get('bias_value')[i]
            }
          }
          params <- list()
          # Do the prediction by sampling from the posterior
          if('offset' %in% names(model$biodiversity[[1]]) ){
            params$off <- model$offset[, 3]
          }

          suppressWarnings(
            post_stan <- do.call( rstanarm::posterior_predict, c(
              list(object = fit_stan, newdata = full,
                   draws = ifelse(self$stan_param$iter >1000,1000,self$stan_param$iter)), params )
            )
          )
          # Summarize the posterior
          pred_stan <- data.frame(
            mean = colMeans(post_stan),
            sd = apply(post_stan, 2, sd),
            Q0.05 = apply(post_stan, 2, function(x) quantile(x, 0.05)),
            Q0.5 = apply(post_stan, 2, function(x) quantile(x, 0.5)),
            Q0.95 = apply(post_stan, 2, function(x) quantile(x, 0.95)),
            mode = apply(post_stan, 2, raster::modal)
          )
          # suppressWarnings(
          #   pred_stan <- rstanarm:::predict.stanreg(object = fit_stan,
          #                            newdata = full,se.fit = TRUE,
          #                            type = 'response')
          # )
          # Fill output
          prediction <- fill_rasters(post = pred_stan,background = prediction)
          prediction <- raster::mask(prediction, model$background) # Mask with background

        } else {
          prediction <- NULL
        }

        # Compute end of computation time
        settings$set('end.time', Sys.time())

        # Define output
        # Create output
        out <- bdproto(
          "STAN-Model",
          DistributionModel,
          id = new_id(),
          model = model,
          settings = settings,
          fits = list(
            "fit_best" = fit_stan,
            "fit_best_equation" = equation,
            "prediction" = prediction
          ),
          # Project function
          project = function(self, newdata){
            return(NULL)
          },
          # Partial effect
          partial = function(self, x.vars, constant = NULL, variable_length = 100, plot = FALSE){
            return(NULL)
          },
          # Spatial partial effect plots
          spartial = function(self, x.vars, constant = NULL, plot = TRUE,...){
            return(NULL)
          },
          # Spatial latent effect
          plot_spatial = function(self, plot = TRUE){
            return(NULL)
          }
        )
        # Return
        return(out)

      } # End of Train
    )
  ) # End of engine definition
}
