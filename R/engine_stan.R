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
#' @param algorithm Mode used to sample from the posterior. See [`cmdstanr`] package. Default: "sampling"
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
  stan_check_cmd(install = TRUE)

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
  algorithm <- match.arg(algorithm, c("sampling", "optimize", "variational"),several.ok = FALSE)

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
          nrow(model$predictors) == ncell(self$get_data('template'))
        )
        # Check that all stan parameters are appropriately set
        assertthat::assert_that(
          is.numeric(self$stan_param$chains),
          is.numeric(self$stan_param$iter),
          is.numeric(self$stan_param$warmup)
        )
        # Set cores
        options(mc.cores = self$stan_param$cores)

        # Stan procedure - First add integration points to all poipo datasets
        # FIXME: Possibly outsoure this across methods
        for(i in 1:length(model$biodiversity)){
          # Add pseudo-absence points if necessary
          # Include nearest predictor values for each
          if('poipo' == model$biodiversity[[i]]$type){

            # Get background layer
            bg <- x$engine$get_data('template')
            assertthat::assert_that(!is.na(cellStats(bg,min)))

            # Sample pseudo absences
            abs <- create_pseudoabsence(
              env = model$predictors,
              presence = model$biodiversity[[i]]$observations,
              bias = settings$get('bias_variable'),
              template = bg,
              npoints = ifelse(ncell(bg)<10000,ncell(bg),10000), # FIXME: Ideally query this from settings
              replace = TRUE
            )
            abs$intercept <- 1 # Add dummy intercept
            # Combine absence and presence and save
            abs_observations <- abs[,c('x','y')]; abs_observations[['observed']] <- 0
            # Furthermore rasterize observed presences
            pres <- raster::rasterize(model$biodiversity[[i]]$predictors[,c('x','y')], bg,
                                      fun = 'count', background = 0)
            # If family is not poisson, assume factor distribution
            # FIXME: Ideally this is better organized through family
            if(model$biodiversity[[i]]$family != 'poisson') pres[] <- ifelse(pres[]==1,1,0)
            obs <- cbind( data.frame(observed = raster::extract(pres, model$biodiversity[[i]]$observations[,c('x','y')])),
                          model$biodiversity[[i]]$observations[,c('x','y')] )
            model$biodiversity[[i]]$observations <- rbind(obs, abs_observations)

            # Format out
            df <- rbind(model$biodiversity[[i]]$predictors,
                        abs[,c('x','y','intercept', model$biodiversity[[i]]$predictors_names)]) %>%
              subset(., complete.cases(.) )

            # Preprocessing security checks
            assertthat::assert_that( all( model$biodiversity[[i]]$observations[['observed']] >= 0 ),
                                     any(!is.na(rbind(obs, abs_observations)[['observed']] )),
                                     nrow(df) == nrow(model$biodiversity[[i]]$observations)
            )
            # Add offset if existent
            if(!is.Waiver(x$offset)) df[[x$get_offset()]] <- raster::extract(x$offset, df[,c('x','y')])

            # Define expectation as very small vector following Renner et al.
            w <- ppm_weights(df = df,
                             pa = model$biodiversity[[i]]$observations[['observed']],
                             bg = bg,
                             weight = 1e-6
            )
            df$w <- w # Also add as column

            model$biodiversity[[i]]$predictors <- df
            model$biodiversity[[i]]$expect <- w
          }
        }
        # --- #
        # Messenger
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Building stan code.')
        sm_code <- vector("list", 7)
        names(sm_code) <- c("functions","data","transformed_data","parameters","transformed_parameters",
                            "model","generated_quantities")

        # Any spatial or other functions needed?
        if(!is.null(self$get_equation_latent_spatial())){
          # Stan functions for CAR and GP models
          ir <- readLines( system.file("src/spatial_functions.stan",package = "ibis.iSDM") )
          assertthat::assert_that(length(ir)>0)
          for(i in ir) sm_code$functions <- append(sm_code$functions, i)
        }

        # Load all the data parameters
        ir <- readLines( system.file("src/data_parameters.stan",package = "ibis.iSDM") )
        assertthat::assert_that(length(ir)>0)
        for(i in ir) sm_code$data <- append(sm_code$data, i)

        # Transformed data

        # Transformed parameters

        # Add (default) priors to model likelihood if set
        if(!is.Waiver(model$priors)){
          # Parameters
          sm_code$parameters <- append(sm_code$parameters, "vector[K] beta;")

        } else {
          # Add regularized horseshoe prior
          # See brms::horseshoe
          ir <- readLines( system.file("src/prior_functions.stan",package = "ibis.iSDM") )
          assertthat::assert_that(length(ir)>0)
          for(i in ir) sm_code$functions <- append(sm_code$functions, i)

          sm_code$data <- append(sm_code$data,"
                                  // data for the horseshoe prior
                                  real<lower=0> hs_df;  // local degrees of freedom
                                  real<lower=0> hs_df_global;  // global degrees of freedom
                                  real<lower=0> hs_df_slab;  // slab degrees of freedom
                                  real<lower=0> hs_scale_global;  // global prior scale
                                  real<lower=0> hs_scale_slab;  // slab prior scale"
          )

          # Parameters and transformed parameters
          sm_code$parameters <- append(sm_code$parameters,"
          // local parameters for horseshoe prior
          vector[K] zb;
          vector<lower=0>[K] hs_local;
          // horseshoe shrinkage parameters
          real<lower=0> hs_global;  // global shrinkage parameters
          real<lower=0> hs_slab;  // slab regularization parameter
                                       ")
          sm_code$transformed_parameters <- append(sm_code$transformed_parameters,"
            vector[K] beta;  // population-level effects
            // compute actual regression coefficients
            beta = horseshoe(zb, hs_local, hs_global, hs_scale_slab^2 * hs_slab);
          ")

          # Finally add priors to model
          sm_code$model <- append(sm_code$model, "
          // priors including constants
          target += std_normal_lpdf(zb);
          target += student_t_lpdf(hs_local | hs_df, 0, 1)
          - rows(hs_local) * log(0.5);
          target += student_t_lpdf(hs_global | hs_df_global, 0, hs_scale_global)
          - 1 * log(0.5);
          target += inv_gamma_lpdf(hs_slab | 0.5 * hs_df_slab, 0.5 * hs_df_slab);
          ")

          # sm_code$model <- append(sm_code$model,
          # "for (j in 1:K){
          #   target += normal_lpdf(b[j] | 0, 5);
          # }"
          # )
        }

        # Now add the model depending on the type
        if(length(model$biodiversity)>1){
          # For integrated model
          stop("TBD")

        } else if(model$biodiversity[[1]]$type == "poipo" && model$biodiversity[[1]]$family == "poisson"){
          # For poisson process model add likelihood
          ir <- readLines( system.file("src/poipo_ll_poisson.stan",package = "ibis.iSDM") )
          assertthat::assert_that(length(ir)>0)
          for(i in ir) sm_code$model <- append(sm_code$model, i)

        } else if(model$biodiversity[[1]]$type == "poipa" && model$biodiversity[[1]]$family == "binomial"){
          # For logistic regression
          ir <- readLines( system.file("src/poipa_ll_bernoulli.stan",package = "ibis.iSDM") )
          assertthat::assert_that(length(ir)>0)
          for(i in ir) sm_code$model <- append(sm_code$model, i)
        } else {
          # Else
          stop("Model as of now not implemented for Stan!")
        }

        # Wrap list entries in model code and save in model object
        self$set_data("stancode", wrap_stanmodel(sm_code))

        # --- #
        # Messenger
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Engine setup.')

        # Return modified model object
        return(model)
      },
      train = function(self, model, settings, ...){
        # Messenger
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Starting fitting...')

        # Define an algorithm for MCMC sampling
        # an be "sampling" for MCMC (the default), "optimize" for optimization,
        # "variational" for variational inference with independent normal distributions,
        settings$set('algorithm', self$stan_param$algorithm)

        # Set a model seed for reproducibility
        # FIXME: Ideally this is passed on better and earlier
        settings$set('seed', 31337)

        # --- #
        # Collect data for stan modelling
        if(length(model$biodiversity)>1){
          stop("done")
        } else {
          # Format datalist
          dl <- list(
            N = nrow( model$biodiversity[[1]]$observations),
            observed = model$biodiversity[[1]]$observations[["observed"]],
            X = as.matrix( model$biodiversity[[1]]$predictors[, model$biodiversity[[1]]$predictors_names] ),
            K = length( model$biodiversity[[1]]$predictors_names ),
            offset_exposure = log(model$biodiversity[[1]]$expect),
            has_intercept = attr(terms(model$biodiversity[[1]]$equation), "intercept"),
            has_spatial = ifelse(is.null(self$get_equation_latent_spatial()), 0, 1),
            # Horseshoe prior default parameters
            # FIXME: Allow passing this one via a parameter
            hs_df = 1,
            hs_df_global = 1, hs_df_slab = 4,
            hs_scale_global = 1, hs_scale_slab = 2
          )
        }

        # Get output raster and new data
        prediction <- self$get_data('template')
        # Full data for prediction
        full <- model$predictors
        full <- subset(full, select = c('x','y',model$predictors_names))
        suppressWarnings(
          full <- sp::SpatialPointsDataFrame(coords = full[,c("x","y")],
                                              data = full,
                                              proj4string = sp::CRS( sp::proj4string(as(model$background, "Spatial")) )
          )
        )
        full <- as(full, 'SpatialPixelsDataFrame')
        # Remove missing data
        # full <- subset(full, complete.cases(full@data))
        # full$cellid <- rownames(full) # Add rownames
        # full[is.na(full)] <- 0 # For Prediction set all NA to 0. Mask later

        # Model estimation
        # ---- #
        # Fitting
        fpath_code <- write_stanmodel( self$get_data("stancode") )
        fit_stan <- run_stan(
          model_code = fpath_code,
          data = dl,
          algorithm = settings$get('algorithm'),
          cores = self$stan_param$cores,
          chains = self$stan_param$chains,
          iter = self$stan_param$iter,
          warmup = self$stan_param$warmup,
          path = getwd(),
          force = TRUE # Force recompile
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
          # Add offset
          if('offset' %in% names(model$biodiversity[[1]]) ){
            params$off <- model$offset[, 3]
          }

          # Do the prediction by sampling from the posterior
          pred_stan <- posterior_predict_stanfit(obj = fit_stan,
                                         form = to_formula(paste0("observed ~ ", paste(model$biodiversity[[1]]$predictors_names,collapse = " + "))),
                                         newdata = full@data,
                                         offset = NULL
          )

          # Convert full to raster
          prediction <- raster::stack(full)
          # Fill output
          prediction <- fill_rasters(post = pred_stan, background = prediction)
          prediction <- raster::mask(prediction, model$background) # Mask with background
          # plot(prediction, col = ibis.iSDM:::ibis_colours$sdm_colour)

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
