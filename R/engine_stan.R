#' @include bdproto-engine.R bdproto-distributionmodel.R
NULL
#' Use Stan as engine
#'
#' @description Stan is probabilistic programming language that can be used to
#' specify most types of statistical linear and non-linear regression models.
#' Stan provides full Bayesian inference for continuous-variable models through Markov chain Monte Carlo methods
#' such as the No-U-Turn sampler, an adaptive form of Hamiltonian Monte Carlo sampling.
#' Stan code has to be written separately and this function acts as compiler to
#' build the stan-model.
#' **Requires the [cmdstanr] package to be installed!**
#' @details
#' By default the posterior is obtained through sampling, however stan also supports
#' approximate inference forms through penalized maximum likelihood estimation (see Carpenter et al. 2017).
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param chains A positive [`integer`] specifying the number of Markov chains (Default: \code{4} chains).
#' @param iter A positive [`integer`] specifying the number of iterations for each chain (including warmup). (Default: \code{2000}).
#' @param warmup A positive [`integer`] specifying the number of warmup (aka burnin) iterations per chain.
#' If step-size adaptation is on (Default: \code{TRUE}), this also controls the number of iterations for which
#' adaptation is run (and hence these warmup samples should not be used for inference).
#' The number of warmup iterations should be smaller than \code{iter} and the default is \code{iter/2}.
#' @param cores If set to NULL take values from specified ibis option \code{getOption('ibis.nthread')}.
#' @param init Initial values for parameters (Default: \code{'random'}). Can also be specified as [list] (see: [`rstan::stan`])
#' @param algorithm Mode used to sample from the posterior. Available options are \code{"sampling"}, \code{"optimize"},
#' or \code{"variational"}.
#' See [`cmdstanr`] package for more details. (Default: \code{"sampling"}).
#' @param control See [`rstan::stan`] for more details on specifying the controls.
#' @param type The mode used for creating posterior predictions. Either summarizing the linear \code{"predictor"} or \code{"response"} (Default: \code{"response"}).
#' @param ... Other variables
#' @seealso [rstan], [cmdstanr]
#' @note
#' The function \code{obj$stancode()} can be used to print out the stancode of the model.
#' @references
#' * Jonah Gabry and Rok Češnovar (2021). cmdstanr: R Interface to 'CmdStan'. https://mc-stan.org/cmdstanr, https://discourse.mc-stan.org.
#' * Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B., Betancourt, M., ... & Riddell, A. (2017). Stan: A probabilistic programming language. Journal of statistical software, 76(1), 1-32.
#' * Piironen, J., & Vehtari, A. (2017). Sparsity information and regularization in the horseshoe and other shrinkage priors. Electronic Journal of Statistics, 11(2), 5018-5051.
#' @family engine
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
                        type = "response",
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
  # Match algorithm and posterior prediction type
  algorithm <- match.arg(algorithm, c("sampling", "optimize", "variational"),several.ok = FALSE)
  type <- match.arg(type, c("response", "predictor"), several.ok = FALSE)
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
        control = control,
        type = type
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
            # Get cell ids
            ce <- raster::cellFromXY(pres, model[['biodiversity']][[i]]$observations[,c('x','y')])

            # If family is not poisson, assume factor distribution
            # FIXME: Ideally this is better organized through family
            if(model$biodiversity[[i]]$family != 'poisson') pres[] <- ifelse(pres[]==1,1,0)
            # Get new presence data
            obs <- cbind(
              data.frame(observed = raster::values(pres)[ce],
                         raster::xyFromCell(pres, ce) # Center of cell
              )
            ) |> unique() # Unique to remove any duplicate values (otherwise double counted cells)

            # Re-extract counts environment variables
            envs <- get_ngbvalue(coords = obs[,c('x','y')],
                                 env =  model$predictors[,c("x","y", model[['predictors_names']])],
                                 longlat = raster::isLonLat(self$get_data("template")),
                                 field_space = c('x','y')
            )
            envs$intercept <- 1
            # Overwrite observations
            model$biodiversity[[i]]$observations <- rbind(obs, abs_observations)

            # Format out
            df <- rbind(envs,
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

        # Add (gaussian) priors to model likelihood if set
        if((!is.Waiver(model$priors) || settings$get(what='varsel') == "none")){
          # Parameters
          sm_code$parameters <- append(sm_code$parameters, "vector[K] beta;")

          # Add priors for each variable for which it is set to the model
          sm_code$model <- append(sm_code$model, "// priors including constants")
          # Now add for each one a normal effect
          for(i in 1:length(model$predictors_names)){
            if(!is.Waiver(model$priors)){
              if(model$predictors_names[i] %in% model$priors$varnames()) {
                sm_code$model <- append(sm_code$model, paste0(
                  "target += normal_lpdf(beta[",i,"] | ",model$priors$get(model$predictors_names[8])[1],", ",model$priors$get(model$predictors_names[8])[2],");"
                ))
              } else {
                # Default gaussian prior
                sm_code$model <- append(sm_code$model, paste0("target += normal_lpdf(beta[",i,"] | 0, 2);"))
              }
            } else {
              # Default gaussian prior
              sm_code$model <- append(sm_code$model, paste0("target += normal_lpdf(beta[",i,"] | 0, 2);"))
            }
          }
        } else
          if( settings$get(what='varsel') == "reg" ){
          if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Adding regularized Bayesian priors.')
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
        settings$set('cores', self$stan_param$cores)
        settings$set('chains', self$stan_param$chains)
        settings$set('iter', self$stan_param$iter)
        settings$set('warmup', self$stan_param$warmup)
        settings$set('type', self$stan_param$type)
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
            offsets = log(model$biodiversity[[1]]$expect),
            has_intercept = attr(terms(model$biodiversity[[1]]$equation), "intercept"),
            has_spatial = ifelse(is.null(self$get_equation_latent_spatial()), 0, 1),
            # Horseshoe prior default parameters
            # FIXME: Allow passing this one via a parameter
            hs_df = 1,
            hs_df_global = 1, hs_df_slab = 4,
            hs_scale_global = 1, hs_scale_slab = 2
          )
          # If any additional offset is set, simply to the existing one in sum
          if(!is.Waiver(model$offset)) dl$offsets <- dl$offsets + log( model$biodiversity[[1]]$offset[,3] )
        }

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

          # Prepare prediction dataset
          prediction <- self$get_data('template') # Get output raster and new data
          # Full data for prediction
          full <- subset(model$predictors, select = c('x','y',model$predictors_names))

          # If poipo, add w to prediction container
          bd_poipo <- sapply(model$biodiversity, function(x) x$type) == "poipo"
          if(any(bd_poipo)){
            # FIXME: Bit hackish. See if works for other projections
            full$w <- unique(model$biodiversity[[which(bd_poipo)]]$expect)[2] # Absence location being second unique value
          }
          # Add offset if set
          if(!is.Waiver(model$offset)) {
            # Offsets are simply added linearly (albeit transformed)
            if(hasName(full,"w")) full$w <- full$w + model$offset[,3] else full$w <- model$offset[,3]
            # full[[colnames(model$offset)[3]]] <- model$offset[,3]
          }
          suppressWarnings(
            full <- sp::SpatialPointsDataFrame(coords = full[,c("x","y")],
                                               data = full,
                                               proj4string = sp::CRS( sp::proj4string(as(model$background, "Spatial")) )
            )
          )
          full <- as(full, 'SpatialPixelsDataFrame')

          # Set target variables to bias_value for prediction if specified
          if(!is.Waiver(settings$get('bias_variable'))){
            for(i in 1:length(settings$get('bias_variable'))){
              if(settings$get('bias_variable')[i] %notin% names(full)) next()
              full[[settings$get('bias_variable')[i]]] <- settings$get('bias_value')[i]
            }
          }

          # For Integrated model, follow poisson
          fam <- ifelse(length(model$biodiversity)>1, "poisson", model$biodiversity[[1]]$family)

          # Do the prediction by sampling from the posterior
          pred_stan <- posterior_predict_stanfit(obj = fit_stan,
                                         form = to_formula(paste0("observed ~ ", paste(model$biodiversity[[1]]$predictors_names,collapse = " + "))),
                                         newdata = full@data,
                                         offset = (full$w),
                                         family = fam,
                                         mode = self$stan_param$type # Simulated response
          )

          # Convert full to raster
          prediction <- raster::stack(full)
          # Fill output
          prediction <- fill_rasters(post = pred_stan, background = prediction)
          prediction <- raster::mask(prediction, model$background) # Mask with background
          # plot(prediction$mean, col = ibis.iSDM:::ibis_colours$sdm_colour)

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
            "prediction" = prediction,
            "sm_code" = self$get_data("sm_code")
          ),
          # Project function
          project = function(self, newdata, offset = NULL, type = NULL){
            assertthat::assert_that(
              nrow(newdata) > 0,
              all( c("x", "y") %in% names(newdata) ),
              is.null(offset) || is.numeric(offset),
              is.character(type) || is.null(type)
            )
            # Check that fitted model exists
            obj <- self$get_data("fit_best")
            model <- self$model
            if(is.null(type)) type <- self$settings$get("type")
            assertthat::assert_that(inherits(obj, "stanfit"),
                                    all(model$predictors_names %in% colnames(newdata)))

            # For Integrated model, follow poisson
            fam <- ifelse(length(model$biodiversity)>1, "poisson", model$biodiversity[[1]]$family)

            # Build prediction stack
            full <- subset(newdata, select = c('x','y', model$predictors_names))

            # If poipo, add w to prediction container
            bd_poipo <- sapply(model$biodiversity, function(x) x$type) == "poipo"
            if(any(bd_poipo)){
              # FIXME: Bit hackish. See if works for other projections
              full$w <- unique(model$biodiversity[[which(bd_poipo)]]$expect)[2] # Absence location being second unique value
            }
            # Add offset if set
            if(!is.null(offset)) {
              # Offsets are simply added linearly (albeit transformed)
              if(hasName(full,"w")) full$w <- full$w + offset else full$w <- offset
              # full[[colnames(model$offset)[3]]] <- model$offset[,3]
            }
            suppressWarnings(
              full <- sp::SpatialPointsDataFrame(coords = full[,c("x","y")],
                                                 data = full,
                                                 proj4string = sp::CRS( sp::proj4string(as(model$background, "Spatial")) )
              )
            )
            full <- as(full, 'SpatialPixelsDataFrame')

            # Do the prediction by sampling from the posterior
            pred_stan <- posterior_predict_stanfit(obj = obj,
                                                   form = to_formula(paste0("observed ~ ", paste(model$predictors_names,collapse = " + "))),
                                                   newdata = full@data,
                                                   offset = (full$w),
                                                   family = fam,
                                                   mode = type # Linear predictor
            )

            # Fill output with summaries of the posterior
            prediction <- emptyraster( raster::stack(full) ) # Background
            prediction <- fill_rasters(pred_stan, prediction)

            return(prediction)

          },
          # Partial effect
          partial = function(self, x.var, constant = NULL, length.out = 100, plot = FALSE, type = NULL){
            mod <- self$get_data('fit_best')
            model <- self$model
            if(is.null(type)) type <- self$settings$get("type")
            assertthat::assert_that(inherits(mod,'stanfit'),
                                    is.character(x.var),
                                    is.numeric(length.out),
                                    is.null(constant) || is.numeric(constant)
            )
            # Check that given variable is in x.var
            assertthat::assert_that(x.var %in% model$predictors_names)
            # Calculate
            rr <- sapply(model$predictors, function(x) range(x, na.rm = TRUE)) |> as.data.frame()
            df_partial <- list()
            # Add all others as constant
            if(is.null(constant)){
              for(n in names(rr)) df_partial[[n]] <- rep( mean(model$predictors[[n]], na.rm = TRUE), length.out )
            } else {
              for(n in names(rr)) df_partial[[n]] <- rep( constant, length.out )
            }
            df_partial[[x.var]] <- seq(rr[1,x.var], rr[2,x.var], length.out = length.out)
            df_partial <- df_partial %>% as.data.frame()

            # For Integrated model, follow poisson
            fam <- ifelse(length(model$biodiversity)>1, "poisson", model$biodiversity[[1]]$family)

            # Simulate from the posterior
            pred_part <- posterior_predict_stanfit(obj = mod,
                                                   form = to_formula(paste0("observed ~ ", paste(model$predictors_names,collapse = " + "))),
                                                   newdata = df_partial,
                                                   offset = NULL,
                                                   family = fam,
                                                   mode = type # Linear predictor
            )
            # Also attach the partial variable
            pred_part <- cbind("partial_effect" = df_partial[[x.var]], pred_part)
            if(plot){
              o <- pred_part
              pm <- ggplot2::ggplot(data = o, ggplot2::aes(x = partial_effect, y = mean,
                                                           ymin = mean-sd,
                                                           ymax = mean+sd) ) +
                ggplot2::theme_classic() +
                ggplot2::geom_ribbon(fill = "grey90") +
                ggplot2::geom_line() +
                ggplot2::labs(x = x.var, y = "Partial effect")
              print(pm)
            }
            return(pred_part) # Return the partial data
          },
          # Spatial partial effect plots
          spartial = function(self, x.var, constant = NULL, plot = TRUE,type = "predictor", ...){
            # Get model object and check that everything is in order
            mod <- self$get_data('fit_best')
            model <- self$model
            assertthat::assert_that(inherits(mod,'stanfit'),
                                    'model' %in% names(self),
                                    is.character(x.var),
                                    is.null(constant) || is.numeric(constant)
            )

            # Match variable name
            x.var <- match.arg(x.var, model$predictors_names, several.ok = FALSE)

            # Make spatial container for prediction
            suppressWarnings(
              df_partial <- sp::SpatialPointsDataFrame(coords = model$predictors[,c('x', 'y')],
                                                       data = model$predictors[, names(model$predictors) %notin% c('x','y')],
                                                       proj4string = sp::CRS( sp::proj4string(as(model$background, "Spatial")) )
              )
            )
            df_partial <- as(df_partial, 'SpatialPixelsDataFrame')

            # Add all others as constant
            if(is.null(constant)){
              for(n in names(df_partial)) if(n != x.var) df_partial[[n]] <- mean(model$predictors[[n]], na.rm = TRUE)
            } else {
              for(n in names(df_partial)) if(n != x.var) df_partial[[n]] <- constant
            }
            bd_poipo <- sapply(model$biodiversity, function(x) x$type) == "poipo"
            if(any(bd_poipo)){
              # FIXME: Bit hackish. See if works for other projections
              df_partial$w <- unique(model$biodiversity[[which(bd_poipo)]]$expect)[2] # Absence location being second unique value
            }

            # For Integrated model, follow poisson
            fam <- ifelse(length(model$biodiversity)>1, "poisson", model$biodiversity[[1]]$family)

            # Simulate from the posterior
            pred_part <- posterior_predict_stanfit(obj = mod,
                                                   form = to_formula(paste0("observed ~ ", paste(model$predictors_names,collapse = " + "))),
                                                   newdata = df_partial@data,
                                                   offset = NULL,
                                                   family = fam,
                                                   mode = type # Linear predictor
            )

            prediction <- emptyraster( self$get_data('prediction') ) # Background
            prediction <- fill_rasters(pred_part, prediction)

            # Do plot and return result
            if(plot){
              plot(prediction[[c("mean","sd")]], col = ibis_colours$viridis_orig)
            }
            return(prediction)
          },
          # Spatial latent effect
          plot_spatial = function(self, plot = TRUE){
            return(NULL)
          },
          # Custom function to show stan code
          show_code = function(self){
            message(
              self$get_data("sm_code")
            )
          }
        )
        # Return
        return(out)

      } # End of Train
    )
  ) # End of engine definition
}
