#' @include bdproto-engine.R utils-spatial.R bdproto-distributionmodel.R
NULL

#' Engine for Bayesian regularized regression models
#'
#' @description
#' Efficient MCMC algorithm for linear regression models that makes use of
#' 'spike-and-slab' priors for some modest regularization on the amount of posterior
#' probability for a subset of the coefficients.
#' @details TBD
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param iter [`numeric`] on the number of MCMC iterations to run (Default: \code{10000}).
#' @param nthreads [`numeric`] on the number of CPU-threads to use for data augmentation.
#' @param type The mode used for creating posterior predictions. Either making \code{"link"} or \code{"response"} (Default: \code{"response"}).
#' @param ... Other none specified parameters passed on to the model.
#' @references
#' * Nguyen, K., Le, T., Nguyen, V., Nguyen, T., & Phung, D. (2016, November). Multiple kernel learning with data augmentation. In Asian Conference on Machine Learning (pp. 49-64). PMLR.
#' * Steven L. Scott (2021). BoomSpikeSlab: MCMC for Spike and Slab Regression. R package version 1.2.4. https://CRAN.R-project.org/package=BoomSpikeSlab
#' @family engine
#' @name engine_breg
NULL
#' @rdname engine_breg
#' @export

engine_breg <- function(x,
                           iter = 10000,
                           nthread = getOption('ibis.nthread'),
                           type = "response",
                           ...) {

  # Check whether xgboost package is available
  check_package('BoomSpikeSlab')
  if(!("BoomSpikeSlab" %in% loadedNamespaces()) || ('BoomSpikeSlab' %notin% sessionInfo()$otherPkgs) ) {
    try({requireNamespace('BoomSpikeSlab');attachNamespace("BoomSpikeSlab")},silent = TRUE)
  }

  # assert that arguments are valid
  assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                          inherits(x$background,'sf'),
                          is.numeric(iter),
                          is.character(type),
                          is.numeric(nthread)
  )
  type <- match.arg(type, choices = c("link", "response"),several.ok = FALSE)

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

  # Set up the parameter list
  params <- list(
    iter = iter,
    nthread = nthread,
    type = type,
    ...
  )

  # Print a message in case there is already an engine object
  if(!is.Waiver(x$engine)) myLog('[Setup]','yellow','Replacing currently selected engine.')

  # Set engine in distribution object
  x$set_engine(
    bdproto(
      "BREG-Engine",
      Engine,
      name = "<BREG>",
      data = list(
        'template' = template,
        'params' = params
      ),
      # Dummy function for spatial latent effects
      calc_latent_spatial = function(self, type = NULL, priors = NULL){
        new_waiver()
      },
      # Dummy function for getting the equation of latent effects
      get_equation_latent_spatial = function(self, method){
        new_waiver()
      },
      # Function to respecify the control parameters
      set_control = function(self,
                             params
      ){
        assertthat::assert_that(is.list(params))
        # Overwrite existing
        self$data$params <- params
        invisible()
      },
      # Setup function
      setup = function(self, model, settings = NULL, ...){
        # Simple security checks
        assertthat::assert_that(
          assertthat::has_name(model, 'background'),
          assertthat::has_name(model, 'biodiversity'),
          inherits(settings,'Settings') || is.null(settings),
          nrow(model$predictors) == ncell(self$get_data('template')),
          !is.Waiver(self$get_data("params")),
          length(model$biodiversity) == 1 # Only works with single likelihood. To be processed separately
        )
        # Messenger
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Engine setup.')

        # Get parameters
        params <- self$data$params

        # Distribution specific procedure
        fam <- model$biodiversity[[1]]$family

        # If a poisson family is used, weight the observations by their exposure
        if(fam == "poisson"){
          # Get background layer
          bg <- self$get_data("template")
          assertthat::assert_that(!is.na(cellStats(bg,min)))

          abs <- create_pseudoabsence(
            env = model$predictors,
            presence = model$biodiversity[[1]]$observations,
            bias = settings$get('bias_variable'),
            template = bg,
            npoints = ifelse(ncell(bg)<10000,ncell(bg),10000),
            replace = TRUE
          )
          # Combine absence and presence and save
          abs$intercept <- 1
          abs_observations <- abs[,c('x','y')]; abs_observations[['observed']] <- 0

          # Rasterize observed presences
          pres <- raster::rasterize(model$biodiversity[[1]]$predictors[,c('x','y')],
                                    bg, fun = 'count', background = 0)

          # Get cell ids
          ce <- raster::cellFromXY(pres, model[['biodiversity']][[1]]$observations[,c('x','y')])

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
          model$biodiversity[[1]]$observations <- rbind(obs, abs_observations)

          # Format out
          df <- rbind(envs,
                      abs[,c('x','y','intercept', model$biodiversity[[1]]$predictors_names)]) %>%
            subset(., complete.cases(.) )

          # Pre-processing security checks
          assertthat::assert_that( all( model$biodiversity[[1]]$observations[['observed']] >= 0 ),
                                   any(!is.na(rbind(obs, abs_observations)[['observed']] )),
                                   nrow(df) == nrow(model$biodiversity[[1]]$observations)
          )

          # Define expectation as very small vector following Renner et al.
          w <- ppm_weights(df = df,
                           pa = model$biodiversity[[1]]$observations[['observed']],
                           bg = bg,
                           weight = 1e-6
          )
          assertthat::assert_that(length(w) == nrow(df))

          model$biodiversity[[1]]$predictors <- df
          model$biodiversity[[1]]$expect <- w

          # Get for the full dataset
          w_full <- ppm_weights(df = model$predictors,
                                pa = pres[],
                                bg = bg,
                                weight = 1 # Set those to 1 so that absences become ratio of pres/abs
          )

          # Add exposure to full model predictor
          model$exposure <- w_full

        } else if(fam == "binomial"){
          # calculating the case weights (equal weights)
          # the order of weights should be the same as presences and backgrounds in the training data
          prNum <- as.numeric(table(model$biodiversity[[1]]$observations[['observed']])["1"]) # number of presences
          bgNum <- as.numeric(table(model$biodiversity[[1]]$observations[['observed']])["0"]) # number of backgrounds
          w <- ifelse(model$biodiversity[[1]]$observations[['observed']] == 1, 1, prNum / bgNum)
          model$biodiversity[[1]]$expect <- w
          # Convert to numeric
          model$biodiversity[[1]]$observations$observed <- as.numeric( model$biodiversity[[1]]$observations$observed )
        }

        # Instead of invisible return the model object
        return( model )
      },
      # Training function
      train = function(self, model, settings, ...){
        assertthat::assert_that(
          inherits(settings,'Settings'),
          is.list(model),length(model)>1,
          # Check that model id and setting id are identical
          settings$modelid == model$id
        )
        # Messager
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Starting fitting...')

        # Verbosity
        verbose <- settings$get("verbose")

        # seed
        seed <- settings$get("seed")
        if(is.Waiver(seed)) seed <- 1337

        # Get output raster
        prediction <- self$get_data('template')

        # Get parameters control
        params <- self$get_data('params')

        # Check only linear and reset to linear booster then
        if(settings$data$only_linear) params$booster <- "gblinear" else params$booster <- "gbtree"

        # All other needed data for model fitting
        fam <- model$biodiversity[[1]]$family
        form <- model$biodiversity[[1]]$equation
        df <- cbind(model$biodiversity[[1]]$predictors,
                    data.frame(observed = model$biodiversity[[1]]$observations[,'observed'])
                    )
        w <- model$biodiversity[[1]]$expect # The expected exposure
        # Get full prediction container
        full <- model$predictors
        w_full <- model$exposure

        # Priors
        if(!is.Waiver(model$priors)){
          print("TBD")
        } else { pp <- NULL }

        assertthat::assert_that(
          is.null(w) || length(w) == nrow(df)
        )

        # --- #
        # Pass this parameter possibly on from upper levels
        if(settings$get('varsel') == "reg"){
          if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Starting hyperparameters search.')

          # Create combinations of random hyper parameters
          set.seed(20)
          if(params$booster == 'gblinear'){
            parameters_df <- expand.grid(
              lambda = seq(0,7,0.25), alpha = seq(0,1, 0.1),
              eval = NA
            )
          } else {
            parameters_df <- expand.grid(
              # overfitting
              max_depth = 1:6,
              gamma = 0:5,
              min_child_weight = 0:6,
              # Randomness
              subsample = runif(1, .7, 1),
              colsample_bytree = runif(1, .6, 1),
              eval = NA
            )
          }

          # Progressbar
          pb <- progress::progress_bar$new(total = nrow(parameters_df))
          # TODO: Could be parallized
          for (row in 1:nrow(parameters_df)){
            test_params <- list(
              booster = params$booster,
              objective = params$objective
            )
            if(test_params$booster=='gblinear'){
              test_params$lambda <- parameters_df$lambda[row]
              test_params$alpha <- parameters_df$alpha[row]
            } else {
              test_params$gamma <- parameters_df$gamma[row]
              test_params$max_depth <- parameters_df$max_depth[row]
              test_params$min_child_weight <- parameters_df$min_child_weight[row]
              test_params$subsample <- parameters_df$subsample[row]
              test_params$colsample_bytree <- parameters_df$colsample_bytree[row]
            }

            suppressMessages(
              test_xgb <- xgboost::xgboost(
                params = test_params,
                data = df_train,
                nrounds = 100,
                verbose = ifelse(verbose, 1, 0)
              )
            )
            if(verbose) pb$tick()
            parameters_df$eval[row] <- min(test_xgb$evaluation_log[,2])
          }
          # Get the one with minimum error and replace params values
          p <- parameters_df[which.min(parameters_df$eval),]
          p$eval <- NULL
          for(i in names(p)){ params[[i]] <- as.numeric(p[i]) }

          # Find the optimal number of rounds using 5-fold cross validation
          if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Crossvalidation for determining early stopping rule.')

          fit_cv <- xgboost::xgb.cv(
            params = params,
            data = df_train,
            verbose = ifelse(verbose, 1, 0),
            print_every_n = 100,
            nrounds = nrounds, nfold = 5,
            showsd = TRUE,      # standard deviation of loss across folds
            stratified = TRUE,  # sample is unbalanced; use stratified samplin
            maximize = FALSE,
            early_stopping_rounds = 10
          )
          # Set new number of rounds
          nround <- fit_cv$best_iteration
        }

        # Fit the model depending on the family
        if(fam == "poisson"){
          # Fitting poisson model
          fit_breg <- BoomSpikeSlab::poisson.spike(
            formula = form,
            exposure = w,
            niter = params$iter,
            data = df,
            prior = pp,
            nthreads = params$nthread,
            ping = ifelse( settings$get("verbose"), params$iter / 10 , 0),
            seed = seed
          )
        } else if(fam == "binomial"){
          fit_breg <- BoomSpikeSlab::logit.spike(
            formula = form,
            niter = params$iter,
            data = df,
            prior = pp,
            nthreads = params$nthread,
            ping = ifelse( settings$get("verbose"), params$iter / 10 , 0),
            seed = seed
          )
        } else {
          if(getOption('ibis.setupmessages')) myLog('[Estimation]','yellow','Non supported family: ', fam)
          fit_breg <- BoomSpikeSlab::lm.spike(
            formula = form,
            niter = params$iter,
            data = df,
            prior = pp,
            nthreads = params$nthread,
            ping = ifelse( settings$get("verbose"), params$iter / 10 , 0),
            seed = seed
          )
        }
        # --- #

        # Predict spatially
        if(!settings$get('inference_only')){
          # Messager
          if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Starting prediction...')

          # Make a prediction
          if(fam == "poisson"){
            suppressWarnings(
              pred_breg <- BoomSpikeSlab::predict.poisson.spike(
                object = fit_breg,
                newdata = full,
                exposure = w_full,
                burn = ceiling(params$iter*0.1),
                type = params$type,
                mean.only = FALSE # Return full posterior
              )
            )
          } else if(fam == "binomial"){
            suppressWarnings(
              pred_breg <- BoomSpikeSlab::predict.logit.spike(
                object = fit_breg,
                newdata = full,
                burn = ceiling(params$iter*0.1),
                type = params$type,
                mean.only = FALSE # Return full posterior
              )
            )
          } else {
            suppressWarnings(
              pred_breg <- BoomSpikeSlab::predict.lm.spike(
                object = fit_breg,
                newdata = full,
                burn = ceiling(params$iter*0.1),
                type = params$type,
                mean.only = FALSE # Return full posterior
              )
            )
          }
          # Summarize the posterior
          preds <- cbind(
            matrixStats::rowMeans2(pred_breg, na.rm = TRUE),
            matrixStats::rowSds(pred_breg, na.rm = TRUE),
            matrixStats::rowQuantiles(pred_breg, probs = c(.05,.5,.95), na.rm = TRUE),
            apply(pred_breg, 1, mode)
          ) %>% as.data.frame()
          names(preds) <- c("mean", "sd", "q05", "q50", "q95", "mode")
          preds$cv <- preds$mean / preds$sd

          # Fill output with summaries of the posterior
          prediction <- fill_rasters(preds, prediction)
          prediction <- raster::mask(prediction, self$get_data("template"))

        } else {
          # No prediction done
          prediction <- NULL
        }
        # Compute end of computation time
        settings$set('end.time', Sys.time())

        # Definition of BREG Model object ----
        # Create output
        out <- bdproto(
          "BREG-Model",
          DistributionModel,
          id = new_id(),
          model = model,
          settings = settings,
          fits = list(
            "fit_best" = fit_breg,
            "prediction" = prediction
          ),
          # Partial effects
          partial = function(self, x.var = NULL, plot = TRUE, ...){
            assertthat::assert_that(is.character(x.var) || is.null(x.var))
            mod <- self$get_data('fit_best')
            df <- self$model$biodiversity[[length( self$model$biodiversity )]]$predictors
            df <- subset(df, select = mod$feature_names)

            # Match x.var to argument
            if(is.null(x.var)){
              x.var <- colnames(df)
            } else {
              x.var <- match.arg(x.var, mod$feature_names, several.ok = FALSE)
            }

            # Check that variables are in
            assertthat::assert_that(all( x.var %in% colnames(df) ),
                                    msg = 'Variable not in predicted model.')

            pp <- data.frame()
            pb <- progress::progress_bar$new(total = length(x.var))
            for(v in x.var){
              p1 <- pdp::partial(mod, pred.var = v, ice = FALSE, center = TRUE,
                                 plot = FALSE, rug = TRUE, train = df)
              names(p1) <- c("partial", "yhat")
              p1$variable <- v
              pp <- rbind(pp, p1)
              if(length(x.var) > 1) pb$tick()
            }

            if(plot){
              # Make a plot
              ggplot2::ggplot(data = pp, ggplot2::aes(x = partial, y = yhat)) +
                ggplot2::theme_classic(base_size = 18) +
                ggplot2::geom_line() +
                ggplot2::labs(x = "", y = expression(hat(y))) +
                ggplot2::facet_wrap(~variable,scales = 'free')
            }
            # Return the data
            return(pp)
          },
          # Spatial partial dependence plot
          spartial = function(self, x.var, constant = NULL){
            assertthat::assert_that(is.character(x.var) || is.null(x.var),
                                    "model" %in% names(self))
            stop("TBD")
            # Get data
            mod <- self$get_data('fit_best')
            model <- self$model
            x.var <- match.arg(x.var, colnames(df), several.ok = FALSE)

            # Get predictor
            df <- subset(model$predictors, select = mod$feature_names)
            # Convert all non x.vars to the mean
            # Make template of target variable(s)
            template <- raster::rasterFromXYZ(cbind(model$predictors$x,model$predictors$y),
                                              crs = raster::projection(model$background))

            # Set all variables other the target variable to constant
            if(is.null(constant)){
              # Calculate mean
              # FIXME: for factor use mode!
              constant <- apply(df, 2, function(x) mean(x, na.rm=T))
              for(v in mod$feature_names[ mod$feature_names %notin% x.var]){
                if(v %notin% names(df) ) next()
                df[!is.na(df[v]),v] <- as.numeric( constant[v] )
              }
            } else {
              df[!is.na(df[,x.var]), mod$feature_names[ mod$feature_names %notin% x.var]] <- constant
            }
            df <- xgboost::xgb.DMatrix(data = as.matrix(df))

            # Spartial prediction
            suppressWarnings(
              pp <- xgboost:::predict.xgb.Booster(
                object = mod,
                newdata = df
              )
            )
            # Fill output with summaries of the posterior
            template[] <- pp
            names(template) <- 'mean'
            template <- raster::mask(template, model$background)

            # Quick plot
            raster::plot(template, col = ibis_colours$viridis_plasma, main = paste0(x.var, collapse ='|'))
            # Also return spatial
            return(template)
          },
          # Engine-specific projection function
          project = function(self, newdata){
            assertthat::assert_that(!missing(newdata),
                                    is.data.frame(newdata) || inherits(newdata, "xgb.DMatrix") )
            stop("TBD")

            mod <- self$get_data('fit_best')
            if(!inherits(newdata, "xgb.DMatrix")){
              assertthat::assert_that(
                all( mod$feature_names %in% colnames(newdata) )
              )
              newdata <- subset(newdata, select = mod$feature_names)
              newdata <- xgboost::xgb.DMatrix(as.matrix(newdata))
            }

            # Make a prediction
            suppressWarnings(
              pred_xgb <- xgboost:::predict.xgb.Booster(
                object = mod,
                newdata = newdata
              )
            )

            # Fill output with summaries of the posterior
            prediction <- emptyraster( self$get_data('prediction') ) # Background
            prediction[] <- pred_xgb
            prediction <- raster::mask(prediction, self$get_data('prediction') )

            return(prediction)
          }
        )
        return(out)
      }
    )
  ) # End of bdproto object
} # End of function
