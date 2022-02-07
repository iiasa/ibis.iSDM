#' @include bdproto-engine.R utils-spatial.R bdproto-distributionmodel.R
NULL

#' Engine for extreme gradient boosting (XGBoost)
#'
#' @description Allows to estimate eXtreme gradient descent boosting for tree-based or linear boosting regressions.
#' The XGBoost engine is a flexible, yet powerful engine with many customization options,
#' supporting multiple options to perform single and multiclass regression
#' and classification tasks. For a full list of options users are advised to have a look at the
#' [xgboost::xgb.train] help file and [https://xgboost.readthedocs.io](https://xgboost.readthedocs.io).
#' @details
#' XGBoost supports the specification of monotonic contraints on certain variables. Within
#' ibis this is possible via [`XGBPrior`].
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param booster A [`character`] of the booster to use. Either "gbtree" or "gblinear" (Default: \code{gblinear})
#' @param learning_rate [`numeric`] value indicating the learning rate (eta).
#' Lower values generally being better but also computationally more costly.
#' @param nrounds [`numeric`] value giving the the maximum number of boosting iterations for cross-validation.
#' @param gamma [`numeric`] A regularization parameter in the model. Lower values for better estimates (Default: 3)
#' Also see [reg_lambda] parameter for the L2 regularization on the weights
#' @param reg_lambda [`numeric`] L2 regularization term on weights (Default: \code{0}).
#' @param reg_alpha [`numeric`] L1 regularization term on weights (Default: \code{0}).
#' @param max_depth [`numeric`] The Maximum depth of a tree (Default: \code{3}).
#' @param subsample [`numeric`] The ratio used for subsampling to prevent overfitting (Default: \code{0.8}).
#' @param ... Other none specificed parameters.
#' @seealso [xgboost::xgb.train]
#' @references
#' * Tianqi Chen and Carlos Guestrin, "XGBoost: A Scalable Tree Boosting System", 22nd SIGKDD Conference on Knowledge Discovery and Data Mining, 2016, https://arxiv.org/abs/1603.02754
#' @family engine
#' @name engine_xgboost
NULL
#' @rdname engine_xgboost
#' @export

engine_xgboost <- function(x,
                        booster = "gbtree",
                        learning_rate = 0.01,
                        nrounds = 10000,
                        gamma = 3,
                        reg_lambda = 0,
                        reg_alpha = 0,
                        max_depth = 2,
                        subsample = 0.8,
                        colsample_bytree = 0.4,
                        nthread = getOption('ibis.nthread'),
                        ...) {

  # Check whether xgboost package is available
  check_package('xgboost')
  if(!("xgboost" %in% loadedNamespaces()) || ('xgboost' %notin% sessionInfo()$otherPkgs) ) {
    try({requireNamespace('xgboost');attachNamespace("xgboost")},silent = TRUE)
  }

  # assert that arguments are valid
  assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                          inherits(x$background,'sf'),
                          is.character(booster) && booster %in% c("gbtree","gblinear"),
                          is.numeric(nrounds),
                          is.numeric(learning_rate) && (learning_rate > 0 && learning_rate < 1),
                          is.numeric(max_depth),
                          is.numeric(subsample) && (subsample > 0 && subsample <= 1),
                          is.numeric(colsample_bytree),
                          is.numeric(nthread)
  )

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
    booster = booster,
    nrounds = nrounds,
    eta = learning_rate,
    gamma = gamma,
    lambda = reg_lambda,
    alpha = reg_alpha,
    max_depth = max_depth,
    subsample = subsample,
    colsample_bytree = colsample_bytree,
    nthread = nthread,
    ...
    )

  # Print a message in case there is already an engine object
  if(!is.Waiver(x$engine)) myLog('[Setup]','yellow','Replacing currently selected engine.')

  # Set engine in distribution object
  x$set_engine(
    bdproto(
      "XGBOOST-Engine",
      Engine,
      name = "<XGBOOST>",
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
        fam <- switch(model$biodiversity[[1]]$family,
                      "poisson" = "count:poisson",
                      "binomial" = "binary:logistic",
                      model$biodiversity[[1]]$family
        )
        # If a poisson family is used, weight the observations by their exposure
        if(fam == "count:poisson"){
          # Get background layer
          bg <- self$get_data("template")
          assertthat::assert_that(!is.na(cellStats(bg,min)))

          abs <- create_pseudoabsence(
            env = model$predictors_object,
            presence = model$biodiversity[[1]]$observations,
            bias = settings$get('bias_variable'),
            template = bg,
            npoints = ifelse(ncell(bg)<10000,ncell(bg),10000),
            replace = TRUE
          )
          # Combine absence and presence and save
          abs$Intercept <- 1
          abs_observations <- abs[,c('x','y')]; abs_observations[['observed']] <- 0

          # Rasterize observed presences
          pres <- raster::rasterize(model$biodiversity[[1]]$observations[,c("x","y")],
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
          envs <- get_rastervalue(coords = obs[,c('x','y')],
                                  env = model$predictors_object$get_data(df = FALSE),
                                  rm.na = FALSE)
          envs$Intercept <- 1

          # Overwrite observations
          df <- rbind(envs[,c('x','y','Intercept', model$biodiversity[[1]]$predictors_names)],
                      abs[,c('x','y','Intercept', model$biodiversity[[1]]$predictors_names)])
          any_missing <- which(apply(df, 1, function(x) any(is.na(x))))
          if(length(any_missing)>0) abs_observations <- abs_observations[-any_missing,]
          df <- subset(df, complete.cases(df))
          # Overwrite observations
          model$biodiversity[[1]]$observations <- rbind(obs, abs_observations)

          # Pre-processing security checks
          assertthat::assert_that( all( model$biodiversity[[1]]$observations[['observed']] >= 0 ),
                                   any(!is.na(rbind(obs, abs_observations)[['observed']] )),
                                   nrow(df) == nrow(model$biodiversity[[1]]$observations)
          )

          # Define expectation as very small vector following Renner et al.
          w <- ppm_weights(df = df,
                           pa = model$biodiversity[[1]]$observations[['observed']],
                           bg = bg,
                           weight = 1 # Set those to 1 so that absences become ratio of pres/abs
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

        } else if(fam == "binary:logistic"){
          # calculating the case weights (equal weights)
          # the order of weights should be the same as presences and backgrounds in the training data
          prNum <- as.numeric(table(model$biodiversity[[1]]$observations[['observed']])["1"]) # number of presences
          bgNum <- as.numeric(table(model$biodiversity[[1]]$observations[['observed']])["0"]) # number of backgrounds
          w <- ifelse(model$biodiversity[[1]]$observations[['observed']] == 1, 1, prNum / bgNum)
          model$biodiversity[[1]]$expect <- w
          # Convert to numeric
          model$biodiversity[[1]]$observations$observed <- as.numeric( model$biodiversity[[1]]$observations$observed )
        }

        # Get Preds and convert to sparse matrix with set labels
        # FIXME: Support manual provision of data via xgb.DMatrix.save to save preprocessing time?
        train_cov <- as.matrix(
          model$biodiversity[[1]]$predictors[,model$biodiversity[[1]]$predictors_names]
        )
        labels <- model$biodiversity[[1]]$observations$observed

        # Create the sparse matrix
        df_train <- xgboost::xgb.DMatrix(data = train_cov,
                                         label = labels)
        # Prediction container
        pred_cov <- model$predictors[,model$biodiversity[[1]]$predictors_names]
        # Set target variables to bias_value for prediction if specified
        if(!is.Waiver(settings$get('bias_variable'))){
          for(i in 1:length(settings$get('bias_variable'))){
            if(settings$get('bias_variable')[i] %notin% colnames(pred_cov)) next()

            pred_cov[[settings$get('bias_variable')[i]]] <- settings$get('bias_value')[i]
          }
        }
        df_pred <- xgboost::xgb.DMatrix(data = as.matrix(pred_cov))

        if(fam == "count:poisson"){
          # Specifically for count poisson data we will set the areas
          # as an exposure offset for the base_margin
          xgboost::setinfo(df_train, "base_margin", log(w))
          xgboost::setinfo(df_pred, "base_margin", log(w_full))
          params$eval_metric <- "logloss"
        } else if(fam == 'binary:logistic'){
          params$eval_metric <- "logloss"
        }

        # Process and add priors if set
        if(!is.Waiver(model$priors)){
          assertthat::assert_that(
            all( model$priors$varnames() %in% model$biodiversity[[1]]$predictors_names )
          )
          # Match position of variables with monotonic constrains
          mc <- rep(0, ncol(train_cov))
          names(mc) <- colnames(train_cov)
          for(v in model$priors$varnames()){
            mc[v] <- switch (model$priors$get(v),
              'increasing' = 1, 'positive' = 1,
              'decreasing' = -1, 'negative' = -1,
              0
            )
          }
          # Save the monotonic constrain
          params$monotone_constraints <- mc
        }

        if(!is.Waiver(model$offset) ){
          # For the offset we simply add the (log-transformed) offset to the existing
          # Add exp at the end to backtransform
          of_train <- xgboost::getinfo(df_train, "base_margin") |> exp()
          of_pred <- xgboost::getinfo(df_pred, "base_margin") |> exp()
          # Add offset to full prediction and load vector

          # Respecify offset
          # (Set NA to 1 so that log(1) == 0)
          of <- model$offset; of[, "spatial_offset" ] <- ifelse(is.na(of[, "spatial_offset" ]), 1, of[, "spatial_offset"])
          of1 <- get_ngbvalue(coords = model$biodiversity[[1]]$observations[,c("x","y")],
                              env = of,
                              longlat = raster::isLonLat(self$get_data("template")),
                              field_space = c('x','y')
          )
          assertthat::assert_that(nrow(of1) == length(of_train),
                                  nrow(of) == length(of_pred))
          of_train <- of_train + of1[,"spatial_offset"]
          of_pred <- of_pred + of[,"spatial_offset"]

          # Set the new offset
          xgboost::setinfo(df_train, "base_margin", of_train)
          xgboost::setinfo(df_pred, "base_margin", of_pred)
        }

        # --- #
        # Save both training and predicting data in the engine data
        self$set_data("df_train", df_train)
        self$set_data("df_pred", df_pred)
        # --- #

        # Set objective
        params$objective <- fam

        self$set_control( params = params )

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

        # Get output raster
        prediction <- self$get_data('template')

        # Get parameters control
        params <- self$get_data('params')
        # Check only linear and reset to linear booster then
        if(settings$data$only_linear) params$booster <- "gblinear" else params$booster <- "gbtree"

        # All other needed data for model fitting
        df_train <- self$get_data("df_train")
        df_pred <- self$get_data("df_pred")
        w <- model$biodiversity[[1]]$expect # The expected weight

        assertthat::assert_that(
          is.null(w) || length(w) == nrow(df_train),
          all( colnames(df_train) %in% colnames(df_pred) )
        )

        # Get number of rounds from parameters
        nrounds <- params$nrounds;params$nrounds <- NULL

        # --- #
        # Pass this parameter possibly on from upper level
        # This implements a simple grid search for optimal parameter values
        # Using the training data only (!)
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

        # Fit the model.
        fit_xgb <- xgboost::xgboost(
          params = params,
          data = df_train,
          nrounds = nrounds,
          verbose = ifelse(verbose, 1, 0),
          print_every_n = 100
        )
        # --- #

        # Predict spatially
        if(!settings$get('inference_only')){
          # Messager
          if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Starting prediction...')

          # Make a prediction
          suppressWarnings(
            pred_xgb <- xgboost:::predict.xgb.Booster(
              object = fit_xgb,
              newdata = df_pred
              )
          )
          # Fill output with summaries of the posterior
          prediction[] <- pred_xgb
          names(prediction) <- 'mean'
          prediction <- raster::mask(prediction, self$get_data("template") )

        } else {
          # No prediction done
          prediction <- NULL
        }
        # Compute end of computation time
        settings$set('end.time', Sys.time())
        # Also append boosting control option to settings
        for(entry in names(params)) settings$set(entry, params[entry])

        # Definition of XGBOOST Model object ----
        # Create output
        out <- bdproto(
          "XGBOOST-Model",
          DistributionModel,
          id = new_id(),
          model = model,
          settings = settings,
          fits = list(
            "fit_best" = fit_xgb,
            "prediction" = prediction
          ),
          # Partial effects
          partial = function(self, x.var = NULL, plot = TRUE, ...){
            assertthat::assert_that(is.character(x.var) || is.null(x.var))
            check_package("pdp")
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
          },
          # Save the model object
          save = function(self, fname, what = "fit_best"){
            assertthat::assert_that(is.character(fname))
            xgboost::xgb.save( self$get_data(what), fname)
          }
        )
        return(out)
      }
    )
  ) # End of bdproto object
} # End of function
