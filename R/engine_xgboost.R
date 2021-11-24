#' @include bdproto-engine.R utils-spatial.R bdproto-distributionmodel.R
NULL

#' Engine for use of extreme gradient boosting (XGBoost)
#'
#' @description eXtreme gradient descent boosting for tree-based or linear boosting.
#' Flexible, yet powerful engine with many customization options.
#' The XGBoost engine supports multiple options to perform single and multiclass regression
#' and classification tasks. For a full list of options users are advised to have a look at the
#' [xgboost::xgb.train] help file and [https://xgboost.readthedocs.io](https://xgboost.readthedocs.io)
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param booster A [`character`] of the booster to use. Either "gbtree" or "gblinear" (Default: gblinear)
#' @param learning_rate [`numeric`] value indicating the learning rate (eta).
#' Lower values generally being better but more costly
#' @param nrounds [`numeric`] value giving the the maximum number of boosting iterations for cross-validation
#' @param gamma [`numeric`] A regularization parameter in the model. Lower values for better estimates (Default: 3)
#' Also see [reg_lambda] parameter for the L2 regularization on the weights
#' @param reg_lambda [`numeric`] L2 regularization term on weights (Default: 0)
#' @param reg_alpha [`numeric`] L1 regularization term on weights (Default: 0)
#' @param max_depth [`numeric`] The Maximum depth of a tree. (Default: 3)
#' @param subsample [`numeric`] The ratio used for subsampling to prevent overfitting. (Default: 0.8)
#' @param ... Other none specificed parameters. For constraints see [`XGBoostPrior`]
#' @seealso [xgboost::xgb.train]
#' @references  Tianqi Chen and Carlos Guestrin, "XGBoost: A Scalable Tree Boosting System", 22nd SIGKDD Conference on Knowledge Discovery and Data Mining, 2016, https://arxiv.org/abs/1603.02754
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
            env = model$predictors,
            presence = model$biodiversity[[1]]$observations,
            bias = settings$get('bias_variable'),
            template = bg,
            npoints = ifelse(ncell(bg)<10000,ncell(bg),10000),
            replace = TRUE
          )
          # Combine absence and presence and save
          abs_observations <- abs[,c('x','y')]; abs_observations[['observed']] <- 0

          # Rasterize observed presences
          pres <- raster::rasterize(model$biodiversity[[1]]$predictors[,c('x','y')], bg, fun = 'count', background = 0)
          # If family is not poisson, assume factor distribution
          # FIXME: Ideally this is better organized through family
          if(model$biodiversity[[1]]$family != 'poisson') pres[] <- ifelse(pres[]==1,1,0)
          obs <- cbind( data.frame(observed = raster::extract(pres, model$biodiversity[[1]]$observations[,c('x','y')])),
                        model$biodiversity[[1]]$observations[,c('x','y')] )
          model$biodiversity[[1]]$observations <- rbind(obs, abs_observations)

          # Format and add predictors
          abs <- subset(abs, select = c('x','y', model$biodiversity[[1]]$predictors_names) )
          df <- rbind(model$biodiversity[[1]]$predictors, abs) %>%
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

        } else
          if(fam == "binary:logistic"){
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
        pred_cov <- as.matrix(
          model$predictors[,model$biodiversity[[1]]$predictors_names]
        )
        df_pred <- xgboost::xgb.DMatrix(data = pred_cov)

        if(fam == "count:poisson"){
          # Specifically for count poisson data we will set the areas
          # as an exposure offset for the base_margin
          setinfo(df_train, "base_margin", log(w))
          setinfo(df_pred, "base_margin", log(w_full))
          params$eval_metric <- "logloss"
        } else if(fam == 'binary:logistic'){
          params$eval_metric <- "logloss"
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

        # Get output raster
        prediction <- self$get_data('template')

        # Get parameters control
        params <- self$get_data('params')
        # Check only linear and reset to linear booster then
        if(settings$data$only_linear) params$booster <- "gblinear"

        verbose <- settings$get("verbose")

        # All other needed data for model fitting
        df_train <- self$get_data("df_train")
        df_pred <- self$get_data("df_pred")
        w <- model$biodiversity[[1]]$expect # The expected weight

        assertthat::assert_that(
          is.null(w) || length(w) == nrow(df_train),
          all( colnames(df_train) %in% colnames(df_pred) )
        )

        if('offset' %in% names(model$biodiversity[[1]]) ){
          # Add offset to full prediction and load vector
          stop("TBD")
          # TBD
        }
        # Get number of rounds from parameters
        nrounds <- params$nrounds;params$nrounds <- NULL

        # --- #
        # Find the optimal number of rounds using 5-fold cross validation
        # TODO: Implement proper parameter tuning
        # https://stats.stackexchange.com/questions/171043/how-to-tune-hyperparameters-of-xgboost-trees
        fit_cv <- xgboost::xgb.cv(
          params = params,
          data = df_train,
          verbose = 0,
          nrounds = nrounds, nfold = 5,
          showsd = TRUE,      # standard deviation of loss across folds
          stratified = TRUE,  # sample is unbalanced; use stratified samplin
          maximize = FALSE,
          early_stopping_rounds = 10
        )

        # Fit the model.
        fit_xgb <- xgboost::xgboost(
          params = params,
          data = df_train,
          nrounds = fit_cv$best_iteration,
          verbose = verbose,
          print_every_n = 100
        )
        # --- #

        # Predict spatially
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

        # Create output
        out <- bdproto(
          "XGBOOST-Model",
          DistributionModel,
          id = new_id(),
          model = model,
          settings = settings,
          fits = list(
            "fit_best" = fit_xgb,
            "fit_cv" = fit_cv,
            "prediction" = prediction
          ),
          # Partial effects
          partial = function(self, x.vars, doplot = TRUE, ...){
            assertthat::assert_that(is.character(x.vars),
                                    !missing(x.vars))
            check_package("pdp")
            mod <- self$get_data('fit_best')
            df <- self$model$biodiversity[[length( self$model$biodiversity )]]$predictors
            df <- subset(df, select = mod$feature_names)

            # Check that variables are in
            assertthat::assert_that(x.vars %in% colnames(df),
                                    msg = 'Variable not in predicted model' )

            p1 <- pdp::partial(mod, pred.var = x.vars, ice = FALSE, center = TRUE,
                          plot = FALSE, rug = TRUE, train = df)
            names(p1) <- c("partial", "yhat")

            if(doplot){
              # Make a plot
              ggplot2::ggplot(data = p1, ggplot2::aes(x = partial, y = yhat)) +
                ggplot2::geom_line()+
                ggplot2::labs(x = x.vars, y = expression(hat(y)))
            }

            # Return the data
            return(p1)

          },
          # Spatial partial dependence plot option from embercardo
          spartial = function(self, predictors, x.vars = NULL, equal = FALSE, smooth = 1, transform = TRUE){
            stop("TBD")
            model <- self$get_data('fit_best')
            assertthat::assert_that(x.vars %in% attr(model$fit$data@x,'term.labels'),
                                    msg = 'Variable not in predicted model' )

            if( self$model$biodiversity[[1]]$family != 'binomial' && transform) warning('Check whether transform should not be set to False!')

            # Calculate
            p <- bart_partial_space(model, predictors, x.vars, equal, smooth, transform)

            cols <- c("#000004FF","#1B0C42FF","#4B0C6BFF","#781C6DFF","#A52C60FF","#CF4446FF","#ED6925FF","#FB9A06FF","#F7D03CFF","#FCFFA4FF")
            plot(p, col = cols, main = paste0(x.vars, collapse ='|'))
            # Also return spatial
            return(p)
          },
          # Engine-specific projection function
          project = function(self, newdata, summary = 'mean'){
            assertthat::assert_that(!missing(newdata),
                                    is.data.frame(newdata))

            # Define rowids as those with no missing data
            newdata$rowid <- rownames(newdata)
            newdata <- subset(newdata, complete.cases(newdata))
            # Make a prediction
            suppressWarnings(
              pred_bart <- dbarts:::predict.bart(object = self$get_data('fit_best'),
                                                 newdata = newdata,
                                                 type = 'response')
            )
            # Fill output with summaries of the posterior
            prediction <- emptyraster(self$fits$prediction) # Background
            prediction[as.numeric(newdata$rowid)] <- apply(pred_bart, 2, mean)

            # TODO: Generalize uncertainty prediction
            # # Summarize quantiles and sd from posterior
            # ms <- as.data.frame(
            #   cbind( matrixStats::colQuantiles(pred_bart, probs = c(.05,.5,.95)),
            #          matrixStats::colSds(pred_bart)
            #   )
            # )
            # names(ms) <- c('0.05ci','0.5ci','0.95ci','sd')
            # # Add them
            # for(post in names(ms)){
            #   prediction2 <- self$get_data('template')
            #   prediction2[as.numeric(full$cellid)] <- ms[[post]]; names(prediction2) <- post
            #   prediction <- raster::addLayer(prediction, prediction2)
            #   rm(prediction2)
            # }
            return(prediction)
          }
        )
        return(out)
      },
      # Save the model object
      save = function(self, fname, what = "fit_best"){
        assertthat::assert_that(is.character(fname))
        xgboost::xgb.save( self$get_data(what), fname)
      }
    )
  ) # End of bdproto object
} # End of function
