#' @include class-engine.R class-distributionmodel.R
NULL

#' Engine for extreme gradient boosting (XGBoost)
#'
#' @description Allows to estimate eXtreme gradient descent boosting for
#' tree-based or linear boosting regressions. The XGBoost engine is a
#' flexible, yet powerful engine with many customization options, supporting
#' multiple options to perform single and multi-class regression and
#' classification tasks. For a full list of options users are advised to have
#' a look at the [xgboost::xgb.train] help file and
#' [https://xgboost.readthedocs.io](https://xgboost.readthedocs.io).
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param booster A [`character`] of the booster to use. Either \code{"gbtree"}
#' or \code{"gblinear"} (Default: \code{gblinear})
#' @param learning_rate [`numeric`] value indicating the learning rate (eta).
#' Lower values generally being better but also computationally more costly.
#' (Default: \code{1e-3})
#' @param iter [`numeric`] value giving the the maximum number of boosting
#' iterations for cross-validation (Default: \code{8e3L}).
#' @param gamma [`numeric`] A regularization parameter in the model. Lower
#' values for better estimates (Default: \code{3}). Also see \code{"reg_lambda"}
#' parameter for the L2 regularization on the weights
#' @param reg_lambda [`numeric`] L2 regularization term on weights (Default: \code{0}).
#' @param reg_alpha [`numeric`] L1 regularization term on weights (Default: \code{0}).
#' @param max_depth [`numeric`] The Maximum depth of a tree (Default: \code{3}).
#' @param subsample [`numeric`] The ratio used for subsampling to prevent overfitting.
#' Also used for creating a random tresting dataset (Default: \code{0.75}).
#' @param colsample_bytree [`numeric`] Sub-sample ratio of columns when constructing
#' each tree (Default: \code{0.4}).
#' @param min_child_weight [`numeric`] Broadly related to the number of instances
#' necessary for each node (Default: \code{3}).
#' @param nthread [`numeric`] on the number of CPU-threads to use.
#' @param ... Other none specified parameters.
#'
#' @details The default parameters have been set relatively conservative as to
#' reduce overfitting.
#'
#' XGBoost supports the specification of monotonic constraints and interaction
#' constraints on certain variables. Within ibis these are possible via
#' [`XGBPrior`] and [`XGBInteractionPrior`], respectively. However constraints
#' are available only for the \code{"gbtree"} baselearners.
#'
#' @note
#' *'Machine learning is statistics minus any checking of models and assumptions‘* ~ Brian D. Ripley, useR! 2004, Vienna
#'
#' @returns An [Engine].
#'
#' @references
#' * Tianqi Chen and Carlos Guestrin, "XGBoost: A Scalable Tree Boosting System",
#' 22nd SIGKDD Conference on Knowledge Discovery and Data Mining, 2016,
#' https://arxiv.org/abs/1603.02754
#'
#' @seealso [xgboost::xgb.train]
#' @family engine
#'
#' @examples
#' \dontrun{
#' # Add xgboost as an engine
#' x <- distribution(background) |> engine_xgboost(iter = 4000)
#' }
#'
#' @name engine_xgboost
NULL

#' @rdname engine_xgboost
#' @export
engine_xgboost <- function(x,
                        booster = "gbtree",
                        iter = 8e3L,
                        learning_rate = 1e-3,
                        gamma = 6,
                        reg_lambda = 0,
                        reg_alpha = 0,
                        max_depth = 2,
                        subsample = 0.75,
                        colsample_bytree = 0.4,
                        min_child_weight = 3,
                        nthread = getOption('ibis.nthread'),
                        ...) {

  # Check whether xgboost package is available
  check_package('xgboost')
  if(!("xgboost" %in% loadedNamespaces()) || ('xgboost' %notin% utils::sessionInfo()$otherPkgs) ) {
    try({requireNamespace('xgboost');attachNamespace("xgboost")},silent = TRUE)
  }

  # assert that arguments are valid
  assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                          inherits(x$background,'sf'),
                          is.character(booster) && booster %in% c("gbtree","gblinear"),
                          is.numeric(iter),
                          is.numeric(learning_rate) && (learning_rate > 0 && learning_rate < 1),
                          is.numeric(max_depth),
                          is.numeric(subsample) && (subsample > 0 && subsample <= 1),
                          is.numeric(colsample_bytree),
                          is.numeric(nthread)
  )

  # Create a background raster
  template <- create_background(x)

  # mask template where all predictor layers are NA; change na.rm = FALSE for complete.cases
  if (!is.Waiver(x$predictors)){
    if(x$predictors$is_spatial()){
      template <- terra::mask(template, sum(x$predictors$get_data(), na.rm = TRUE))
    }
  }

  # Set up the parameter list
  params <- list(
    booster = booster,
    nrounds = iter,
    eta = learning_rate,
    gamma = gamma,
    lambda = reg_lambda,
    alpha = reg_alpha,
    max_depth = max_depth,
    subsample = subsample,
    colsample_bytree = colsample_bytree,
    min_child_weight = min_child_weight,
    nthread = nthread,
    ...
    )

  # Print a message in case there is already an engine object
  if(!is.Waiver(x$engine)) myLog('[Setup]','yellow','Replacing currently selected engine.')

  # Define new engine object of class
  eg <- Engine

  # Dummy function for spatial latent effects
  eg$set("public", "calc_latent_spatial", function(type = NULL, priors = NULL){new_waiver()},overwrite = TRUE)

  # Dummy function for getting the equation of latent effects
  eg$set("public", "get_equation_latent_spatial", function(method){new_waiver()},overwrite = TRUE)

  # Function to respecify the control parameters
  eg$set("public", "set_control", function(params){
    assertthat::assert_that(is.list(params))
    # Overwrite existing
    self$data$params <- params
    invisible()
  },overwrite = TRUE)

  # Setup function
  eg$set("public", "setup", function(model, settings = NULL, ...){
    # Simple security checks
    assertthat::assert_that(
      assertthat::has_name(model, 'background'),
      assertthat::has_name(model, 'biodiversity'),
      inherits(settings,'Settings') || is.null(settings),
      nrow(model$predictors) == terra::ncell(self$get_data('template')),
      !is.Waiver(self$get_data("params")),
      length(model$biodiversity) == 1 # Only works with single likelihood. To be processed separately
    )
    # Messenger
    if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','green','Engine setup.')

    # Get parameters
    params <- self$data$params

    # Distribution specific procedure
    fam <- switch(model$biodiversity[[1]]$family,
                  "poisson" = "count:poisson",
                  "binomial" = "binary:logistic",
                  model$biodiversity[[1]]$family
    )

    # Change the number of variables included if custom equation is used
    if(!is.Waiver(model$biodiversity[[1]]$equation)){
      form <- model$biodiversity[[1]]$equation
      varn <- model$biodiversity[[1]]$predictors_names[which( model$biodiversity[[1]]$predictors_names %in% formula_terms(form) )]
      assertthat::assert_that(length(varn)>0)
      # Match to existing ones and remove those not covered
      model$biodiversity[[1]]$predictors_names <- model$biodiversity[[1]]$predictors_names[match(varn, model$biodiversity[[1]]$predictors_names)]
      model$biodiversity[[1]]$predictors_types <- subset(model$biodiversity[[1]]$predictors_types,
                                                         predictors %in% model$biodiversity[[1]]$predictors_names)
    }

    # If a poisson family is used, weight the observations by their exposure
    if(fam == "count:poisson" && model$biodiversity[[1]]$type == "poipo"){
      # Get background layer
      bg <- self$get_data("template")
      assertthat::assert_that(!is.na(terra::global(bg, "min", na.rm = TRUE)))

      # Add pseudo-absence points
      suppressMessages(
        presabs <- add_pseudoabsence(df = model$biodiversity[[1]]$observations,
                                     field_occurrence = 'observed',
                                     template = bg,
                                     settings = model$biodiversity[[1]]$pseudoabsence_settings)
      )
      if(inherits(presabs, 'sf')) presabs <- presabs |> sf::st_drop_geometry()

      # Sample environmental points for absence only points
      abs <- subset(presabs, observed == 0)
      # Re-extract environmental information for absence points
      envs <- get_rastervalue(coords = abs[,c('x','y')],
                              env = model$predictors_object$get_data(df = FALSE),
                              rm.na = FALSE)
      if(assertthat::has_name(model$biodiversity[[1]]$predictors, "Intercept")){ envs$Intercept <- 1}

      # Format out
      df <- rbind(model$biodiversity[[1]]$predictors[,c('x','y','Intercept', model$biodiversity[[1]]$predictors_names)],
                  envs[,c('x','y','Intercept', model$biodiversity[[1]]$predictors_names)] )
      any_missing <- which(apply(df, 1, function(x) any(is.na(x))))
      if(length(any_missing)>0){
        presabs <- presabs[-any_missing,] # This works as they are in the same order
        model$biodiversity[[1]]$expect <- model$biodiversity[[1]]$expect[-any_missing]
      }
      df <- subset(df, stats::complete.cases(df))
      assertthat::assert_that(nrow(presabs) == nrow(df))

      # Check that expect matches
      if(length(model$biodiversity[[1]]$expect)!=nrow(df)){
        # Fill the absences with 1 as multiplier. This works since absences follow the presences
        model$biodiversity[[1]]$expect <- c( model$biodiversity[[1]]$expect,
                                             rep(1, nrow(presabs)-length(model$biodiversity[[1]]$expect) ))
      }

      # Overwrite observation data
      model$biodiversity[[1]]$observations <- presabs

      # Will expectations with 1 for rest of data points
      if(length(model$biodiversity[[1]]$expect)!= nrow(model$biodiversity[[1]]$observations)){
        model$biodiversity[[1]]$expect <- c(model$biodiversity[[1]]$expect,
                                            rep(1, nrow(model$biodiversity[[1]]$observations) - length(model$biodiversity[[1]]$expect))
        )
      }
      # Preprocessing security checks
      assertthat::assert_that( all( model$biodiversity[[1]]$observations[['observed']] >= 0 ),
                               any(!is.na(presabs[['observed']])),
                               length(model$biodiversity[[1]]$expect)==nrow(model$biodiversity[[1]]$observations),
                               nrow(df) == nrow(model$biodiversity[[1]]$observations)
      )

      # Add offset if existent
      if(!is.Waiver(model$offset)){
        # ofs <- get_ngbvalue(coords = df[,c('x','y')],
        #                     env =  model$offset,
        #                     longlat = terra::is.lonlat(bg),
        #                     field_space = c('x','y')
        # )
        ofs <- get_rastervalue(coords = df[,c('x','y')],
                               env = model$offset_object,
                               rm.na = FALSE)
        # Rename to spatial offset
        names(ofs)[which(names(ofs)==names(model$offset_object))] <- "spatial_offset"
        model$biodiversity[[1]]$offset <- ofs
      }

      # Define expectation as very small vector following Renner et al.
      w <- ppm_weights(df = df,
                       pa = model$biodiversity[[1]]$observations[['observed']],
                       bg = bg,
                       weight = 1e-6 # Set those to 1 so that absences become ratio of pres/abs
      )
      assertthat::assert_that(length(w) == nrow(df))

      model$biodiversity[[1]]$predictors <- df
      model$biodiversity[[1]]$expect <- w * (1/model$biodiversity[[1]]$expect)

      # Get for the full dataset
      pres <- terra::rasterize(x = guess_sf( model$biodiversity[[1]]$observations[,c("x","y")] ),
                               y = bg, fun = 'length', background = 0)
      w_full <- ppm_weights(df = model$predictors,
                            pa = pres[],
                            bg = bg,
                            weight = 1 # Set those to 1 so that absences become ratio of pres/abs
      )
      # Multiply with first weight value
      w_full <- w_full * (1/unique(model$biodiversity[[1]]$expect)[1])
      assertthat::assert_that(
        !anyNA(w_full), all(is.finite(log(w_full))),
        !anyNA(w_full),
        length(w_full) == nrow(model$predictors)
      )

    } else if(fam == "binary:logistic"){
      # calculating the case weights (equal weights)
      # the order of weights should be the same as presences and backgrounds in the training data
      prNum <- as.numeric(table(model$biodiversity[[1]]$observations[['observed']])["1"]) # number of presences
      bgNum <- as.numeric(table(model$biodiversity[[1]]$observations[['observed']])["0"]) # number of backgrounds
      w <- ifelse(model$biodiversity[[1]]$observations[['observed']] == 1, 1, prNum / bgNum)
      model$biodiversity[[1]]$expect <- w * model$biodiversity[[1]]$expect
      # Convert to numeric
      model$biodiversity[[1]]$observations$observed <- as.numeric( model$biodiversity[[1]]$observations$observed )
    }

    # Get Preds and convert to sparse matrix with set labels
    # FIXME: Support manual provision of data via xgb.DMatrix.save to save preprocessing time?
    train_cov <- model$biodiversity[[1]]$predictors[,model$biodiversity[[1]]$predictors_names]
    # Check if there any factors, if yes split up
    if(any(model$biodiversity[[1]]$predictors_types$type=='factor')){
      vf <- model$biodiversity[[1]]$predictors_types$predictors[which(model$biodiversity[[1]]$predictors_types$type == "factor")]
      # explode all present factors
      for (i in 1:length(vf)) {
        z <- explode_factor(train_cov[[vf[i]]], name = vf[i])
        # Remove variables from train_cov
        train_cov[[vf[i]]] <- NULL
        train_cov <- cbind(train_cov, z)
        # append to types data.frame
        model$biodiversity[[1]]$predictors_types <- rbind(model$biodiversity[[1]]$predictors_types,
                                                          data.frame(predictors = colnames(z), type = "numeric"))
      }
      # update containers
      model$biodiversity[[1]]$predictors <- train_cov # Save new in model object
      model$biodiversity[[1]]$predictors_types <- dplyr::filter(model$biodiversity[[1]]$predictors_types,
                                                                !predictors %in% vf)
      # Update biodiversity-level predictor names to match exploded training columns
      model$biodiversity[[1]]$predictors_names <- colnames(train_cov)
    }
    # Save the feature names used for training (after factor explosion)
    # model$predictors_names is updated later after model$predictors is also exploded
    train_feature_names <- colnames(train_cov)
    train_cov <- as.matrix( train_cov )
    labels <- model$biodiversity[[1]]$observations$observed

    # ---- #
    # Create the subsample based on the subsample parameter for all presence data
    # if(model$biodiversity[[1]]$type == "poipo"){
    #   ind <- sample(which(labels>0), size = params$subsample * length(which(labels>0)) )
    #   ind2 <- which( which(labels>0) %notin% ind )
    #   ind_ab <- which(labels==0)
    #   ind_train <- c(ind, ind_ab); ind_test <- c(ind2, ind_ab)
    # } else {
    #   ind_train <- sample(1:length(labels), size = params$subsample * length(labels) )
    #   ind_test <- which((1:length(labels)) %notin% ind_train )
    # }
    # Create the sparse matrix for training and testing data
    df_train <- xgboost::xgb.DMatrix(data = train_cov,
                                     label = labels#[ind_train]
    )
    # df_test <- xgboost::xgb.DMatrix(data = train_cov[c(ind_test),],
    #                                 label = labels[c(ind_test)]
    # )
    # --- #
    # If predictions are to be used, define here
    if(!settings$get('inference_only')){
      # Prediction container - use model$predictors_names (pre-explosion) since model$predictors
      # still has the original factor columns at this point
      pred_cov <- model$predictors[,model$predictors_names]
      if(any(model$predictors_types$type=='factor')){
        vf <- model$predictors_types$predictors[which(model$predictors_types$type == "factor")]
        # Get factors
        for (i in 1:length(vf)) {
          z <- explode_factor(pred_cov[[vf[i]]], name = vf[i])
          # Remove variables from train_cov and append
          pred_cov[[vf[i]]] <- NULL
          pred_cov <- cbind(pred_cov, z)
          model$predictors_types <- rbind(model$predictors_types,
                                          data.frame(predictors = colnames(z), type = "numeric"))
        }

        model$predictors <- pred_cov # Save new in model object

        model$predictors_names <- colnames(pred_cov)
        assertthat::assert_that(all( train_feature_names %in% colnames(pred_cov) ))
      }
      pred_cov <- as.matrix( pred_cov )
      # Ensure that the column names are identical for both
      pred_cov <- pred_cov[, colnames(train_cov)]

      # Clamp?
      if( settings$get("clamp") ) pred_cov <- clamp_predictions(model, pred_cov)

      # Set target variables to bias_value for prediction if specified
      if(!is.Waiver(settings$get('bias_variable'))){
        for(i in 1:length(settings$get('bias_variable'))){
          if(settings$get('bias_variable')[i] %notin% colnames(pred_cov)){
            if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','red','Did not find bias variable in prediction object!')
            next()
          }
          pred_cov[,settings$get('bias_variable')[i]] <- settings$get('bias_value')[i]
        }
      }
      df_pred <- xgboost::xgb.DMatrix(data = as.matrix(pred_cov))
      assertthat::assert_that(all(colnames(df_train) == colnames(df_pred)))

    } else {
      # inference_only: still update model$predictors_names to match the
      # exploded training features so that later projection uses correct names
      model$predictors_names <- train_feature_names
    }

    if(fam == "count:poisson"){
      # Specifically for count poisson data we will set the areas
      assertthat::assert_that(all(is.finite(log(w))))
      # as an exposure offset for the base_margin
      xgboost::setinfo(df_train, "base_margin", log(w))
      # xgboost::setinfo(df_test, "base_margin", log(w[ind_test]))
      if(!settings$get('inference_only')){
        assertthat::assert_that(nrow(df_pred) == length(w_full),
                                all(is.finite(log(w_full))) )
        xgboost::setinfo(df_pred, "base_margin", log(w_full))
      }
      params$eval_metric <- "logloss"
    } else if(fam == 'binary:logistic'){
      params$eval_metric <- "logloss"
    }

    # Process and add XGBoost priors if set
    xgb_prior_params <- format_xgboost_priors(
      priors = model$priors,
      feature_names = colnames(train_cov),
      # The final booster is selected in train() from the only_linear setting.
      booster = "gbtree"
    )
    for(entry in names(xgb_prior_params)) params[[entry]] <- xgb_prior_params[[entry]]

    if(!is.Waiver(model$offset) ){
      # Set offset to 1 (log(0)) in case nothing is found
      if(is.null(xgboost::getinfo(df_train, "base_margin"))) {
        of_train <- rep(1, nrow(model$biodiversity[[1]]$observations[,c("x","y")]))
        if(!settings$get('inference_only')) of_pred <- rep(1, nrow(model$offset))
      } else {
        # For the offset we simply add the (log-transformed) offset to the existing one
        # given that for example log(2*3) == log(2) + log(3)
        of_train <- xgboost::getinfo(df_train, "base_margin")
        # of_test <- xgboost::getinfo(df_test, "base_marginfit_xgb") |> exp()
        if(!settings$get('inference_only')){
          of_pred <- xgboost::getinfo(df_pred, "base_margin")
        }
      }
      # -- Add offset to full prediction and load vector --

      # Respecify offset
      # (Set NA to 1 so that log(1) == 0)
      of <- model$offset; of[, "spatial_offset" ] <- ifelse(is.na(of[, "spatial_offset" ]), 1, of[, "spatial_offset"])
      of1 <- get_rastervalue(coords = model$biodiversity[[1]]$observations[,c("x","y")],
                             env = model$offset_object,
                             rm.na = FALSE
      )
      names(of1)[which(names(of1)==names(model$offset_object))] <- "spatial_offset"
      # of2 <- get_rastervalue(coords = model$biodiversity[[1]]$observations[ind_test,c("x","y")],
      #                        env = model$offset_object,
      #                        rm.na = FALSE
      #                     # longlat = terra::is.lonLat(self$get_data("template")),
      #                     # field_space = c('x','y')
      # )
      # names(of2)[which(names(of2)==names(model$offset_object))] <- "spatial_offset"
      assertthat::assert_that(nrow(of1) == length(of_train)) # nrow(of2) == length(of_test)
      of_train <- of_train + of1[,"spatial_offset"]
      # of_test <- of_test + of2[,"spatial_offset"]
      # Check that values are valid
      assertthat::assert_that(all(is.finite(of_train)), !anyNA(of_train), )

      # Set the new offset
      xgboost::setinfo(df_train, "base_margin", ( of_train ))
      # xgboost::setinfo(df_test, "base_margin", of_test)
      if(!settings$get('inference_only')){
        assertthat::assert_that(nrow(of) == length(of_pred))
        of_pred <- of_pred + of[,"spatial_offset"]

        assertthat::assert_that(
          !anyNA(of_pred), all(is.finite(of_pred))
        )
        # Set the new offset
        xgboost::setinfo(df_pred, "base_margin", ( of_pred ))
      }
    }

    # --- #
    # Save both training and predicting data in the engine data
    self$set_data("df_train", df_train)
    # self$set_data("df_test", df_test)
    if(!settings$get('inference_only')) self$set_data("df_pred", df_pred)
    # --- #

    # Set objective
    params$objective <- fam

    self$set_control( params = params )

    # Instead of invisible return the model object
    return( model )
  },overwrite = TRUE)

  # Training function
  eg$set("public", "train", function(model, settings, ...){
    assertthat::assert_that(
      inherits(settings,'Settings'),
      is.list(model),length(model)>1,
      # Check that model id and setting id are identical
      settings$modelid == model$id
    )
    # Get name
    name <- model$biodiversity[[1]]$name

    # Messenger
    if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','green', paste0('Starting fitting: ', name))

    # Verbosity
    verbose <- settings$get("verbose")

    # Get output raster
    prediction <- self$get_data('template')

    # Get parameters control
    params <- self$get_data('params')
    assertthat::assert_that(is.list(params))

    # Check only linear and reset to linear booster then
    if(settings$get("only_linear")) params$booster <- "gblinear" else params$booster <- "gbtree"
    # Check that link function and objective is changed if needed
    li <- model$biodiversity[[1]]$link
    if(!is.null(li)){
      if(model$biodiversity[[1]]$family=="binomial"){
        li <- match.arg(li, c("logit", "cloglog"),several.ok = FALSE)
        if(li=="cloglog") params$objective <- "binary:logitraw"
      } else {
        if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','red',paste0("Package does not support custom link functions. Ignored!"))
      }
    }

    # All other needed data for model fitting
    df_train <- self$get_data("df_train")
    # df_test <- self$get_data("df_test")
    df_pred <- self$get_data("df_pred")
    w <- model$biodiversity[[1]]$expect # The expected weight

    if(!settings$get('inference_only')){
      assertthat::assert_that(
        !is.null(df_pred),
        all( colnames(df_pred) %in% colnames(df_train) )
      )
    }

    # Get number of rounds from parameters
    nrounds <- params$nrounds;params$nrounds <- NULL

    # make sure to only include priors of current model
    if (!is.null(params$monotone_constraints)) {
      params$monotone_constraints <- params$monotone_constraints[model$biodiversity[[1]]$predictors_names]
    }

    # --- #
    # Pass this parameter possibly on from upper level
    # This implements a simple grid search for optimal parameter values
    # Using the training data only (!)
    if(settings$get('optim_hyperparam')){
      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','green','Starting hyperparameters search...')

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
          subsample = stats::runif(1, .7, 1),
          colsample_bytree = stats::runif(1, .6, 1),
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
          test_xgb <- xgboost::xgb.train(
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
      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','green','Crossvalidation for determining early stopping rule.')

      fit_cv <- xgboost::xgb.cv(
        params = params,
        data = df_train,
        verbose = ifelse(verbose, 1, 0),
        print_every_n = 100,
        nrounds = nrounds, nfold = 5,
        showsd = TRUE,      # standard deviation of loss across folds
        stratified = TRUE,  # sample is unbalanced; use stratified sampling
        maximize = FALSE,
        early_stopping_rounds = 10
      )
      # Set new number of rounds
      nround <- fit_cv$best_iteration
    }

    # Remove unneeded parameters
    if(params$booster == "gblinear"){
      tree_constraints <- c("monotone_constraints", "interaction_constraints")
      existing_tree_constraints <- intersect(tree_constraints, names(params))
      if(length(existing_tree_constraints) > 0){
        warning(
          "XGBoost monotone and interaction constraints are only supported for tree boosters. Ignoring xgboost priors for gblinear.",
          call. = FALSE
        )
        params[existing_tree_constraints] <- NULL
      }
      params[c("colsample_bytree", "gamma", "max_depth", "min_child_weight", "subsample")] <- NULL
    }
    # Fit the model.
    # watchlist <- list(train = df_train,test = df_test)
    fit_xgb <- xgboost::xgb.train(
      params = params,
      data = df_train,
      # watchlist = watchlist,
      nrounds = nrounds,
      verbose = ifelse(verbose, 1, 0),
      print_every_n = 100
    )
    # --- #

    # Predict spatially
    if(!settings$get('inference_only')){
      # Messager
      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','green','Starting prediction...')

      # Make a prediction
      suppressWarnings(
        pred_xgb <- predict(
          object = fit_xgb,
          newdata = df_pred
        )
      )
      if(params$objective=="binary:logitraw") pred_xgb <- ilink(pred_xgb, "cloglog")

      # Fill output with summaries of the posterior
      prediction[] <- pred_xgb
      names(prediction) <- 'mean'
      prediction <- terra::mask(prediction, self$get_data("template") )
    } else {
      # No prediction done
      prediction <- NULL
    }
    # Compute end of computation time
    settings$set('end.time', Sys.time())
    # Also append boosting control option to settings
    for(entry in names(params)) settings$set(entry, params[entry])

    # Definition of XGBOOST Model object ----
    obj <- DistributionModel # Make a copy to set new functions

    # Partial effects
    obj$set("public", "partial", function(x.var = NULL, constant = NULL, variable_length = 100, values = NULL,
                           newdata = NULL, plot = TRUE, type = "response"){
      assertthat::assert_that(is.character(x.var) || is.null(x.var))
      if(!is.null(constant)) message("Constant is ignored for xgboost!")
      # variable_length = 100; values = NULL; newdata = NULL; plot = TRUE; type = "response"

      # Settings
      settings <- self$settings
      mod <- self$get_data('fit_best')
      model <- self$model

      # Set type
      if(is.null(type)) type <- self$settings$get("type")
      type <- match.arg(type, c("link", "response"), several.ok = FALSE)
      settings$set("type", type)

      # Match x.var to argument
      if(is.null(x.var)){
        x.var <- model$predictors_names
      } else {
        x.var <- match.arg(x.var, model$predictors_names, several.ok = TRUE)
      }

      df <- model$biodiversity[[length( model$biodiversity )]]$predictors
      df <- subset(df, select = model$biodiversity[[1]]$predictors_names)
      # Ensure all columns are numeric (factors should already be exploded during training)
      fac_cols <- vapply(df, is.factor, logical(1))
      if(any(fac_cols)){
        for(fc in names(df)[fac_cols]){
          df[[fc]] <- as.numeric(as.character(df[[fc]]))
        }
      }
      if(!is.null(newdata)){
        newdata <- subset(newdata, select = names(df))
        assertthat::assert_that(nrow(newdata)>1,ncol(newdata)>1,
                                all(names(df) %in% names(newdata)))
      }

      # Calculate range of predictors (only for numeric columns)
      rr <- apply(df, 2, function(x) range(as.numeric(x), na.rm = TRUE)) |> as.data.frame()

      if(is.null(newdata)){
        # if values are set, make sure that they cover the data.frame
        if(!is.null(values)){
          df2 <- list()
          df2[[x.var]] <- values
          # Then add the others
          for(var in colnames(df)){
            if(var %in% x.var) next()
            df2[[var]] <- mean(df[[var]], na.rm = TRUE)
          }
          df2 <- df2 |> as.data.frame()
          df2 <- df2[, x.var, drop = FALSE]
        } else {
          df2 <- list()
          for(i in x.var) {
            df2[[i]] <- base::as.data.frame(seq(rr[1,i],rr[2,i], length.out = variable_length))
          }
          df2 <- do.call(cbind, df2); names(df2) <- x.var
        }
      } else {
        # Assume that newdata container has all the variables for the grid
        df2 <- dplyr::select(newdata, dplyr::any_of(x.var))
      }

      # Get offset if set
      if(!is.Waiver(model$offset)){
        of <- model$offset$spatial_offset
      } else of <- new_waiver()

      # Inverse link function
      ilf <- switch(settings$get('type'),
                     "link" = NULL,
                     "response" = if(model$biodiversity[[1]]$family == 'poisson') exp else logistic
      )

      # Prediction function for xgboost
      xgb_predict_fun <- function(object, newdata) {
        dm <- xgboost::xgb.DMatrix(data = data.matrix(newdata))
        as.numeric(predict(object, dm))
      }

      pp <- data.frame()
      pb <- progress::progress_bar$new(total = length(x.var))
      for(v in x.var){
        # New wrapper specified in utils
        p1 <- compute_partial_dependence(
          object = mod, pred.var = v, pred.grid = df2, train = df,
          predict_fun = xgb_predict_fun,
          inv.link = ilf,
          offset = if(!is.Waiver(of)) of else NULL
        )
        names(p1) <- c("partial_effect", "mean")
        p1 <- cbind(variable = v, p1)
        pp <- rbind(pp, p1)
        if(length(x.var) > 1) pb$tick()
      }

      if(plot){
        # Make a plot
        g <- ggplot2::ggplot(data = pp, ggplot2::aes(x = partial_effect)) +
          ggplot2::theme_classic() +
          ggplot2::geom_line(ggplot2::aes(y = mean)) +
          ggplot2::facet_wrap(. ~ variable, scales = "free") +
          ggplot2::labs(x = "Variable", y = "Partial effect")
        print(g)
      }
      # Return the data
      return(pp)
    },overwrite = TRUE)

    # Spatial partial dependence plot
    obj$set("public", "spartial", function(x.var, constant = NULL, newdata = NULL, plot = TRUE, ...){
      assertthat::assert_that(is.character(x.var) || is.null(x.var),
                              "model" %in% names(self),
                              is.null(newdata) || is.data.frame(newdata))

      # Get data
      mod <- self$get_data('fit_best')
      model <- self$model
      settings <- self$settings
      # Feature names from the model (xgboost 3.x no longer stores feature_names on the booster)
      fn <- model$predictors_names
      x.var <- match.arg(x.var, fn, several.ok = FALSE)

      # Get predictor data (subset to training features only)
      df <- subset(model$predictors, select = fn)
      assertthat::assert_that(is.data.frame(df), ncol(df) > 0,
                              msg = "Could not obtain predictor data from model.")

      # Make template of target variable(s)
      template <- model_to_background(model)

      # Set all variables other than the target variable to constant
      if(!is.null(constant)){
        other_vars <- fn[fn %notin% x.var]
        df[!is.na(df[,x.var]), other_vars] <- constant
      }

      # Convert to properly aligned numeric matrix for xgb.DMatrix
      dm <- xgboost::xgb.DMatrix(data = data.matrix(df))

      # Spatial prediction contributions: predcontrib = TRUE calculates
      # contributions of each feature to individual predictions.
      # For "gblinear" booster, these are linear terms (beta * value).
      # For "gbtree" booster, these are SHAP values.
      pp <- predict(object = mod,
                    newdata = dm,
                    predcontrib = TRUE) |>
        as.data.frame()
      # Get only target variable
      pp <- subset(pp, select = x.var)
      # if(settings$get('objective')[[1]]=="binary:logitraw") pp <- ilink(pp, "cloglog")
      assertthat::assert_that(terra::ncell(template) == nrow(pp))

      # Fill output with summaries of the posterior
      template <- fill_rasters(as.data.frame(pp), template)
      names(template) <- 'mean'
      template <- terra::mask(template, model$background)

      if(plot){
        # Quick plot
        terra::plot(template, col = ibis_colours$ohsu_palette,
                    main = paste0(x.var, collapse ='|'))
      }
      # Also return spatial
      return(template)
    },overwrite = TRUE)

    # Engine-specific projection function
    obj$set("public", "project", function(newdata, layer = "mean"){
      assertthat::assert_that(!missing(newdata),
                              is.data.frame(newdata) || inherits(newdata, "xgb.DMatrix") )

      mod <- self$get_data('fit_best')
      # Get model object
      model <- self$model
      # Internal newdata copy
      newdata_copy <- newdata

      # Also get settings for bias values
      settings <- self$settings

      if(!inherits(newdata, "xgb.DMatrix")){
        # Feature names from the model (xgboost 3.x no longer stores feature_names on the booster)
        fn <- model$predictors_names

        # Explode factor variables in newdata to match training (if factors were exploded during setup)
        if(!is.null(model$biodiversity[[1]]$predictors_types)){
          fac_types <- model$biodiversity[[1]]$predictors_types
          # Check if any exploded factor columns are expected (pattern: "name.level")
          # by testing if fn contains names not in newdata but that look like exploded factors
          missing_fn <- fn[!fn %in% colnames(newdata)]
          if(length(missing_fn) > 0){
            # Try to identify original factor columns that need explosion
            for(col in colnames(newdata)){
              if(is.factor(newdata[[col]]) || is.character(newdata[[col]])){
                potential_names <- paste(col, levels(as.factor(newdata[[col]])), sep = ".")
                if(any(potential_names %in% missing_fn)){
                  z <- explode_factor(as.factor(newdata[[col]]), name = col)
                  newdata[[col]] <- NULL
                  newdata <- cbind(newdata, z)
                }
              }
            }
          }
        }

        # Check for predictor mismatches and provide informative warnings
        missing_in_newdata <- fn[!fn %in% colnames(newdata)]
        extra_in_newdata <- colnames(newdata)[!colnames(newdata) %in% c("x", "y", fn)]
        if(length(missing_in_newdata) > 0){
          warning(paste0("Predictors expected by model but missing in newdata: ",
                         paste(missing_in_newdata, collapse = ", ")))
        }
        if(length(extra_in_newdata) > 0 && getOption('ibis.setupmessages', default = TRUE)){
          myLog('[Estimation]','yellow',
                paste0('Extra predictors in newdata ignored: ', paste(extra_in_newdata, collapse = ", ")))
        }
        assertthat::assert_that(
          all( fn %in% colnames(newdata) ),
          msg = paste0("Feature name mismatch between model and newdata. Missing: ",
                       paste(fn[!fn %in% colnames(newdata)], collapse = ", "))
        )
        newdata <- subset(newdata, select = fn)

        # Clamp?
        if( settings$get("clamp") ) newdata <- clamp_predictions(model, newdata)

        if(!is.Waiver(settings$get('bias_variable'))){
          for(i in 1:length(settings$get('bias_variable'))){
            if(settings$get('bias_variable')[i] %notin% colnames(newdata)){
              if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','red','Did not find bias variable in prediction object!')
              next()
            }
            newdata[,settings$get('bias_variable')[i]] <- settings$get('bias_value')[i]
          }
        }
        newdata <- xgboost::xgb.DMatrix(data.matrix(newdata))
      } else {cli::cli_abort("Not implemented. Supply a data.frame as newdata!")}

      # Make a prediction
      suppressWarnings(
        pred_xgb <- predict(
          object = mod,
          newdata = newdata
        )
      )

      # Fill output with summaries of the posterior
      if(nrow(newdata)==nrow(model$predictors)){
        prediction <- try({model_to_background(model)}, silent = TRUE)
        prediction[] <- pred_xgb
        # Create mask using only numeric columns to avoid issues with factor/character columns
        mask_df <- newdata_copy[, intersect(c("x", "y", fn), colnames(newdata_copy))]
        for(mc in names(mask_df)){
          if(!is.numeric(mask_df[[mc]])) mask_df[[mc]] <- as.numeric(as.factor(mask_df[[mc]]))
        }
        mask_tmp <- sum(terra::rast(mask_df, crs = terra::crs(prediction)))
        prediction <- terra::mask(prediction, mask_tmp) # make sure to only use nonNA cells
      } else {
        assertthat::assert_that(utils::hasName(newdata_copy,"x")&&utils::hasName(newdata_copy,"y"),
                                msg = "Projection data.frame has no valid coordinates or differs in grain!")
        prediction <- try({
          terra::rast(newdata_copy[,c("x", "y")],
                      crs = terra::crs(model$background),
                      type = "xyz") |>
            emptyraster()
        }, silent = TRUE)
        prediction[] <- pred_xgb
      }

      return(prediction)
    },overwrite = TRUE)

    # Model convergence check
    obj$set("public", "has_converged", function(){
      fit <- self$get_data("fit_best")
      if(is.Waiver(fit)) return(FALSE)
      # Get evaluation log if existing
      if(utils::hasName(fit, "evaluation_log")){
        evl <- fit$evaluation_log
        if(fit$best_iteration >= (nrow(evl)-(nrow(evl)*.01))) return(FALSE)
      }
      return(TRUE)
    },overwrite = TRUE)

    # Residual function
    obj$set("public", "get_residuals", function(){
      # Get best object
      obj <- self$get_data("fit_best")
      if(is.Waiver(obj)) return(obj)
      message("Not yet implemented!")
      return(new_waiver())
      # Get residuals
      model <- self$model
      pred <- model$biodiversity[[length(model$biodiversity)]]
      predf <- pred$predictors |> subset(select = model$predictors_names)
      newdata <- xgboost::xgb.DMatrix(data.matrix(predf))

      fn <- predict(obj, newdata,type = "class")
      return(fn)
    },overwrite = TRUE)

    # Get coefficients
    obj$set("public", "get_coefficients", function(){
      # Returns a vector of the coefficients with direction/importance
      obj <- self$get_data('fit_best')
      # Simply use the weights from the importance estimates
      cofs <- xgboost::xgb.importance(model = obj) |>
        as.data.frame()
      cofs$Sigma <- NA
      if(!self$settings$get("only_linear")){
        cofs <- subset(cofs, select = c("Feature", "Gain", "Sigma"))
      }
      names(cofs) <- c("Feature", "Beta", "Sigma")
      return(cofs)
    },overwrite = TRUE)

    # Save the model object
    obj$set("public", "save", function(fname, what = "fit_best"){
      assertthat::assert_that(is.character(fname))
      xgboost::xgb.save( self$get_data(what), fname)
    },overwrite = TRUE)

    out <- obj$new(name = "XGBOOST-Model")
    out$id <- model$id
    out$model <- model
    out$settings <- settings
    out$fits <- list(
      "fit_best" = fit_xgb,
      "prediction" = prediction
    )

    return(out)
  },overwrite = TRUE)

  # Define engine object
  eg <- eg$new(engine = "XGBOOST-Engine", name = "<XGBOOST>")
  eg$data <- list(
    'template' = template,
    'params' = params
  )

  # Set engine in distribution object
  y <- x$clone(deep = TRUE)
  return( y$set_engine(eg) )
} # End of function
