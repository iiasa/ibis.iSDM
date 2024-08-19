#' @include class-engine.R class-distributionmodel.R
NULL

#' Use of Gradient Descent Boosting for model estimation
#'
#' @description Gradient descent boosting is an efficient way to optimize any
#' loss function of a generalized linear or additive model (such as the GAMs
#' available through the \code{"mgcv"} R-package). It furthermore automatically
#' regularizes the fit, thus the resulting model only contains the covariates
#' whose baselearners have some influence on the response. Depending on the type
#' of the \code{add_biodiversity} data, either poisson process models or
#' logistic regressions are estimated. If the \code{"only_linear"} term in
#' [train] is set to \code{FALSE}, splines are added to the estimation, thus
#' providing a non-linear additive inference.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param iter An [`integer`] giving the number of boosting iterations (Default: \code{2e3L}).
#' @param learning_rate A bounded [`numeric`] value between \code{0} and \code{1}
#' defining the shrinkage parameter.
#' @param empirical_risk method for empirical risk calculation. Available options
#' are \code{'inbag'}, \code{'oobag'} and \code{'none'}. (Default: \code{'inbag'}).
#' @param type The mode used for creating posterior predictions. Either making
#' \code{"link"}, \code{"response"} or \code{"class"} (Default: \code{"response"}).
#' @param ... Other variables or control parameters
#'
#' @details: This package requires the \code{"mboost"} R-package to be
#' installed. It is in philosophy somewhat related to the [engine_xgboost] and
#' \code{"XGBoost"} R-package, however providing some additional desirable
#' features that make estimation quicker and particularly useful for spatial
#' projections. Such as for instance the ability to specifically add spatial
#' baselearners via [add_latent_spatial] or the specification of monotonically
#' constrained priors via [GDBPrior].
#'
#' @note
#' The coefficients resulting from gdb with poipa data (Binomial) are only 0.5
#' of the typical coefficients of a logit model obtained via glm (see Binomial).
#'
#' @returns An engine.
#'
#' @references
#' * Hofner, B., Mayr, A., Robinzonov, N., & Schmid, M. (2014). Model-based boosting
#' in R: a hands-on tutorial using the R package mboost. Computational statistics, 29(1-2), 3-35.
#' * Hofner, B., Müller, J., Hothorn, T., (2011). Monotonicity-constrained species
#' distribution models. Ecology 92, 1895–901.
#' * Mayr, A., Hofner, B. and Schmid, M. (2012). The importance of knowing when
#' to stop - a sequential stopping rule for component-wise gradient boosting.
#' Methods of Information in Medicine, 51, 178–186.
#'
#' @family engine
#'
#' @examples
#' \dontrun{
#' # Add GDB as an engine
#' x <- distribution(background) |> engine_gdb(iter = 1000)
#' }
#'
#' @name engine_gdb
NULL

#' @rdname engine_gdb
#' @export
engine_gdb <- function(x,
                       iter = 2000,
                       learning_rate = 0.1,
                       empirical_risk = 'inbag',
                       type = "response",
                        ...) {
  # Check whether mboost package is available
  check_package('mboost')
  if(!("mboost" %in% loadedNamespaces()) || ('mboost' %notin% utils::sessionInfo()$otherPkgs) ) {
    try({requireNamespace('mboost');attachNamespace("mboost")},silent = TRUE)
    }

  # assert that arguments are valid
  assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                          inherits(x$background,'sf'),
                          is.numeric(iter),
                          is.numeric(learning_rate),
                          is.character(empirical_risk),
                          is.character(type),
                          empirical_risk %in% c('inbag','oobag','none')
                          )
  # Match type
  type <- match.arg(type, choices = c("link", "response", "class"), several.ok = FALSE)
  # Get background
  background <- x$background

  # Create a background raster
  if(is.Waiver(x$predictors)){
    # Create from background
    template <- terra::rast(
      ext = terra::ext(background),
      crs = terra::crs(background),
      res = c(diff( (sf::st_bbox(background)[c(1,3)]) ) / 100, # Simplified assumption for resolution
              diff( (sf::st_bbox(background)[c(1,3)]) ) / 100
             )
      )
  } else {
    # If predictor existing, use them
    template <- emptyraster(x$predictors$get_data() )
  }

  # Burn in the background
  template <- terra::rasterize(background, template, field = 0)

  # Set up boosting control
  bc <- mboost::boost_control(mstop = iter,
                              nu = learning_rate,
                              risk = empirical_risk
                              )

  # Set up the parameter list
  params <- list(
    type = type,
    ...
  )

  # Print a message in case there is already an engine object
  if(!is.Waiver(x$engine)) message('Replacing currently selected engine.')

  # Define new engine object of class
  eg <- Engine

  # Function to respecify the control parameters
  eg$set("public", "set_control", function(
    iter = 20,
    learning_rate = 0.1, # Set relatively low to not regularize too much
    empirical_risk = 'inbag',
    verbose = TRUE
  ){
    # Set up boosting control
    bc <- mboost::boost_control(mstop = iter,
                                nu = learning_rate,
                                risk = empirical_risk,
                                trace = verbose
    )
    # Overwrite existing
    self$data$bc <- bc

  },overwrite = TRUE)

  # Dummy function for latent factors
  eg$set("public", "calc_latent_spatial", function(...){ invisible()},overwrite = TRUE)

  # Get equation for spatial effect
  eg$set("public", "get_equation_latent_spatial", function(spatial_field = c('x','y'),
                                            df = 6, knots = 4,...){
    return(
      paste0(
        'bspatial(',spatial_field[1],',',spatial_field[2],', center = TRUE, df = ',df,', knots = ',knots,')',
        ' + ',
        'bols(',spatial_field[1],')', '+', 'bols(',spatial_field[2],')', '+', 'bols(',spatial_field[1],',',spatial_field[2],')'
      )
    )
  },overwrite = TRUE)

  # Setup function
  eg$set("public", "setup", function(model, settings = NULL, ...){
    # Simple security checks
    assertthat::assert_that(
      assertthat::has_name(model, 'background'),
      assertthat::has_name(model, 'biodiversity'),
      inherits(settings,'Settings') || is.null(settings),
      # Check that all predictors are present
      nrow(model$predictors) == terra::ncell(self$get_data('template')),
      length(model$biodiversity) == 1 # Only works with single likelihood. To be processed separately
    )
    # Add in case anything needs to be further prepared here
    # Messenger
    if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','green','Engine setup.')

    # Add pseudo-absence points if necessary
    # Include nearest predictor values for each
    if('poipo' == model$biodiversity[[1]]$type && model$biodiversity[[1]]$family == 'poisson') {

      # Get background layer
      bg <- self$get_data("template")
      assertthat::assert_that(is.Raster(bg), !is.na(terra::global(bg, "min", na.rm = TRUE)[,1]))

      # Add pseudo-absence points
      presabs <- add_pseudoabsence(df = model$biodiversity[[1]]$observations,
                                   field_occurrence = 'observed',
                                   template = bg,
                                   settings = model$biodiversity[[1]]$pseudoabsence_settings)
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
      if(length(any_missing)>0) {
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

      # Check that factors have been correctly set if any
      if(any(model$predictors_types$type=="factor")){
        df[,model$predictors_types$predictors[model$predictors_types$type=="factor"]] <-
          lapply(df[,model$predictors_types$predictors[model$predictors_types$type=="factor"]], as.factor)
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
        ofs <- get_rastervalue(coords = df[,c('x','y')],
                               env = model$offset_object,
                               rm.na = FALSE)
        # Rename to spatial offset
        names(ofs)[which(names(ofs)==names(model$offset_object))] <- "spatial_offset"
        # ofs <- get_ngbvalue(coords = df[,c('x','y')],
        #                     env =  model$offset,
        #                     longlat = terra::is.lonlat(bg),
        #                     field_space = c('x','y')
        #                     )
        model$biodiversity[[1]]$offset <- ofs
      }
      # Define expectation as very small vector following Renner et al.
      w <- ppm_weights(df = df,
                       pa = model$biodiversity[[1]]$observations[['observed']],
                       bg = bg,
                       weight = 1e-6
      )
      df$w <- w * (1/model$biodiversity[[1]]$expect) # Also add as column

      model$biodiversity[[1]]$predictors <- df
      model$biodiversity[[1]]$expect <- df$w

      # Rasterize observed presences
      pres <- terra::rasterize( guess_sf( model$biodiversity[[1]]$observations[,c("x","y")] ),
                                bg, fun = 'count', background = 0)
      # Get for the full dataset
      w_full <- ppm_weights(df = model$predictors,
                            pa = pres[],
                            bg = bg,
                            weight = 1 # Set those to 1 so that absences become ratio of pres/abs
      )

      # Add exposure to full model predictor
      model$exposure <- w_full * (1/unique(model$biodiversity[[1]]$expect)[1])

    } else if(model$biodiversity[[1]]$family != 'poisson'){
      # calculating the case weights (equal weights)
      # the order of weights should be the same as presences and backgrounds in the training data
      prNum <- as.numeric(table(model$biodiversity[[1]]$observations[['observed']])["1"]) # number of presences
      bgNum <- as.numeric(table(model$biodiversity[[1]]$observations[['observed']])["0"]) # number of backgrounds
      w <- ifelse(model$biodiversity[[1]]$observations[['observed']] == 1, 1, prNum / bgNum)
      # If family is not poisson, assume factor distribution for response
      model$biodiversity[[1]]$observations[['observed']] <- factor(model$biodiversity[[1]]$observations[['observed']])

      # Add offset if existent
      if(!is.Waiver(model$offset)){
        ofs <- get_rastervalue(coords = model$biodiversity[[1]]$observations[,c('x','y')],
                               env = model$offset_object,
                               rm.na = FALSE)
        # Rename to spatial offset
        names(ofs)[which(names(ofs)==names(model$offset_object))] <- "spatial_offset"
        model$biodiversity[[1]]$offset <- ofs
      }

      model$biodiversity[[1]]$expect <- w * model$biodiversity[[1]]$expect
    }

    # ---- #
    # Detect and format the family
    fam <- model$biodiversity[[1]]$family
    li <- model$biodiversity[[1]]$link
    if(is.null(li) && fam == "binomial") li <- "logit"
    fam <- switch (fam,
                   "poisson" = mboost::Poisson(),
                   "binomial" = mboost::Binomial(type = "glm", link = li),
                   "gaussian" = mboost::Gaussian(),
                   "hurdle" = mboost::Hurdle(nuirange = c(0,100))
    )
    self$data[['family']] <- fam
    assertthat::assert_that(inherits(fam,'boost_family'),msg = 'Family misspecified.')

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
    if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','green',paste0('Starting fitting: ', name))

    # Get output raster
    prediction <- self$get_data('template')
    # Get boosting control and family data
    bc <- self$get_data('bc')
    bc$trace <- settings$get('verbose') # Reset trace as specified
    params <- self$get_data('params')
    fam <- self$get_data('family')

    # All other needed data for model fitting
    equation <- model$biodiversity[[1]]$equation
    data <- cbind(model$biodiversity[[1]]$predictors,
                  data.frame(observed = model$biodiversity[[1]]$observations[,'observed', drop = TRUE]) )
    w <- model$biodiversity[[1]]$expect

    # Select predictors
    full <- model$predictors
    full <- subset(full, select = c('x','y',model$biodiversity[[1]]$predictors_names))
    full$cellid <- rownames(full) # Add row.names
    full$w <- model$exposure
    full$Intercept <- 1
    full <- subset(full, stats::complete.cases(full))

    # Clamp?
    if( settings$get("clamp") ) full <- clamp_predictions(model, full)

    # Rescale exposure
    check_package('scales')
    w <- scales::rescale(w, to = c(1e-6, 1))
    full$w <- scales::rescale(full$w, to = c(1e-6, 1))
    if(anyNA(w)){
      w[is.na(w)] <- 1e-6
      full$w[is.na(full$w)] <- 1e-6
    }

    assertthat::assert_that(
      is.null(w) || length(w) == nrow(data),
      is.formula(equation),
      all(model$biodiversity[[1]]$predictors_names %in% names(full)),
      all(names(full[,model$biodiversity[[1]]$predictors_names]) %in% names(data)),
      all( model$biodiversity[[1]]$predictors_names %in% names(full) )
    )

    if(!is.Waiver(model$offset)){
      # Add offset to full prediction and load vector
      n <- data.frame(model$offset[as.numeric(full$cellid), "spatial_offset"], model$offset[as.numeric(full$cellid), "spatial_offset"] )
      names(n) <- c( "spatial_offset", paste0('offset(',"spatial_offset",')') )
      # Add weights
      # n <- n + full$w
      full <- cbind(full, n)
      # And for biodiversity object
      n <- cbind(model$biodiversity[[1]]$offset[,"spatial_offset"],
                 model$biodiversity[[1]]$offset[,"spatial_offset"]) |> as.data.frame()
      names(n) <- c( "spatial_offset", paste0('offset(',"spatial_offset",')') )
      # Add weights
      # n <- n + w
      data <- cbind(data, n)
    }

    # --- #
    # Fit the base model
    # First test for infinitives
    bc2 <- bc;
    bc2$mstop <- 1; bc2$trace <- FALSE; bc2$stopintern <- TRUE
    test <- try({
      mboost::gamboost(
        formula = equation,
        data = data,
        # weights = w,
        family = fam,
        offset = w, # Add exposure as offset
        control = bc2
      )
    },silent = FALSE)
    if(any(is.infinite(test$risk()))){ return("Infinite residuals. Try simplifying model")}
    fit_gdb <- try({
      mboost::gamboost(
        formula = equation,
        data = data,
        # weights = w,
        family = fam,
        offset = w, # Add exposure as offset
        control = bc
      )
    },silent = FALSE)
    # If error, decrease step size by a factor of 10 and try again.
    if(inherits(fit_gdb, "try-error") || length(names(stats::coef(fit_gdb)))< 2){
      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','red','Reducing learning rate by 1/100.')
      bc$nu <- bc$nu * 0.01
      fit_gdb <- try({
        mboost::gamboost(
          formula = equation,
          data = data,
          # weights = w,
          family = fam,
          offset = w, # Add exposure as offset
          control = bc
        )
      },silent = FALSE)
      if(inherits(fit_gdb, "try-error")) {
        myLog('[Estimation]','red','Fitting failed. Check model and alter parameters!')
        return(NULL)
      }
    }

    if(settings$get('optim_hyperparam')){

      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','green','Starting parameter search for optimal stopping.')
      # 5 fold Cross validation to prevent overfitting
      if(getOption("ibis.runparallel")){
        grs <- seq(from = 10, to = max( bc$mstop *5), by = 10)
        cvf <- mboost::cv(stats::model.weights(fit_gdb),B = 5, type = "kfold")

        # Start cluster
        # cl <- parallel::makeCluster( getOption('ibis.nthread') )
        try({cvm <- mboost::cvrisk(fit_gdb,
                                   folds = cvf, grid = grs,
                                   papply = parallel::mclapply,
                                   mc.cores = getOption("ibis.nthread"))
        }, silent = TRUE)
        # parallel::stopCluster(cl)
        rm(cvf, grs)

      } else {
        grs <- seq(from = 10, to = max( bc$mstop *5), by = 10)
        cvf <- mboost::cv(stats::model.weights(fit_gdb),B = 5, type = "kfold")
        try({cvm <- mboost::cvrisk(fit_gdb,
                                   folds = cvf, grid = grs,
                                   papply = parallel::parLapply )
        }, silent = TRUE)
        rm(cvf, grs)
      }
      # Check whether crossvalidation has run through successfully
      if(exists('cvm') && mboost::mstop(cvm) > 0){
        # Set the model to the optimal mstop to limit overfitting
        fit_gdb[mboost::mstop(cvm)]
      } else {cvm <- new_waiver()}
    } else {
      cvm <- new_waiver()
    }

    # Predict spatially
    if(!settings$get('inference_only')){
      # Messager
      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','green','Starting prediction...')

      # Set target variables to bias_value for prediction if specified
      if(!is.Waiver(settings$get('bias_variable'))){
        for(i in 1:length(settings$get('bias_variable'))){
          if(settings$get('bias_variable')[i] %notin% names(full)){
            if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','red','Did not find bias variable in prediction object!')
            next()
          }
          full[[settings$get('bias_variable')[i]]] <- settings$get('bias_value')[i]
        }
      }

      if(getOption("ibis.runparallel",default = FALSE)){
        check_package("doFuture")
        if(!("doFuture" %in% loadedNamespaces()) || ('doFuture' %notin% utils::sessionInfo()$otherPkgs) ) {
          try({requireNamespace('doFuture');attachNamespace("doFuture")},silent = TRUE)
        }

        # Chunk the data
        splits <- chunk_data(full, N = getOption("ibis.nthread",default = 10), index_only = TRUE)

        pred_gdb <- foreach::foreach(s = splits,
                                      .combine = "rbind",
                                      .inorder = TRUE,
                                      .options.future = list(seed = TRUE,
                                                             packages = c("mboost"))
        ) %dofuture% {
          # Make a prediction
          suppressWarnings(
            mboost::predict.mboost(object = fit_gdb,
                                   newdata = full[s,],
                                   type = self$get_data('params')$type,
                                   aggregate = 'sum',
                                   offset = full$w[s])
          )
        }

      } else {
        # Make a prediction
        suppressWarnings(
          pred_gdb <- mboost::predict.mboost(object = fit_gdb, newdata = full,
                                             type = self$get_data('params')$type,
                                             aggregate = 'sum',
                                             offset = full$w)
        )
      }

      # Fill output
      prediction[as.numeric(full$cellid)] <- pred_gdb[,1]
      names(prediction) <- 'mean'
      rm(pred_gdb)
    } else {
      prediction <- NULL
    }

    # Compute end of computation time
    settings$set('end.time', Sys.time())
    # Also append boosting control option to settings
    for(entry in names(params)) settings$set(entry, params[[entry]])
    for(entry in names(bc)) settings$set(entry, bc[[entry]])

    # Definition of GDB Model object ----
    obj <- DistributionModel # Make a copy to set new functions

    # Project function
    obj$set("public", "project", function(newdata, type = NULL, layer = "mean"){
      assertthat::assert_that('fit_best' %in% names(self$fits),
                              is.data.frame(newdata) || is.matrix(newdata),
                              assertthat::has_name(newdata,c('x','y'))
      )
      # Get model
      mod <- self$get_data('fit_best')
      model <- self$model
      settings <- self$settings
      assertthat::assert_that(inherits(mod,'mboost'),msg = 'No model found!')
      if(is.null(type)) type <- settings$get('type')
      # Check that all variables are in provided data.frame
      assertthat::assert_that(all( as.character(mboost::extract(mod,'variable.names')) %in% names(newdata) ))

      # Clamp?
      if( settings$get("clamp") ) newdata <- clamp_predictions(model, newdata)

      # Set target variables to bias_value for prediction if specified
      if(!is.Waiver(settings$get('bias_variable'))){
        for(i in 1:length(settings$get('bias_variable'))){
          if(settings$get('bias_variable')[i] %notin% names(newdata)){
            if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','red','Did not find bias variable in prediction object!')
            next()
          }
          newdata[[settings$get('bias_variable')[i]]] <- settings$get('bias_value')[i]
        }
      }

      # Add rowid
      newdata$rowid <- 1:nrow(newdata)
      # Subset to non-missing data
      newdata_sub <- subset(newdata, stats::complete.cases(newdata))

      if(getOption("ibis.runparallel",default = FALSE)){
        check_package("doFuture")
        if(!("doFuture" %in% loadedNamespaces()) || ('doFuture' %notin% utils::sessionInfo()$otherPkgs) ) {
          try({requireNamespace('doFuture');attachNamespace("doFuture")},silent = TRUE)
        }
        # Chunk the data
        splits <- chunk_data(newdata_sub, N = getOption("ibis.nthread",default = 10),
                             index_only = TRUE)

        y <- foreach::foreach(s = splits,
                              .combine = "rbind",
                              .inorder = TRUE,
                              .options.future = list(seed = TRUE,
                                                     packages = c("mboost"))
        ) %dofuture% {
          # Make a prediction
          suppressWarnings(
            mboost::predict.mboost(object = mod,
                                   newdata = newdata_sub[s,],
                                   type = type,
                                   aggregate = 'sum')
          )
        }
      } else {
        # Predict
        y <- suppressWarnings(
          mboost::predict.mboost(object = mod, newdata = newdata_sub,
                                 type = type, aggregate = 'sum')
        )
      }

      # Make empty template
      if(nrow(newdata)==nrow(model$predictors)){
        prediction <- try({model_to_background(model)}, silent = TRUE)
        prediction[as.numeric(newdata_sub$rowid)] <- y[,1]
      } else {
        assertthat::assert_that(utils::hasName(newdata_sub,"x")&&utils::hasName(newdata_sub,"y"),
                                msg = "Projection data.frame has no valid coordinates or differs in grain!")
        prediction <- try({
          terra::rast(newdata_sub[,c("x", "y")],
                      crs = terra::crs(model$background),
                      type = "xyz") |>
            emptyraster()
        }, silent = TRUE)
        prediction[as.numeric(newdata_sub$rowid)] <- y[,1]
      }
      names(prediction) <- "mean" # Rename to mean, layer parameter gets ignored for this engine
      return(prediction)
    },overwrite = TRUE)

    # Partial effect
    obj$set("public", "partial", function(x.var = NULL, constant = NULL, variable_length = 100, values = NULL,
                           newdata = NULL, plot = FALSE, type = NULL){
      # Assert that variable(s) are in fitted model
      assertthat::assert_that(is.character(x.var) || is.null(x.var),
                              inherits(self$get_data('fit_best'), 'mboost'),
                              is.numeric(variable_length),
                              is.null(newdata) || is.data.frame(newdata))

      # Unlike the effects function, build specific predictor for target variable(s) only
      variables <- mboost::extract(self$get_data('fit_best'),'variable.names')

      # Match x.var to argument
      if(is.null(x.var)) {
        x.var <- variables
      } else {
        x.var <- match.arg(x.var, variables, several.ok = TRUE)
      }

      # extract settings and model
      settings <- self$settings
      if(is.null(type)) type <- settings$get('type')
      model <- self$model

      # create dummy data
      if(is.null(newdata)) {
        dummy <- base::as.data.frame(matrix(nrow = variable_length, ncol = length(variables)))
        names(dummy) <- variables

        if(is.null(constant)){
          # FIXME: What about factors?
          nn <- model$predictors_types$predictors[which(model$predictors_types$type=='numeric')]
          constant <- as.list(apply(model$predictors[, nn], 2, function(x) mean(x, na.rm = TRUE)))
          dummy[, nn] <- constant[nn]
        } else {
          dummy[, variables] <- constant
        }

      } else {
        # Assume all present
        dummy <- dplyr::select(newdata, dplyr::any_of(as.character(variables)))
      }

      # create out list
      out <- vector(mode = "list", length = length(x.var))
      names(out) <- x.var

      # calc effect for each x.var
      for(v in x.var){

        dummy_temp <- dummy

        # create partial effect range
        if(v %in% model$predictors_types$predictors[model$predictors_types$type=="factor"]){
          range_temp <- levels(model$predictors[,v])
          dummy_temp[, v] <- factor(range_temp)
        } else {
          if(!is.null(values)){
            variable_length <- length(values)
            assertthat::assert_that(length(values) >= 1)
            dummy_temp[, v] <- values
          } else {
            range_temp <- range(model$predictors[[v]], na.rm = TRUE)
            dummy_temp[, v] <- seq(range_temp[1], range_temp[2], length.out = variable_length)
          }
        }

        # Now predict with model
        suppressWarnings(pp <- mboost::predict.mboost(object = self$get_data('fit_best'),
                                                      newdata = dummy_temp, which = v,
                                                      type = type, aggregate = 'sum'))

        # If bbs is present and non-linear, use bbs estimate. If model is fitted
        # to non-linear then always the linear (bols, fitted by default) and smooth (bbs)
        # is stored in bbs. In most cases the linear effect is stable and regularized out
        # though.
        # Propose to always take non-linear when found in settings
        # MH: Especially the bbs column is named slightly different depending on knots
        # MH: Thus, only searching for bbs(x) or bols(x)
        # MH: IF prior is present, we need to extract bmono
        if(!settings$get("only_linear")){
          # Combine with
          out[[v]] <- data.frame(variable = v, partial_effect = dummy_temp[, v],
                                 mean = pp[, grep("bbs|bmono", colnames(pp))]) #pp[,grep(paste0("bbs\\(", v,"\\)"), colnames(pp))
        } else {
          # Combine with
          out[[v]] <- data.frame(variable = v, partial_effect = dummy_temp[, v],
                                 mean = pp[, grep("bols|bmono", colnames(pp))]) #pp[,grep(v, colnames(pp))]
        }
      }

      # bind to one data.frame
      out <- do.call(what = rbind, args = c(out, make.row.names = FALSE))

      # If plot, make plot, otherwise
      if(plot){
        g <- ggplot2::ggplot(data = out, ggplot2::aes(x = partial_effect)) +
          ggplot2::theme_classic() +
          ggplot2::geom_line(ggplot2::aes(y = mean)) +
          ggplot2::facet_wrap(. ~ variable, scales = "free") +
          ggplot2::labs(x = "Variable",
                        y = paste0("Partial effect (",ifelse(settings$get("only_linear"),"linear","smooth"),")"))
        print(g)
      }
      return(out)
    },overwrite = TRUE)

    # Spatial partial effect plots
    obj$set("public", "spartial", function(x.var, constant = NULL, newdata = NULL,
                                           plot = TRUE, type = NULL, ...){
      assertthat::assert_that('fit_best' %in% names(self$fits),
                              is.character(x.var), length(x.var) == 1)
      # Get model and make empty template
      mod <- self$get_data('fit_best')
      model <- self$model
      settings <- self$settings
      # Also check that what is present in coefficients of model
      variables <- as.character( mboost::extract(mod,'variable.names') )
      assertthat::assert_that(x.var %in% variables,
                              msg = "Variable not found in model! Regularized out?" )

      if(is.null(type)) type <- self$settings$get("type")
      type <- match.arg(type, c("link", "response"), several.ok = FALSE)
      settings$set("type", type)

      # Make template of target variable(s)
      template <- model_to_background(model)

      # Get target variables and predict
      if(!is.null(newdata)){
        df <- newdata
      } else {
        df <- model$predictors
      }
      assertthat::assert_that(x.var %in% colnames(df),
                              msg = "Variable not found in provided data.")

      # Set all variables other the target variable to constant
      if(is.null(constant)){
        # Calculate mean
        nn <- model$predictors_types$predictors[which(model$predictors_types$type=='numeric')]
        constant <- apply(df[,nn], 2, function(x) mean(x, na.rm=T))
        for(v in variables[ variables %notin% x.var]){
          if(v %notin% colnames(df) ) next()
          df[!is.na(df[v]),v] <- constant[v]
        }
      } else {
        df[!is.na(df[,x.var]), variables] <- constant
      }
      df$rowid <- as.numeric( rownames(df) )
      assertthat::assert_that(nrow(df)==ncell(template))

      pp <- suppressWarnings(
        mboost::predict.mboost(mod, newdata = df, which = x.var,
                               type = settings$get('type'), aggregate = 'sum')
      )
      # If both linear and smooth effects are in model
      if(length(df$rowid[which(!is.na(df[[x.var]]))] ) == length(pp[,ncol(pp)])){
        template[ df$rowid[which(!is.na(df[[x.var]]))] ] <- pp[,ncol(pp)]
      } else { template[] <- pp[, ncol(pp) ]}
      names(template) <- paste0('partial__',x.var)

      if(plot){
        # Plot both partial spatial partial
        par.ori <- graphics::par(no.readonly = TRUE)
        graphics::par(mfrow = c(1,3))
        terra::plot(template, main = expression(f[partial]),
                    col = ibis_colours$ohsu_palette)
        mboost::plot.mboost(mod,which = x.var)
        graphics::par(par.ori)
      }
      return(template)
    },overwrite = TRUE)

    # Model convergence check
    obj$set("public", "has_converged", function(){
      fit <- self$get_data("fit_best")
      if(is.Waiver(fit)) return(FALSE)
      # Get risks
      evl <- fit$risk()
      if(fit$mstop() == length(evl)) return(FALSE)
      return(TRUE)
    },overwrite = TRUE)

    # Residual function
    obj$set("public", "get_residuals", function(){
      # Get best object
      obj <- self$get_data("fit_best")
      if(is.Waiver(obj)) return(obj)
      # Get residuals
      rd <- obj$resid()
      assertthat::assert_that(length(rd)>0)
      return(rd)
    },overwrite = TRUE)

    # Get coefficients
    obj$set("public", "get_coefficients", function(){
      # Returns a vector of the coefficients with direction/importance
      cofs <- self$summary()
      if(nrow(cofs)==0) return(NULL)
      # Sanitize and remove base learners from object
      cofs$variable <- gsub('bols\\(|bbs\\(|bmono\\(', '', cofs$variable)
      cofs$variable <- sapply(strsplit(cofs$variable, ","), function(z) z[[1]])
      cofs$variable <- gsub('\\)', '', cofs$variable)
      cofs <- cofs |> dplyr::select(variable, beta)
      names(cofs) <- c("Feature", "Beta")
      return(cofs)
    },overwrite = TRUE)

    # Spatial latent effect
    obj$set("public", "plot_spatial", function(plot = TRUE){
      assertthat::assert_that('fit_best' %in% names(self$fits) )
      # Get model and make empty template
      mod <- self$get_data('fit_best')
      model <- self$model

      # Also check that what is present in coefficients of model
      vars <- as.character( mboost::extract(mod,'bnames') )
      assertthat::assert_that(length(grep('bspatial',vars))>0,
                              msg = 'No spatial effect found in model!')

      # Make template of target variable(s)
      template <- model_to_background(model)

      # Get target variables and predict
      target <- self$model$predictors[,c('x','y')]

      y <- suppressWarnings(
        mboost::predict.mboost(mod, newdata = target, which = c('x','y'))
      )
      assertthat::assert_that(nrow(target)==nrow(y))
      template[] <- y[,2]
      names(template) <- paste0('partial__','space')
      # Mask with background
      template <- terra::mask(template, model$background )

      # Plot both partial spatial partial
      if(plot){
        terra::plot(template, main = expression(f[partial]), col = ibis_colours$divg_bluegreen )
      }
      return(template)

    },overwrite = TRUE)

    out <- obj$new(name = "GDB-Model")
    out$id <- model$id
    out$model <- model
    out$settings <- settings
    out$fits = list(
      "fit_best" = fit_gdb,
      "fit_cv" = cvm,
      "fit_best_equation" = equation,
      "prediction" = prediction
    )

    return(out)
  },overwrite = TRUE)

  # Set engine in distribution object
  eg <- eg$new(engine =  "GDB-Engine", name = "<GDB>")
  eg$data <- list(
    'template' = template,
    'bc' = bc,
    'params' = params
  )

  y <- x$clone(deep = TRUE)
  return( y$set_engine(eg) )
} # End of function
