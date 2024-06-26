#' @include class-engine.R class-distributionmodel.R
NULL

#' Engine for regularized regression models
#'
#' @description This engine allows the estimation of linear coefficients using
#' either ridge, lasso or elastic net regressions techniques. Backbone of this
#' engine is the \pkg{glmnet} R-package which is commonly used in SDMs,
#' including the popular \code{'maxnet'} (e.g. Maxent) package. Ultimately this
#' engine is an equivalent of [engine_breg], but in a "frequentist" setting. If
#' user aim to emulate a model that most closely resembles maxent within the
#' ibis.iSDM modelling framework, then this package is the best way of doing so.
#' Compared to the \code{'maxnet'} R-package, a number of efficiency settings
#' are implemented in particular for cross-validation of alpha and lambda
#' values.
#'
#' Limited amount of prior information can be specified for this engine,
#' specifically via offsets or as [`GLMNETPrior`], which allow to specify priors
#' as regularization constants.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param alpha A [`numeric`] giving the elasticnet mixing parameter, which has
#' to be between \code{0} and \code{1}. \code{alpha=1} is the lasso penalty,
#' and \code{alpha=0} the ridge penalty (Default: \code{0}).
#' @param nlambda A [`numeric`] giving the number of lambda values to be used (Default: \code{100}).
#' @param lambda A [`numeric`] with a user supplied estimate of lambda. Usually
#' best to let this parameter be determined deterministically (Default: \code{NULL}).
#' @param type The mode used for creating posterior predictions. Either making
#' \code{"link"} or \code{"response"} (Default: \code{"response"}).
#' @param ... Other parameters passed on to glmnet.
#'
#' @details Regularized regressions are effectively GLMs that are fitted with
#' ridge, lasso or elastic-net regularization. Which of them is chosen is
#' critical dependent on the alpha value: [*] For \code{alpha} equal to \code{0}
#' a ridge regularization is used. Ridge regularization has the property that it
#' doesn't remove variables entirely, but instead sets their coefficients to
#' \code{0}. [*] For \code{alpha} equal to \code{1} a lasso regularization is
#' used. Lassos tend to remove those coefficients fully from the final model
#' that do not improve the loss function. [*] For \code{alpha} values between
#' \code{0} and \code{1} a elastic-net regularization is used, which is
#' essentially a combination of the two. The optimal lambda parameter can be
#' determined via cross-validation. For this option set \code{"varsel"} in
#' `train()` to \code{"reg"}.
#'
#' @returns An [Engine].
#'
#' @references
#' * Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010). Regularization Paths
#' for Generalized Linear Models via Coordinate Descent. Journal of Statistical Software,
#' 33(1), 1-22. URL https://www.jstatsoft.org/v33/i01/.
#' * Renner, I.W., Elith, J., Baddeley, A., Fithian, W., Hastie, T., Phillips, S.J.,
#'  Popovic, G. and Warton, D.I., 2015. Point process models for presence‐only analysis.
#'   Methods in Ecology and Evolution, 6(4), pp.366-379.
#' * Fithian, W. & Hastie, T. (2013) Finite-sample equivalence in statistical models
#' for presence-only data. The Annals of Applied Statistics 7, 1917–1939
#'
#' @family engine
#'
#' @examples
#' \dontrun{
#' # Add GLMNET as an engine
#' x <- distribution(background) |> engine_glmnet(iter = 1000)
#' }
#'
#' @name engine_glmnet
NULL

#' @rdname engine_glmnet
#' @export
engine_glmnet <- function(x,
                          alpha = 0,
                          nlambda = 100,
                          lambda = NULL,
                          type = "response",
                          ...) {

  # Check whether glmnet package is available
  check_package('glmnet')
  if(!("glmnet" %in% loadedNamespaces()) || ('glmnet' %notin% utils::sessionInfo()$otherPkgs) ) {
    try({requireNamespace('glmnet');attachNamespace("glmnet")},silent = TRUE)
  }
  check_package('glmnetUtils') # glmnetUtils is a helper functions for formulas
  if(!("glmnetUtils" %in% loadedNamespaces()) || ('glmnetUtils' %notin% utils::sessionInfo()$otherPkgs) ) {
    try({requireNamespace('glmnetUtils');attachNamespace("glmnetUtils")},silent = TRUE)
  }

  # assert that arguments are valid
  assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                          inherits(x$background,'sf'),
                          is.numeric(alpha),
                          is.numeric(lambda) || is.null(lambda),
                          is.numeric(nlambda),
                          is.character(type)
  )
  assertthat::assert_that(alpha>=0, alpha <=1)
  type <- match.arg(type, choices = c("predictor","link", "response"),several.ok = FALSE)
  if(type=="predictor") type <- "link" # Convenience conversion

  # Create a background raster
  if(is.Waiver(x$predictors)){
    # Create from background
    template <- terra::rast(
      ext = terra::ext(x$background),
      crs = terra::crs(x$background),
      res = c(diff( (sf::st_bbox(x$background)[c(1,3)]) ) / 100, # Simplified assumption for resolution
              diff( (sf::st_bbox(x$background)[c(1,3)]) ) / 100
      )
    )
  } else {
    # If predictor existing, use them
    template <- emptyraster(x$predictors$get_data() )
  }

  # Burn in the background
  template <- terra::rasterize(x$background, template, field = 0)

  # Set up the parameter list
  params <- list(
    alpha = alpha,
    lambda = lambda,
    nlambda = nlambda,
    type = type,
    ...
  )

  # Print a message in case there is already an engine object
  if(!is.Waiver(x$engine)) myLog('[Setup]','yellow','Replacing currently selected engine.')

  # Define new engine object of class
  eg <- Engine

  # Dummy function for spatial latent effects
  eg$set("public", "calc_latent_spatial", function(type = NULL, priors = NULL){
    new_waiver()
  },overwrite = TRUE)

  # Dummy function for getting the equation of latent effects
  eg$set("public", "get_equation_latent_spatial", function(method){
    new_waiver()
  },overwrite = TRUE)

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
    settings$set('type', params$type)

    # Distribution specific procedure
    fam <- model$biodiversity[[1]]$family
    form <- model$biodiversity[[1]]$equation

    # # -- #
    # # Respecify the predictor names if not matching
    # te <- formula_terms(form)
    # if(length(te) != nrow(model$biodiversity[[1]]$predictors_types)){
    #   model$biodiversity[[1]]$predictors_names <-
    #     model$biodiversity[[1]]$predictors_names[model$biodiversity[[1]]$predictors_names %in% te]
    #   model$biodiversity[[1]]$predictors_types <-
    #     model$biodiversity[[1]]$predictors_types |> dplyr::filter(
    #     predictors %in% te)
    # }
    # # -- #

    # If a poisson family is used, weight the observations by their exposure
    if(fam == "poisson"){
      # Get background layer
      bg <- self$get_data("template")
      assertthat::assert_that(!is.na( terra::global(bg, "min", na.rm = TRUE)[,1]))

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

      # Overwrite observation data
      model$biodiversity[[1]]$observations <- presabs

      # Preprocessing security checks
      assertthat::assert_that( all( model$biodiversity[[1]]$observations[['observed']] >= 0 ),
                               any(!is.na(presabs[['observed']])),
                               length(model$biodiversity[[1]]$expect) == nrow(model$biodiversity[[1]]$observations),
                               nrow(df) == nrow(model$biodiversity[[1]]$observations)
      )

      # Add offset if existent
      if(!is.Waiver(model$offset)){
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
                       weight = 1e-6, # Arbitrary small weight
                       type = "DWPR" # Weights for down-weighted Poisson regression
      )
      assertthat::assert_that(length(w) == nrow(df))

      model$biodiversity[[1]]$predictors <- df
      model$biodiversity[[1]]$expect <- w * (1/model$biodiversity[[1]]$expect) # Multiply with prior weight

      # Rasterize observed presences
      pres <- terra::rasterize( guess_sf(model$biodiversity[[1]]$observations[,c("x","y")]),
                                bg, fun = 'count', background = 0)
      # Get for the full dataset
      w_full <- ppm_weights(df = model$predictors,
                            pa = pres[],
                            bg = bg,
                            weight = 1 # Set those to 1 so that absences become ratio of pres/abs
      )

      # Add exposure to full model predictor
      model$exposure <- w_full * (1/unique(model$biodiversity[[1]]$expect)[1]) # Multiply with prior weight (first value)

    } else if(fam == "binomial"){
      # Check that observations are all <=1
      model$biodiversity[[1]]$observations[['observed']] <- ifelse(model$biodiversity[[1]]$observations[['observed']]>=1,1,0)
      # calculating the case weights (equal weights)
      # the order of weights should be the same as presences and backgrounds in the training data
      prNum <- as.numeric(table(model$biodiversity[[1]]$observations[['observed']])["1"]) # number of presences
      bgNum <- as.numeric(table(model$biodiversity[[1]]$observations[['observed']])["0"]) # number of backgrounds
      w <- ifelse(model$biodiversity[[1]]$observations[['observed']] == 1, 1, prNum / bgNum)
      model$biodiversity[[1]]$expect <- w * model$biodiversity[[1]]$expect
      model$biodiversity[[1]]$observations$observed <- as.factor( model$biodiversity[[1]]$observations$observed )
    }

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

    # Verbosity
    verbose <- settings$get("verbose")

    # Set prediction type also for later
    settings$set('type', self$get_data("params")$type)

    # seed
    seed <- settings$get("seed")
    if(is.Waiver(seed)) { settings$set('seed', getOption("ibis.seed", default = 1000)) }

    # Get output raster
    prediction <- self$get_data('template')

    # Get parameters control
    params <- self$get_data('params')
    # Set glmnet control
    glmnet::glmnet.control(factory = TRUE) # Reset to default

    # All other needed data for model fitting
    fam <- model$biodiversity[[1]]$family
    li <- model$biodiversity[[1]]$link
    if(!is.null(li)){
      if(li %in% c("cloglog", "logit", "probit")){
        fam <- stats::binomial(link = li)
      } else {
        if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','red',paste0("Custom link functions not supported!"))
      }
    }

    form <- model$biodiversity[[1]]$equation
    df <- cbind(model$biodiversity[[1]]$predictors,
                data.frame(observed = model$biodiversity[[1]]$observations[,'observed', drop = TRUE])
    )
    df <- subset(df, select = c(model$biodiversity[[1]]$predictors_names, "observed"))
    w <- df$w <- model$biodiversity[[1]]$expect # The expected exposure
    # Get full prediction container
    full <- model$predictors
    w_full <- model$exposure

    # Subset the predictor types to only those present
    te <- formula_terms(form)
    model$biodiversity[[1]]$predictors_types <-
      model$biodiversity[[1]]$predictors_types |> dplyr::filter(predictors %in% te)
    model$biodiversity[[1]]$predictors_names <-  intersect(model$biodiversity[[1]]$predictors_names, te)

    # Get offset and add it to exposure
    if(!is.Waiver(model$offset)){
      # Add offset to full prediction and load vector
      ofs <- model$biodiversity[[1]]$offset[, 'spatial_offset']
      ofs_pred <- model$offset[,'spatial_offset']
    } else { ofs <- NULL; ofs_pred <- NULL }

    # Check priors, e.g penalty factors
    p.fac <- rep(1, sum( model$biodiversity[[1]]$predictors_types$type=="numeric") ) # Added plus 1 for the weight?
    names(p.fac) <- model$biodiversity[[1]]$predictors_names[which(model$biodiversity[[1]]$predictors_types$type=="numeric")]
    # Then add each factor level if set
    if(any(model$predictors_types$type=="factor")){
      fac <- model$biodiversity[[1]]$predictors_names[which(model$biodiversity[[1]]$predictors_types$type=="factor")]
      # return penalty factor for each level of each factor (even if level values are identical across factors)
      p.fac <- c(p.fac, rep(1, sum(apply(df[, fac, drop = FALSE], 2, function(x) length(unique(x))))))
    }
    # Duplicate p.fac container for lower and upper limits
    lowlim <- rep(-Inf, length(p.fac)) |> stats::setNames(names(p.fac))
    upplim <- rep(Inf, length(p.fac)) |> stats::setNames(names(p.fac))

    # Trick for creation for some default lambda values for the regularization multiplier
    if(is.null(params$lambda)){
      reg <- default.regularization(p = df$observed, m = stats::model.matrix(form, df)) * c(1, p.fac) # add 1 for the intercept
      params$lambda <- 10^(seq(4, 0, length.out = 200)) * sum(p.fac)/length(p.fac) * sum(p.fac)/sum(w)
      if(anyNA(params$lambda)) params$lambda <- NULL
    }

    if(!is.Waiver(model$priors)){
      assertthat::assert_that(all(model$priors$varnames() %in% model$predictors_names))
      # Reset those contained in the prior object
      for(v in model$priors$varnames()){
        if(!(v %in% names(p.fac))) next()
        p.fac[v]  <- model$priors$get(v, what = "value")
        lowlim[v] <- model$priors$get(v, what = "lims")[1]
        upplim[v] <- model$priors$get(v, what = "lims")[2]
      }
    }

    # Clamp?
    if( settings$get("clamp") ) full <- clamp_predictions(model, full)

    # -- #
    # Expand predictors if non-linear is specified in settings
    if(settings$get('only_linear') == FALSE){
      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','yellow',
                                                'Non-linearity to glmnet is best introduced by adding derivates. Ignored!')
      # linear_predictors <- attr(stats::terms.formula(form), "term.labels")
      # m <- outer(linear_predictors, linear_predictors, function(x, y) paste(x, y, sep = ":"))
      #
      # form <- stats::update.formula(form,
      #                        paste0(
      #                          ". ~ . +",
      #                          paste0("I(", linear_predictors,"^2)",collapse = " + "),
      #                          " + ",
      #                          paste0(m[lower.tri(m)], collapse = " + ")
      #                        ))
      # # Update penalty factors and limits
      # for(var in attr(stats::terms.formula(form), "term.labels")){
      #   if(!(var %in% p.fac)){
      #     v <- 1 # Take the maximum regularization penalty by default
      #     vlow <- -Inf; vupp <- Inf
      #     names(v) <- var; names(vlow) <- var; names(vupp) <- var
      #     p.fac <- append(p.fac, v)
      #     lowlim <- append(lowlim, vlow); upplim <- append(upplim, vupp)
      #   }
      # }
    }

    assertthat::assert_that(
      is.null(w) || length(w) == nrow(df),
      is.null(ofs) || is.vector(ofs),
      is.null(ofs_pred) || is.vector(ofs_pred),
      length(p.fac) == length(lowlim),
      all(w >= 0,na.rm = TRUE)
    )

    # --- #
    # Determine the optimal lambda through k-fold cross-validation
    if(getOption("ibis.runparallel")){
      if(!foreach::getDoParRegistered()) ibis_future(cores = getOption("ibis.nthread"),
                                                     strategy = getOption("ibis.futurestrategy"))
    }
    # Depending if regularized should be set, specify this separately
    if( (settings$get('optim_hyperparam')) ){
      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','green',
                                                'Finding optimal hyper parameters alpha and lambda.')
      cv_gn <- try({
        glmnetUtils::cva.glmnet(formula = form,
                                data = df,
                                alpha = seq(0,1,.1), # Elastic net parameters to test
                                lambda = params$lambda, # Overwrite lambda
                                weights = w, # Case weights
                                offset = ofs,
                                family = fam,
                                penalty.factor = p.fac,
                                # Option for limiting the coefficients
                                lower.limits = lowlim,
                                upper.limits = upplim,
                                standardize = FALSE, # Don't standardize to avoid doing anything to weights
                                maxit = (10^5)*2, # Increase the maximum number of passes for lambda
                                parallel = getOption("ibis.runparallel"),
                                trace.it = settings$get("verbose"),
                                nfolds = 10  # number of folds for cross-validation
        )
      },silent = FALSE)
    } else {
      cv_gn <- try({
        glmnetUtils::cv.glmnet(formula = form,
                               data = df,
                               alpha = params$alpha, # Elastic net mixing parameter
                               lambda = params$lambda, # Overwrite lambda
                               weights = w, # Case weights
                               offset = ofs,
                               family = fam,
                               penalty.factor = p.fac,
                               # Option for limiting the coefficients
                               lower.limits = lowlim,
                               upper.limits = upplim,
                               standardize = FALSE, # Don't standardize to avoid doing anything to weights
                               maxit = (10^5)*2, # Increase the maximum number of passes for lambda
                               parallel = getOption("ibis.runparallel"),
                               trace.it = settings$get("verbose"),
                               nfolds = 10  # number of folds for cross-validation
        )
      },silent = FALSE)
    }
    if(inherits(cv_gn, "try-error")) stop("Model failed to converge with provided input data!")

    # --- #
    # Predict spatially
    if(!settings$get('inference_only')){
      # Messenger
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
      # Make a subset of non-na values
      full$rowid <- 1:nrow(full)
      full_sub <- subset(full, stats::complete.cases(full))
      w_full_sub <- w_full[full_sub$rowid]
      assertthat::assert_that((nrow(full_sub) == length(w_full_sub)) || is.null(w_full_sub) )

      # Attempt prediction
      if(inherits(cv_gn, "cv.glmnet")){
        out <- predict(object = cv_gn,
                       newdata = full_sub,
                       weights = w_full_sub,
                       newoffset = ofs_pred[full_sub$rowid],
                       s = determine_lambda(cv_gn), # Determine the best lambda value
                       type = params$type
        )
      } else {
        # Assume cva.glmnet
        out <- predict(
          object = cv_gn,
          newdata = full_sub,
          alpha = cv_gn$alpha,
          weights = w_full_sub,
          newoffset = ofs_pred[full_sub$rowid],
          s = determine_lambda(cv_gn), # Determine the best lambda value
          type = params$type
        )
        # Determine best model based on cross-validated loss
        # ind <- which.min( sapply(cv_gn$modlist, function(z) min(z$cvup)) )
        # cv_gn <- cv_gn$modlist[[ind]]
      }


      # Fill output with summaries of the posterior
      prediction[full_sub$rowid] <- out[,1]
      names(prediction) <- "mean"
      prediction <- terra::mask(prediction, self$get_data("template"))
      try({rm(out, full, full_sub)},silent = TRUE)
    } else {
      # No prediction done
      prediction <- NULL
    }
    # Compute end of computation time
    settings$set('end.time', Sys.time())

    # Definition of GLMNET Model object ----
    obj <- DistributionModel # Make a copy to set new functions

    # Partial effects
    obj$set("public", "partial", function(x.var = NULL, constant = NULL, variable_length = 100,
                           values = NULL, newdata = NULL, plot = FALSE, type = NULL, ...){
      assertthat::assert_that(is.character(x.var) || is.null(x.var),
                              is.null(constant) || is.numeric(constant),
                              is.null(type) || is.character(type),
                              is.null(newdata) || is.data.frame(newdata),
                              is.numeric(variable_length)
      )
      # Settings
      settings <- self$settings

      mod <- self$get_data('fit_best')
      model <- self$model
      co <- stats::coef(mod) |> row.names() # Get model coefficient names
      # Set type
      if(is.null(type)) type <- self$settings$get("type")
      type <- match.arg(type, c("link", "response"), several.ok = FALSE)
      settings$set("type", type)

      # Get data
      df <- model$biodiversity[[length( model$biodiversity )]]$predictors
      df <- subset(df, select = all.vars(mod$terms))
      variables <- names(df)

      # Match x.var to argument
      if(is.null(x.var)){
        x.var <- variables
      } else {
        x.var <- match.arg(x.var, variables, several.ok = TRUE)
      }

      # Calculate range of predictors
      rr <- sapply(df[, names(df) %in% model$predictors_types$predictors[model$predictors_types$type=="numeric"]],
                   function(x) range(x, na.rm = TRUE)) |> as.data.frame()

      if(is.null(newdata)){
        # if values are set, make sure that they cover the data.frame
        if(!is.null(values)){
          assertthat::assert_that(length(x.var) == 1)
          df2 <- list()
          df2[[x.var]] <- values
          # Then add the others
          for(var in variables){
            if(var == x.var) next()
            df2[[var]] <- mean(df[[var]], na.rm = TRUE)
          }
          df2 <- df2 |> as.data.frame()
        } else {
          df2 <- list()
          for(i in x.var) {
            df2[[i]] <- as.data.frame(seq(rr[1,i],rr[2,i], length.out = variable_length))
          }
          df2 <- do.call(cbind, df2); names(df2) <- x.var
        }
      } else {
        # Assume all variables are present
        df2 <- dplyr::select(newdata, dplyr::any_of(names(df)))
        assertthat::assert_that(nrow(df2)>1, ncol(df2)>1)
      }

      # Get offset if set
      if(!is.Waiver(model$offset)){
        of <- model$offset$spatial_offset
      } else of <- new_waiver()

      # HACK: Overwrite lambda to make sure pdp uses it.
      mod$lambda.1se <- determine_lambda(mod)
      # Inverse link function
      ilf <- switch (settings$get('type'),
                     "link" = NULL,
                     "response" = ifelse(model$biodiversity[[1]]$family=='poisson',
                                         exp, logistic)
      )

      pp <- data.frame()
      pb <- progress::progress_bar$new(total = length(x.var))
      for(v in x.var){
        if(!is.Waiver(of)){
          # Predict with offset
          p1 <- pdp::partial(mod, pred.var = v, pred.grid = df2,
                             ice = FALSE, center = FALSE,
                             type = "regression", newoffset = of,
                             inv.link = ilf,
                             plot = FALSE, rug = TRUE, train = df)
        } else {
          p1 <- pdp::partial(mod, pred.var = v, pred.grid = df2,
                             ice = FALSE, center = FALSE,
                             type = "regression", inv.link = ilf,
                             plot = FALSE, rug = TRUE, train = df
          )
        }
        p1 <- p1[, c(v, "yhat")]
        names(p1) <- c("partial_effect", "mean")
        p1 <- cbind(variable = v, p1)
        pp <- rbind(pp, p1)
        rm(p1)
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
      return(pp)
    },overwrite = TRUE)

    # Spatial partial dependence plot
    obj$set("public", "spartial", function(x.var, constant = NULL, newdata = NULL,
                                           plot = TRUE, type = NULL){
      assertthat::assert_that(is.character(x.var),
                              "model" %in% names(self),
                              is.null(constant) || is.numeric(constant),
                              is.null(newdata) || is.data.frame(newdata),
                              is.logical(plot),
                              is.character(type) || is.null(type)
      )
      # Settings
      settings <- self$settings
      # Set type
      if(is.null(type)) type <- self$settings$get("type")
      type <- match.arg(type, c("link", "response"), several.ok = FALSE)
      settings$set("type", type)

      mod <- self$get_data('fit_best')
      model <- self$model
      # For Integrated model, take the last one
      fam <- model$biodiversity[[length(model$biodiversity)]]$family

      # If new data is set
      if(!is.null(newdata)){
        df <- newdata
      } else {
        df <- model$predictors
        df$w <- model$exposure
      }
      assertthat::assert_that(all(x.var %in% colnames(df)))
      df$rowid <- 1:nrow(df)
      # Match x.var to argument
      x.var <- match.arg(x.var, names(df), several.ok = FALSE)
      df_sub <- subset(df, stats::complete.cases(df))
      if(!is.Waiver(model$offset)) ofs <- model$offset[df_sub$rowid,3] else ofs <- NULL
      assertthat::assert_that(nrow(df_sub)>0)

      # Add all others as constant
      if(is.null(constant)){
        for(n in names(df_sub)) if(!n %in% c(x.var, "rowid", "w")) df_sub[[n]] <- suppressWarnings( mean(model$predictors[[n]], na.rm = TRUE) )
      } else {
        for(n in names(df_sub)) if(!n %in% c(x.var, "rowid", "w")) df_sub[[n]] <- constant
      }
      # Reclassify factor levels
      if(any(model$predictors_types$type=="factor")){
        lvl <- levels(model$predictors[[model$predictors_types$predictors[model$predictors_types$type=="factor"]]])
        df_sub[[model$predictors_types$predictors[model$predictors_types$type=="factor"]]] <-
          factor(lvl[1], levels = lvl)
      }
      # Predict with lambda
      pred_gn <- predict(
        object = mod,
        newdata = df_sub,
        weights = df_sub$w, # The second entry of unique contains the non-observed variables
        newoffset = ofs,
        s = determine_lambda(mod), # Determine best available lambda
        fam = fam,
        type = type
      ) |> as.data.frame()

      # Now create spatial prediction
      template <- model_to_background(model)
      template[df_sub$rowid] <- pred_gn[,1]
      names(template) <- paste0("spartial_",x.var)

      # Do plot and return result
      if(plot) terra::plot(template, col = ibis_colours$ohsu_palette)
      return(template)
    },overwrite = TRUE)

    # Convergence check
    obj$set("public", "has_converged", function(){
      fit <- self$get_data("fit_best")
      if(is.Waiver(fit)) return(FALSE)
      # Get lambdas
      lmd <- fit$lambda
      if(determine_lambda(fit) == min(lmd)) return(FALSE)
      return(TRUE)
    }, overwrite = TRUE)

    # Residual function
    obj$set("public", "get_residuals", function(){
      # Get best object
      obj <- self$get_data("fit_best")
      if(is.Waiver(obj)) return(obj)
      # Calculate residuals
      model <- self$model$predictors
      # Get fm
      if(inherits(obj, "cv.glmnet")){
        fitted_values <- predict(obj, model, s = 'lambda.1se')
        fitted_min <- predict(obj, model, s = 'lambda.min')
      } else {
        alpha <- sapply(obj$modlist, function(z) min(z$cvup))
        fitted_values <- predict(obj, model,
                                 which = which.min(alpha),
                                 s = 'lambda.1se')
        fitted_min <- predict(obj, model,
                              which = which.min(alpha),
                              s = 'lambda.min')
      }
      rd <- fitted_min[,1] - fitted_values[,1]
      assertthat::assert_that(length(rd)>0)
      return(rd)
    }, overwrite = TRUE)

    # Get coefficients from glmnet
    obj$set("public", "get_coefficients", function(){
      # Returns a vector of the coefficients with direction/importance
      obj <- self$get_data("fit_best")
      cofs <- tidy_glmnet_summary(obj)
      names(cofs) <- c("Feature", "Beta")
      # Remove intercept(s)
      int <- grep("Intercept",cofs$Feature,ignore.case = TRUE)
      if(length(int)>0) cofs <- cofs[-int,]
      return(cofs)
    },overwrite = TRUE)

    # Engine-specific projection function
    obj$set("public", "project", function(newdata, type = NULL, layer = "mean"){
      assertthat::assert_that("model" %in% names(self),
                              nrow(newdata) > 0,
                              all( c("x", "y") %in% names(newdata) ),
                              is.character(type) || is.null(type)
      )

      # Settings
      settings <- self$settings
      # Set type
      if(is.null(type)) type <- self$settings$get("type")
      type <- match.arg(type, c("link", "response"), several.ok = FALSE)
      settings$set("type", type)

      mod <- self$get_data('fit_best')
      model <- self$model
      # For Integrated model, take the last one
      fam <- model$biodiversity[[length(model$biodiversity)]]$family

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

      df <- newdata
      if(nrow(df)!=length(model$exposure)){
        # Assume the maximum (usually 1e6)
        df$w <- model$exposure[which.max(model$exposure)]
      } else {
        df$w <- model$exposure # Also get exposure variable
      }
      # Make a subset of non-na values
      df$rowid <- 1:nrow(df)
      df_sub <- base::subset(df, stats::complete.cases(df))
      if(!is.Waiver(model$offset)) ofs <- model$offset[df_sub$rowid] else ofs <- NULL
      assertthat::assert_that(nrow(df_sub)>0)

      if(inherits(mod, "cv.glmnet")){
        pred_gn <- predict(
          object = mod,
          newdata = df_sub,
          weights = df_sub$w, # The second entry of unique contains the non-observed variables
          newoffset = ofs,
          na.action = "na.pass",
          s = determine_lambda(mod), # Determine best available lambda
          fam = fam,
          type = type
        ) |> as.data.frame()
      } else {
        pred_gn <- predict(
          object = mod,
          newdata = df_sub,
          alpha = mod$alpha,
          weights = df_sub$w, # The second entry of unique contains the non-observed variables
          newoffset = ofs,
          na.action = "na.pass",
          s = determine_lambda(mod), # Determine the best lambda value
          fam = fam,
          type = type
        )
      }
      names(pred_gn) <- layer
      assertthat::assert_that(nrow(pred_gn)>0, nrow(pred_gn) == nrow(df_sub))

      # Now create spatial prediction
      if(nrow(newdata)==nrow(model$predictors)){
        prediction <- try({model_to_background(model)}, silent = TRUE)
        prediction[df_sub$rowid] <- pred_gn[, layer]
      } else {
        assertthat::assert_that(utils::hasName(df_sub,"x")&&utils::hasName(df_sub,"y"),
                                msg = "Projection data.frame has no valid coordinates or differs in grain!")
        prediction <- try({
          terra::rast(df_sub[,c("x", "y")],
                      crs = terra::crs(model$background),
                      type = "xyz") |>
            emptyraster()
        }, silent = TRUE)
        prediction[] <- pred_gn[, layer]
      }

      return(prediction)
    },overwrite = TRUE)

    out <- obj$new(name = "GLMNET-Model")
    out$id <- model$id
    out$model <- model
    out$settings <- settings
    out$fits = list(
      "fit_best" = cv_gn,
      "fit_best_equation" = form,
      "prediction" = prediction
    )

    return(out)
  },overwrite = TRUE)

  # Define engine object
  eg <- eg$new(engine = "GLMNET-Engine", name = "<GLMNET>")
  eg$data <- list(
    'template' = template,
    'params' = params
  )
  # Set engine in distribution object
  y <- x$clone(deep = TRUE)
  return( y$set_engine(eg) )
} # End of function
