#' @include class-engine.R class-distributionmodel.R
NULL

#' Engine for Generalized linear models (GLM)
#'
#' @description
#' This engine implements a basic generalized linear modle (GLM) for creating
#' species distribution models. The main purpose of this engine is to support
#' a basic, dependency-free method for inference and projection that can be used
#' within the package for examples and vignettes. That being said, the engine is
#' fully functional as any other engine.
#'
#' The basic implementation of GLMs here is part of a general class oflinear models
#' and has - with exception of offsets - only minimal options to integrate other
#' sources of information such as priors or joint integration. The general
#' recommendation is to [engine_glmnet()] instead for regularization support.
#' However basic GLMs can in some cases be useful for quick projections or
#' for [ensemble()] of small models (a practice common for rare species).
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param control A [`list`] containing parameters for controlling the fitting
#' process (Default: \code{NULL}).
#' @param type The mode used for creating posterior predictions. Either making
#' \code{"link"} or \code{"response"} (Default: \code{"response"}).
#' @param ... Other parameters passed on to [stats::glm()].
#'
#' @details
#' This engine is essentially a wrapper for [stats::glm.fit()], however with customized
#' settings to support offsets and weights.
#'
#' If \code{"optim_hyperparam"} is set to \code{TRUE} in [`train()`], then a AIC
#' based step-wise (backwards) model selection is performed.
#' Generally however [`engine_glmnet`] should be the preferred package for models
#' with more than \code{>3} covariates.
#'
#' @returns An [Engine].
#'
#' @references
#' * Hastie, T. J. and Pregibon, D. (1992) Generalized linear models. Chapter 6 of
#' Statistical Models in S eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
#'
#' @family engine
#'
#' @examples
#' # Load background
#' background <- terra::rast(system.file('extdata/europegrid_50km.tif',
#' package='ibis.iSDM',mustWork = TRUE))
#'
#' # Add GLM as an engine
#' x <- distribution(background) |> engine_glm()
#' print(x)
#'
#' @name engine_glm
NULL

#' @rdname engine_glm
#' @export
engine_glm <- function(x,
                       control = NULL,
                       type = "response",
                       ...) {

  # Check whether package is available (Default installed)
  check_package('stats')

  # assert that arguments are valid
  assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                          inherits(x$background,'sf'),
                          is.null(control) || is.list(control),
                          is.character(type)
  )
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

  # Specify default control
  if(is.null(control)){
    control <- stats::glm.control()
  }

  # Set up the parameter list
  params <- list(
    control = control,
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
  }, overwrite = TRUE)

  # Dummy function for getting the equation of latent effects
  eg$set("public", "get_equation_latent_spatial", function(method){
    new_waiver()
  }, overwrite = TRUE)

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
    assertthat::assert_that(nrow(df)<=nrow(full)) # Security check?

    # Subset the predictor types to only those present
    te <- formula_terms(form)
    model$biodiversity[[1]]$predictors_types <- dplyr::filter(model$biodiversity[[1]]$predictors_types, predictors %in% te)
    model$biodiversity[[1]]$predictors_names <-  intersect(model$biodiversity[[1]]$predictors_names, te)

    # Get offset and add it to exposure
    if(!is.Waiver(model$offset)){
      # Add offset to full prediction and load vector
      ofs <- model$biodiversity[[1]]$offset[, 'spatial_offset']
      ofs_pred <- model$offset[,'spatial_offset']
      # Add to data.frame and form
      form <- stats::update.formula(form, . ~ . + offset(spatial_offset))
      df$spatial_offset <- ofs
      full$spatial_offset <- ofs_pred
    } else { ofs <- NULL; ofs_pred <- NULL }

    # Clamp?
    if( settings$get("clamp") ) full <- clamp_predictions(model, full)

    # -- #
    # Expand predictors if non-linear is specified in settings
    if(settings$get('only_linear') == FALSE){
      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','yellow',
                                                'Non-linearity to glm is best introduced by adding derivates. Ignored!')
    }

    assertthat::assert_that(
      is.null(w) || length(w) == nrow(df),
      is.null(ofs) || is.vector(ofs),
      is.null(ofs_pred) || is.vector(ofs_pred),
      all(w >= 0,na.rm = TRUE)
    )
    # --- #
    # Fit Base model
    suppressWarnings(
      fit_glm <- try({
        stats::glm(formula = form,
                   data = df,
                   weights = w, # Case weights
                   family = fam,
                   na.action = "na.pass",
                   control = params$control
        )
      },silent = FALSE)
    )
    if(inherits(fit_glm, "try-error")) stop("Model failed to converge with provided input data!")
    if( (settings$get('optim_hyperparam')) ){
      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','green',
                                                                'Running step-wise AIC selection for glm!')
      suppressWarnings(
        fit_glm <- stats::step(fit_glm,
                             direction = "backward",
                             trace = ifelse(getOption('ibis.setupmessages', default = TRUE),1,0)
                             )
      )
    }

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

      # Attempt prediction
      if( getOption('ibis.runparallel',default = FALSE) ){
        check_package("doFuture")
        if(!("doFuture" %in% loadedNamespaces()) || ('doFuture' %notin% utils::sessionInfo()$otherPkgs) ) {
          try({requireNamespace('doFuture');attachNamespace("doFuture")},silent = TRUE)
        }

        # Prediction function
        do_run <- function() {
          # Chunck the data
          splits <- chunk_data(full, N = getOption("ibis.nthread"),index_only = TRUE)

          y <- foreach::foreach(s = splits,
                                .inorder = TRUE
                                # .options.future = list(globals = ))
                                ) %dofuture% {
            stats::predict.glm(object = fit_glm,
                               newdata = full[s,],
                               type = params$type,
                               se.fit = TRUE,
                               na.action = "na.pass",
                               weights = w_full[s,]
            )
          }
          y
        }
        # Run
        result <- do_run()
        # Combine all
        # FIXME: hacky list flattener, but works. Reduce and do.call failed
        out <- list()
        for(k in 1:length(result)){
          out[['fit']] <- c(out[['fit']], result[[k]]$fit)
          out[['se.fit']] <- c(out[['se.fit']], result[[k]]$se.fit)
        }
        # Security check
        assertthat::assert_that(
          length(out$fit) == nrow(full)
        )
      } else {
        out <- try({
          stats::predict.glm(object = fit_glm,
                             newdata = full,
                             type = params$type,
                             se.fit = TRUE,
                             na.action = "na.pass",
                             weights = w_full
          )
        },silent = TRUE)
      }
      if(!inherits(out,"try-error")){
        # Fill output with summaries of the posterior
        prediction <- fill_rasters(out |> as.data.frame(),
                                   background = prediction)[[1:2]]
        names(prediction) <- c("mean", "se")
        prediction <- terra::mask(prediction, self$get_data("template"))

      } else {
        stop("GLM prediction failed!")
      }
      try({rm(out, full)},silent = TRUE)
    } else {
      # No prediction done
      prediction <- NULL
    }
    # Compute end of computation time
    settings$set('end.time', Sys.time())

    # Definition of GLMNET Model object ----
    obj <- DistributionModel # Make a copy to set new functions

    # Partial effect functions
    obj$set("public", "partial",
            function(x.var = NULL, constant = NULL, variable_length = 100,
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
              co <- stats::coefficients(mod) |> names() # Get model coefficient names
              # Set type
              if(is.null(type)) type <- self$settings$get("type")
              type <- match.arg(type, c("link", "response"), several.ok = FALSE)
              settings$set("type", type)

              # Get data
              df <- model$biodiversity[[length(model$biodiversity)]]$predictors
              df <- subset(df, select = attr(mod$terms, "term.labels"))
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
                  for(var in colnames(df)){
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
                df2 <- dplyr::select(newdata, dplyr::any_of(variables))
                assertthat::assert_that(nrow(df2)>1, ncol(df2)>1)
              }

              # Get offset if set
              if(!is.Waiver(model$offset)){
                of <- model$offset$spatial_offset
              } else of <- new_waiver()

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
            }, overwrite = TRUE)

    # Spatial partial dependence plot
    obj$set("public", "spartial", function(x.var, constant = NULL, newdata = NULL, plot = TRUE, type = NULL){
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

      # Add all others as constant
      if(is.null(constant)){
        for(n in names(df)) if(!n %in% c(x.var, "rowid", "w")) df[[n]] <- suppressWarnings( mean(model$predictors[[n]], na.rm = TRUE) )
      } else {
        for(n in names(df)) if(!n %in% c(x.var, "rowid", "w")) df[[n]] <- constant
      }
      # Reclassify factor levels
      if(any(model$predictors_types$type=="factor")){
        lvl <- levels(model$predictors[[model$predictors_types$predictors[model$predictors_types$type=="factor"]]])
        df[[model$predictors_types$predictors[model$predictors_types$type=="factor"]]] <-
          factor(lvl[1], levels = lvl)
      }

      # Predict
      pred_glm <- stats::predict.glm(
        object = mod,
        newdata = df,
        weights = df$w, # The second entry of unique contains the non-observed variables
        se.fit = FALSE,
        na.action = "na.pass",
        fam = fam,
        type = type
      ) |> as.data.frame()
      assertthat::assert_that(nrow(pred_glm)>0, nrow(pred_glm) == nrow(df))

      # Now create spatial prediction
      prediction <- fill_rasters(pred_glm, model_to_background(model))
      names(prediction) <- paste0("spartial_",x.var)

      # Do plot and return result
      if(plot) terra::plot(prediction, col = ibis_colours$ohsu_palette)
      return(prediction)
    },
    overwrite = TRUE)

    # Convergence check
    obj$set("public", "has_converged", function(){
      obj <- self$get_data("fit_best")
      return( obj$converged )
    }, overwrite = TRUE)

    # Residual function
    obj$set("public", "get_residuals", function(type = NULL){
      # Get best object
      obj <- self$get_data("fit_best")
      if(is.Waiver(obj)) return(obj)
      settings <- self$settings
      if(is.null(type)) type <- settings$get('type')
      # Calculate residuals
      rd <- stats::residuals.glm(obj, type = type)
      return(rd)
    }, overwrite = TRUE)

    # Get coefficients from glmnet
    obj$set("public", "get_coefficients", function(){
      # Returns a vector of the coefficients with direction/importance
      obj <- self$get_data("fit_best")
      cofs <- tidy_glm_summary(obj)
      names(cofs)[1:2] <- c("Feature", "Beta")
      return(cofs)
    }, overwrite = TRUE)

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
      df$w <- model$exposure # Also get exposure variable
      df$rowid <- 1:nrow(df)
      if(!is.Waiver(model$offset)) ofs <- model$offset else ofs <- NULL
      assertthat::assert_that(nrow(df)>0)

      # Run in parallel if specified or not.
      if( getOption('ibis.runparallel',default = FALSE) ){
        check_package("doFuture")
        if(!("doFuture" %in% loadedNamespaces()) || ('doFuture' %notin% utils::sessionInfo()$otherPkgs) ) {
          try({requireNamespace('doFuture');attachNamespace("doFuture")},silent = TRUE)
        }

        # Prediction function
        do_run <- function() {
          # Chunck the data
          splits <- chunk_data(df, N = getOption("ibis.nthread",default = 10),index_only = TRUE)

          y <- foreach::foreach(s = splits,
                                .inorder = TRUE
                                # .options.future = list(globals = ))
          ) %dofuture% {
            stats::predict.glm(object = mod,
                               newdata = df[s,],
                               type = type,
                               se.fit = FALSE,
                               na.action = "na.pass",
                               weights = df$w[s]
            )
          }
          y
        }
        # Run
        result <- do_run()
        # Combine all
        # FIXME: hacky list flattener, but works. Reduce and do.call failed
        pred_glm <- list()
        for(k in 1:length(result)){
          pred_glm[['fit']] <- c(pred_glm[['fit']], result[[k]]$fit)
        }
        pred_glm <- as.data.frame(pred_glm)

      } else {
        if(is.null(ofs)){
          pred_glm <- stats::predict.glm(
            object = mod,
            newdata = df,
            weights = df$w, # The second entry of unique contains the non-observed variables
            se.fit = FALSE,
            na.action = "na.pass",
            fam = fam,
            type = type
          ) |> as.data.frame()
        } else {
          pred_glm <- stats::predict.glm(
            object = mod,
            newdata = df,
            weights = df$w, # The second entry of unique contains the non-observed variables
            offset = ofs,
            se.fit = FALSE,
            na.action = "na.pass",
            fam = fam,
            type = type
          ) |> as.data.frame()
        }
      }
      names(pred_glm) <- layer
      assertthat::assert_that(nrow(pred_glm)>0, nrow(pred_glm) == nrow(df))

      # Now create spatial prediction
      if(nrow(newdata)==nrow(model$predictors)){
        prediction <- try({model_to_background(model)}, silent = TRUE)
      } else {
        assertthat::assert_that(utils::hasName(df,"x")&&utils::hasName(df,"y"),
                                msg = "Projection data.frame has no valid coordinates or differs in grain!")
        prediction <- try({
          terra::rast(df[,c("x", "y")],
                      crs = terra::crs(model$background),
                      type = "xyz") |>
            emptyraster()
        }, silent = TRUE)
      }
      prediction <- fill_rasters(pred_glm, prediction)
      return(prediction)
    }, overwrite = TRUE)

    # Now generate output
    out <- obj$new(name = "GLM-Model")
    out$id <- model$id
    out$model <- model
    out$settings <- settings
    out$fits <- list(
      "fit_best" = fit_glm,
      "fit_best_equation" = form,
      "prediction" = prediction
    )

    return(out)
  },overwrite = TRUE)

  # Define engine object and save
  eg <- eg$new(engine = "GLM-Engine", name = "<GLM>")
  eg$data <- list(
    'template' = template,
    'params' = params
  )

  # Set engine in distribution object
  y <- x$clone(deep = TRUE)
  return( y$set_engine(eg) )
} # End of function
