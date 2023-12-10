#' @include bdproto-engine.R utils-spatial.R bdproto-distributionmodel.R
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
#' @details
#' This engine is essentially a wrapper for [stats::glm.fit()], however with customized
#' settings to support offsets and weights.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param control A [`list`] containing parameters for controlling the fitting
#' process (Default: \code{NULL}).
#' @param type The mode used for creating posterior predictions. Either making
#'   \code{"link"} or \code{"response"} (Default: \code{"response"}).
#' @param ... Other parameters passed on to [stats::glm()].
#' @references
#' * Hastie, T. J. and Pregibon, D. (1992) Generalized linear models. Chapter 6 of Statistical Models in S eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
#' @family engine
#' @returns An [Engine].
#' @aliases engine_glm
#' @examples
#' # Load background
#' background <- terra::rast(system.file('extdata/europegrid_50km.tif', package='ibis.iSDM',mustWork = TRUE))
#'
#' # Add GLM as an engine
#' x <- distribution(background) |> engine_glm()
#' x
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

  # Set engine in distribution object
  x$set_engine(
    bdproto(
      "GLM-Engine",
      Engine,
      name = "<GLM>",
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
          nrow(model$predictors) == terra::ncell(self$get_data('template')),
          !is.Waiver(self$get_data("params")),
          length(model$biodiversity) == 1 # Only works with single likelihood. To be processed separately
        )
        # Messenger
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Engine setup.')

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
      },
      # Training function
      train = function(self, model, settings, ...){
        assertthat::assert_that(
          inherits(settings,'Settings'),
          is.list(model),length(model)>1,
          # Check that model id and setting id are identical
          settings$modelid == model$id
        )
        # Get name
        name <- model$biodiversity[[1]]$name

        # Messenger
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','green',paste0('Starting fitting: ', name))

        # Verbosity
        verbose <- settings$get("verbose")

        # Set prediction type also for later
        settings$set('type', self$get_data("params")$type)

        # seed
        seed <- settings$get("seed")
        if(is.Waiver(seed)) { settings$set('seed', getOption("ibis.seed")) }

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
            if(getOption('ibis.setupmessages')) myLog('[Estimation]','red',paste0("Custom link functions not supported!"))
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

        # Clamp?
        if( settings$get("clamp") ) full <- clamp_predictions(model, full)

        # -- #
        # Expand predictors if non-linear is specified in settings
        if(settings$get('only_linear') == FALSE){
          if(getOption('ibis.setupmessages')) myLog('[Estimation]','yellow',
                                                    'Non-linearity to glm is best introduced by adding derivates. Ignored!')
        }

        assertthat::assert_that(
          is.null(w) || length(w) == nrow(df),
          is.null(ofs) || is.vector(ofs),
          is.null(ofs_pred) || is.vector(ofs_pred),
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
          if(getOption('ibis.setupmessages')) myLog('[Estimation]','yellow',
                                                    'No hyperparameter optimization for glm implemented!')
        } else {
          if(!is.null(ofs)){
            # Split since GLM cannot handle NULL offsets
            suppressWarnings(
              fit_glm <- try({
                stats::glm(formula = form,
                           data = df,
                           weights = w, # Case weights
                           offset = ofs,
                           family = fam,
                           na.action = "na.pass",
                           control = params$control
                )
              },silent = FALSE)
            )
          } else {
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
          }
        }
        if(inherits(fit_glm, "try-error")) stop("Model failed to converge with provided input data!")

        # --- #
        # Predict spatially
        if(!settings$get('inference_only')){
          # Messenger
          if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Starting prediction...')
          # Set target variables to bias_value for prediction if specified
          if(!is.Waiver(settings$get('bias_variable'))){
            for(i in 1:length(settings$get('bias_variable'))){
              if(settings$get('bias_variable')[i] %notin% names(full)){
                if(getOption('ibis.setupmessages')) myLog('[Estimation]','red','Did not find bias variable in prediction object!')
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
          out <- try({
            stats::predict.glm(object = fit_glm,
                               newdata = full,
                               type = params$type,
                               se.fit = TRUE,
                               na.action = "na.pass",
                               weights = w_full
                               )
          },silent = TRUE)
          if(!inherits(out,"try-error")){
            # Fill output with summaries of the posterior
            prediction <- fill_rasters(out |> as.data.frame(),
                                       background = prediction)[[1:2]]
            names(prediction) <- c("mean", "se")
            prediction <- terra::mask(prediction, self$get_data("template"))

          } else {
            stop("GLM prediction failed!")
          }
          try({rm(out, full, full_sub)},silent = TRUE)
        } else {
          # No prediction done
          prediction <- NULL
        }
        # Compute end of computation time
        settings$set('end.time', Sys.time())

        # Definition of GLMNET Model object ----
        # Create output
        out <- bdproto(
          "GLM-Model",
          DistributionModel,
          id = model$id,
          model = model,
          settings = settings,
          fits = list(
            "fit_best" = fit_glm,
            "fit_best_equation" = form,
            "prediction" = prediction
          ),
          # Partial effects
          partial = function(self, x.var = NULL, constant = NULL, variable_length = 100,
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
            df <- model$biodiversity[[length( model$biodiversity )]]$predictors

            # Match x.var to argument
            if(is.null(x.var)){
              x.var <- colnames(df)
            } else {
              x.var <- match.arg(x.var, names(df), several.ok = FALSE)
            }

            # Calculate range of predictors
            if(any(model$predictors_types$type=="factor")){
              rr <- sapply(df[model$predictors_types$predictors[model$predictors_types$type=="numeric"]],
                           function(x) range(x, na.rm = TRUE)) |> as.data.frame()
            } else {
              rr <- sapply(df, function(x) range(x, na.rm = TRUE)) |> as.data.frame()
            }

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
              df2 <- newdata |> dplyr::select(dplyr::any_of(names(df)))
              assertthat::assert_that(nrow(df2)>1, ncol(df2)>1)
            }

            # Get offset if set
            if(!is.Waiver(model$offset)){
              of <- model$offset$spatial_offset
            } else of <- new_waiver()

            # Check that variables are in
            assertthat::assert_that(all( x.var %in% colnames(df) ),
                                    msg = 'Variable not in predicted model.')

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
              p1 <- p1[,c(v, "yhat")]
              names(p1) <- c("partial_effect", "mean")
              p1$variable <- v
              pp <- rbind(pp, p1)
              rm(p1)
              if(length(x.var) > 1) pb$tick()
            }

            if(plot){
              # Make a plot
              g <- ggplot2::ggplot(data = pp, ggplot2::aes(x = partial_effect, y = mean)) +
                ggplot2::theme_classic(base_size = 18) +
                ggplot2::geom_line() +
                ggplot2::labs(x = "", y = expression(hat(y))) +
                ggplot2::facet_wrap(~variable,scales = 'free')
              print(g)
            }
            return(pp)
          },
          # Spatial partial dependence plot
          spartial = function(self, x.var, constant = NULL, newdata = NULL, plot = TRUE, type = NULL){
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
          # Convergence check
          has_converged = function(self){
            obj <- self$get_data("fit_best")
            return( obj$converged )
          },
          # Residual function
          get_residuals = function(self, type = NULL){
            # Get best object
            obj <- self$get_data("fit_best")
            if(is.Waiver(obj)) return(obj)
            settings <- self$settings
            if(is.null(type)) type <- settings$get('type')
            # Calculate residuals
            rd <- stats::residuals.glm(obj, type = type)
            return(rd)
          },
          # Get coefficients from glmnet
          get_coefficients = function(self){
            # Returns a vector of the coefficients with direction/importance
            obj <- self$get_data("fit_best")
            cofs <- tidy_glm_summary(obj)
            names(cofs)[1:2] <- c("Feature", "Beta")
            return(cofs)
          },
          # Engine-specific projection function
          project = function(self, newdata, type = NULL, layer = "mean"){
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
                  if(getOption('ibis.setupmessages')) myLog('[Estimation]','red','Did not find bias variable in prediction object!')
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

            names(pred_glm) <- layer
            assertthat::assert_that(nrow(pred_glm)>0, nrow(pred_glm) == nrow(df))

            # Now create spatial prediction
            prediction <- fill_rasters(pred_glm,
                                       model_to_background(model)
                                       )

            return(prediction)
          }
        )
        return(out)
      }
    )
  ) # End of bdproto object
} # End of function
