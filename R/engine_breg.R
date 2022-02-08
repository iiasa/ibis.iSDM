#' @include bdproto-engine.R utils-spatial.R bdproto-distributionmodel.R
NULL

#' Engine for Bayesian regularized regression models
#'
#' @description
#' Efficient MCMC algorithm for linear regression models that makes use of
#' 'spike-and-slab' priors for some modest regularization on the amount of posterior
#' probability for a subset of the coefficients.
#' @details
#' This engine provides efficient Bayesian predictions through the \pkg{Boom} R-package. However note
#' that not all link and models functions are supported and certain functionalities such as offsets are generally
#' not available.
#' This engines allows the estimation of linear and non-linear effects via the \code{"only_linear"} parameter
#' specified in [train]. If set to \code{FALSE}, non-linearity will be approximated through a combination of
#' B- and I-Splines.
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
        settings$set('iter', params$iter)
        settings$set('type', params$type)

        # Distribution specific procedure
        fam <- model$biodiversity[[1]]$family

        # If offset is specified, raise warning
        # if(!is.Waiver(model$offset)){
        #   if(getOption('ibis.setupmessages')) myLog('[Estimation]','red','Engine breg can not use offsets. Ignored for now.')
        # }

        # Check whether regularization parameter is set to none, if yes, raise message
        if(settings$get("varsel") == "none"){
          if(getOption('ibis.setupmessages')) myLog('[Estimation]','yellow','Note: Engine_breg always applies regularization.')
        }

        # -- #
        # Expand predictors if specified in settings
        if(settings$get('only_linear') == FALSE){
          if(getOption('ibis.setupmessages')) myLog('[Estimation]','yellow','Non-linear estimation not added to engine. Suggest to create variable derivatives externally.')
        }
        # -- #

        # If a poisson family is used, weight the observations by their exposure
        if(fam == "poisson"){
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
          # Add offset if existent
          if(!is.Waiver(model$offset)){
            ofs <- get_ngbvalue(coords = df[,c('x','y')],
                                env =  model$offset,
                                longlat = raster::isLonLat(self$get_data('template')),
                                field_space = c('x','y')
            )
            model$biodiversity[[1]]$offset <- ofs
          }

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

        # Set prediction type also for later
        settings$set('type', self$get_data("params")$type)

        # seed
        seed <- settings$get("seed")
        if(is.Waiver(seed)) { seed <- 1337; settings$set('seed', 1337) }

        # Get output raster
        prediction <- self$get_data('template')

        # Get parameters control
        params <- self$get_data('params')

        # All other needed data for model fitting
        fam <- model$biodiversity[[1]]$family
        form <- model$biodiversity[[1]]$equation
        df <- cbind(model$biodiversity[[1]]$predictors,
                    data.frame(observed = model$biodiversity[[1]]$observations[,'observed'])
                    )
        df <- subset(df, select = c(model$biodiversity[[1]]$predictors_names, "observed"))
        w <- model$biodiversity[[1]]$expect # The expected exposure
        # Get full prediction container
        full <- model$predictors
        w_full <- model$exposure

        # Priors
        if(!is.Waiver(model$priors)){
          # Define a family specific Boom prior
          pp <- setup_prior_boom(form = form,
                                 data = df,
                                 prior = model$priors,
                                 family = fam,
                                 exposure = w
                                 )
        } else { pp <- NULL }

        # Get offset and add it to exposure
        if(!is.Waiver(model$offset)){
          # Add offset to full prediction and load vector
          w <- w + model$biodiversity[[1]]$offset[, 'spatial_offset']
          w_full <- w_full + model$offset[,'spatial_offset']
          # negative exposure does not work, so normalize again to range of 1e-6 to 1
          if(any(w < 0,na.rm = TRUE)) {
            w <- scales::rescale(w, to = c(1e-6, 1))
            w_full <- scales::rescale(w_full, to = c(1e-6, 1))
          }
          if(anyNA(w)){
            w[is.na(w)] <- 1e-6
            w_full[is.na(w_full)] <- 1e-6
          }
        }

        assertthat::assert_that(
          is.null(w) || length(w) == nrow(df),
          all(w >= 0,na.rm = TRUE) # Required for engine_breg
        )
        # --- #
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
          # -> external code in utils-boom
          pred_breg <- predict_boom(
            obj = fit_breg,
            newdata = full,
            w = w_full,
            fam = fam,
            params = params
          )

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
          partial = function(self, x.var = NULL, constant = NULL, length.out = 100, plot = FALSE, type = NULL, ...){
            assertthat::assert_that(is.character(x.var) || is.null(x.var),
                                    is.null(constant) || is.numeric(constant),
                                    is.null(type) || is.character(type),
                                    is.numeric(length.out)
            )
            # Settings
            settings <- self$settings
            # Set type
            if(is.null(type)) type <- self$settings$get("type")
            type <- match.arg(type, c("link", "response"), several.ok = FALSE)
            settings$set("type", type)

            mod <- self$get_data('fit_best')
            df <- self$model$biodiversity[[length( self$model$biodiversity )]]$predictors
            df <- subset(df, select = attr(mod$terms, "term.labels"))
            w <- self$model$biodiversity[[1]]$expect # Also get exposure variable

            # Match x.var to argument
            if(is.null(x.var)){
              x.var <- colnames(df)
            } else {
              x.var <- match.arg(x.var, names(df), several.ok = FALSE)
            }

            # Calculate range of predictors
            rr <- sapply(df, function(x) range(x, na.rm = TRUE)) |> as.data.frame()
            df_partial <- list()

            # Add all others as constant
            if(is.null(constant)){
              for(n in names(rr)) df_partial[[n]] <- rep( mean(df[[n]], na.rm = TRUE), length.out )
            } else {
              for(n in names(rr)) df_partial[[n]] <- rep( constant, length.out )
            }

            df_partial[[x.var]] <- seq(rr[1,x.var], rr[2,x.var], length.out = length.out)
            df_partial <- df_partial %>% as.data.frame()

            # For Integrated model, take the last one
            fam <- model$biodiversity[[length(model$biodiversity)]]$family

            pred_breg <- predict_boom(
              obj = mod,
              newdata = df_partial,
              w = unique(w)[2], # The second entry of unique contains the non-observed variables
              fam = fam,
              params = settings$data # Use the settings as list
            )            # Also attach the partial variable

            # Summarize the partial effect
            pred_part <- cbind(
              matrixStats::rowMeans2(pred_breg, na.rm = TRUE),
              matrixStats::rowSds(pred_breg, na.rm = TRUE),
              matrixStats::rowQuantiles(pred_breg, probs = c(.05,.5,.95), na.rm = TRUE),
              apply(pred_breg, 1, mode)
            ) %>% as.data.frame()
            names(pred_part) <- c("mean", "sd", "q05", "q50", "q95", "mode")
            pred_part$cv <- pred_part$mean / pred_part$sd
            # And attach the variable
            pred_part <- cbind("partial_effect" = df_partial[[x.var]], pred_part)

            if(plot){
              # Make a plot
              ggplot2::ggplot(data = pred_part, ggplot2::aes(x = partial_effect, y = q50, ymin = q05, ymax = q95)) +
                ggplot2::theme_classic(base_size = 18) +
                ggplot2::geom_ribbon(fill = 'grey90') +
                ggplot2::geom_line() +
                ggplot2::labs(x = paste0("partial of ",x.var), y = expression(hat(y)))
            }
            # Return the data
            return(pred_part)
          },
          # Spatial partial dependence plot
          spartial = function(self, x.var, constant = NULL, plot = TRUE, type = NULL){
            assertthat::assert_that(is.character(x.var) || is.null(x.var),
                                    "model" %in% names(self),
                                    is.null(constant) || is.numeric(constant),
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
            df <- model$biodiversity[[length( model$biodiversity )]]$predictors
            df <- subset(df, select = attr(mod$terms, "term.labels"))
            w <- model$biodiversity[[1]]$expect # Also get exposure variable

            # Match x.var to argument
            if(is.null(x.var)){
              x.var <- colnames(df)
            } else {
              x.var <- match.arg(x.var, names(df), several.ok = FALSE)
            }

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

            # For Integrated model, take the last one
            fam <- model$biodiversity[[length(model$biodiversity)]]$family

            pred_breg <- predict_boom(
              obj = mod,
              newdata = df_partial@data,
              w = unique(w)[2], # The second entry of unique contains the non-observed variables
              fam = fam,
              params = settings$data # Use the settings as list
            )

            # Summarize the partial effect
            pred_part <- cbind(
              matrixStats::rowMeans2(pred_breg, na.rm = TRUE),
              matrixStats::rowSds(pred_breg, na.rm = TRUE),
              matrixStats::rowQuantiles(pred_breg, probs = c(.05,.5,.95), na.rm = TRUE),
              apply(pred_breg, 1, mode)
            ) %>% as.data.frame()
            names(pred_part) <- c("mean", "sd", "q05", "q50", "q95", "mode")
            pred_part$cv <- pred_part$mean / pred_part$sd

            # Now create spatial prediction
            prediction <- emptyraster( self$get_data('prediction') ) # Background
            prediction <- fill_rasters(pred_part, prediction)

            # Do plot and return result
            if(plot){
              plot(prediction, col = ibis_colours$viridis_orig)
            }
            return(prediction)
          },
          # Engine-specific projection function
          project = function(self, newdata, type = NULL){
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
            df <- model$biodiversity[[length( model$biodiversity )]]$predictors
            df <- subset(df, select = attr(mod$terms, "term.labels"))
            w <- model$biodiversity[[1]]$expect # Also get exposure variable

            # Make spatial container for prediction
            suppressWarnings(
              df_partial <- sp::SpatialPointsDataFrame(coords = model$predictors[,c('x', 'y')],
                                                       data = model$predictors[, names(model$predictors) %notin% c('x','y')],
                                                       proj4string = sp::CRS( sp::proj4string(as(model$background, "Spatial")) )
              )
            )
            df_partial <- as(df_partial, 'SpatialPixelsDataFrame')
            # For Integrated model, take the last one
            fam <- model$biodiversity[[length(model$biodiversity)]]$family

            pred_breg <- predict_boom(
              obj = mod,
              newdata = df_partial@data,
              w = unique(w)[2], # The second entry of unique contains the non-observed variables
              fam = fam,
              params = settings$data # Use the settings as list
            )

            # Summarize the partial effect
            pred_part <- cbind(
              matrixStats::rowMeans2(pred_breg, na.rm = TRUE),
              matrixStats::rowSds(pred_breg, na.rm = TRUE),
              matrixStats::rowQuantiles(pred_breg, probs = c(.05,.5,.95), na.rm = TRUE),
              apply(pred_breg, 1, mode)
            ) %>% as.data.frame()
            names(pred_part) <- c("mean", "sd", "q05", "q50", "q95", "mode")
            pred_part$cv <- pred_part$mean / pred_part$sd

            # Now create spatial prediction
            prediction <- emptyraster( self$get_data('prediction') ) # Background
            prediction <- fill_rasters(pred_part, prediction)

            return(prediction)
          }
        )
        return(out)
      }
    )
  ) # End of bdproto object
} # End of function
