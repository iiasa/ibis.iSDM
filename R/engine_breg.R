#' @include bdproto-engine.R utils-spatial.R bdproto-distributionmodel.R
NULL

#' Engine for Bayesian regularized regression models
#'
#' @description Efficient MCMC algorithm for linear regression models that makes
#' use of 'spike-and-slab' priors for some modest regularization on the amount
#' of posterior probability for a subset of the coefficients.
#' @details This engine provides efficient Bayesian predictions through the
#' \pkg{Boom} R-package. However note that not all link and models functions are
#' supported and certain functionalities such as offsets are generally not
#' available. This engines allows the estimation of linear and non-linear
#' effects via the \code{"only_linear"} option specified in [train].
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param iter [`numeric`] on the number of MCMC iterations to run (Default:
#'   \code{10000}).
#' @param nthread [`numeric`] on the number of CPU-threads to use for data
#'   augmentation.
#' @param type The mode used for creating posterior predictions. Either making
#'   \code{"link"} or \code{"response"} (Default: \code{"response"}).
#' @param ... Other none specified parameters passed on to the model.
#' @references
#' * Nguyen, K., Le, T., Nguyen, V., Nguyen, T., & Phung, D. (2016, November). Multiple kernel learning with data augmentation. In Asian Conference on Machine Learning (pp. 49-64). PMLR.
#' * Steven L. Scott (2021). BoomSpikeSlab: MCMC for Spike and Slab Regression. R package version 1.2.4. https://CRAN.R-project.org/package=BoomSpikeSlab
#' @family engine
#' @returns An [Engine].
#' @aliases engine_breg
#' @examples
#' \dontrun{
#' # Add BREG as an engine
#' x <- distribution(background) |> engine_breg(iter = 1000)
#' }
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
  if(!("BoomSpikeSlab" %in% loadedNamespaces()) || ('BoomSpikeSlab' %notin% utils::sessionInfo()$otherPkgs) ) {
    try({requireNamespace('BoomSpikeSlab');attachNamespace("BoomSpikeSlab")},silent = TRUE)
  }

  # assert that arguments are valid
  assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                          inherits(x$background,'sf'),
                          is.numeric(iter),
                          is.character(type),
                          is.numeric(nthread)
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

        # Check whether regularization parameter is set to none, if yes, raise message
        if(settings$get('optim_hyperparam')){
          if(getOption('ibis.setupmessages')) myLog('[Estimation]','yellow','Note: Engine_breg always applies regularization.')
        }

        # -- #
        # Expand predictors if specified in settings
        if(settings$get('only_linear') == FALSE){
          if(getOption('ibis.setupmessages')) myLog('[Estimation]','yellow','Non-linear estimation not added to engine. Suggest to create variable derivatives externally.')
        }

        # Check if offset present and fam binomial, Raise warning
        if(fam == "binomial" && !is.Waiver(model$offset)){
          if(getOption('ibis.setupmessages')) myLog('[Estimation]','red','Binomial models fitted with BREG do not support offsets. Offsets were ignored!')
        }
        # -- #

        # If a poisson family is used, weight the observations by their exposure
        if(fam == "poisson"){
          # Get background layer
          bg <- self$get_data("template")
          assertthat::assert_that(!is.na( terra::global(bg, "min", na.rm = TRUE)[,1] ))

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
            # )
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
          model$biodiversity[[1]]$expect <- w * (1/model$biodiversity[[1]]$expect)

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
          model$exposure <- w_full * (1/ unique(model$biodiversity[[1]]$expect)[1])

        } else if(fam == "binomial"){
          # calculating the case weights (equal weights)
          # the order of weights should be the same as presences and backgrounds in the training data
          prNum <- as.numeric(table(model$biodiversity[[1]]$observations[['observed']])["1"]) # number of presences
          bgNum <- as.numeric(table(model$biodiversity[[1]]$observations[['observed']])["0"]) # number of backgrounds
          w <- ifelse(model$biodiversity[[1]]$observations[['observed']] == 1, 1, prNum / bgNum)
          model$biodiversity[[1]]$expect <- w * model$biodiversity[[1]]$expect
          # Convert to numeric
          model$biodiversity[[1]]$observations$observed <- as.numeric( model$biodiversity[[1]]$observations$observed )
        }

        # Check for factors and split them up
        train_cov <- model$biodiversity[[1]]$predictors[,model$biodiversity[[1]]$predictors_names]
        # Check if there any factors, if yes split up
        if(any(model$biodiversity[[1]]$predictors_types$type=='factor')){
          vf <- model$biodiversity[[1]]$predictors_types$predictors[which(model$biodiversity[[1]]$predictors_types$type == "factor")]
          # Get factors
          z <- explode_factor(train_cov[[vf]], name = vf)
          # Remove variables from train_cov and append
          train_cov[[vf]] <- NULL
          train_cov <- cbind(train_cov, z)
          model$biodiversity[[1]]$predictors <- train_cov # Save new in model object
          model$biodiversity[[1]]$predictors_types <- rbind(model$biodiversity[[1]]$predictors_types, data.frame(predictors = colnames(z), type = "numeric"))

          # Also update the formula
          model$biodiversity[[1]]$equation <- stats::update.formula(model$biodiversity[[1]]$equation, paste0(". ~ . -", vf))
          model$biodiversity[[1]]$equation <- stats::update.formula(model$biodiversity[[1]]$equation, paste0(". ~ . +", paste0(colnames(z),collapse = "+")))
        }

        # Prediction container
        pred_cov <- model$predictors[,model$biodiversity[[1]]$predictors_names]
        if(any(model$predictors_types$type=='factor')){
          vf <- model$predictors_types$predictors[which(model$predictors_types$type == "factor")]
          # Get factors
          z <- explode_factor(pred_cov[[vf]], name = vf)
          # Remove variables from train_cov and append
          pred_cov[[vf]] <- NULL
          pred_cov <- cbind(pred_cov, z)
          pred_cov <- pred_cov[,colnames(train_cov)]
          model$predictors <- pred_cov# Save new in model object
          model$predictors_types <- rbind(model$predictors_types, data.frame(predictors = colnames(z), type = "numeric"))
          model$biodiversity[[1]]$predictors_names <- colnames(train_cov)
          model$predictors_names <- colnames(pred_cov)
          assertthat::assert_that(all( colnames(train_cov) %in% colnames(pred_cov) ))
        }
        rm(train_cov, pred_cov)

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
        if(is.Waiver(seed)) { seed <- 1337; settings$set('seed', 1337) }

        # Get output raster
        prediction <- self$get_data('template')

        # Get parameters control
        params <- self$get_data('params')

        # All other needed data for model fitting
        fam <- model$biodiversity[[1]]$family
        li <- model$biodiversity[[1]]$link
        if(!is.null(li)) if(getOption('ibis.setupmessages')) myLog('[Estimation]','red',paste0("Package does not support custom link functions. Ignored!"))
        form <- model$biodiversity[[1]]$equation
        df <- cbind(model$biodiversity[[1]]$predictors,
                    data.frame(observed = model$biodiversity[[1]]$observations[,'observed', drop = TRUE])
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
                                 priors = model$priors,
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
            check_package('scales')
            w <- scales::rescale(w, to = c(1e-6, 1))
            w_full <- scales::rescale(w_full, to = c(1e-6, 1))
          }
          if(anyNA(w)){
            w[is.na(w)] <- 1e-6
            w_full[is.na(w_full)] <- 1e-6
          }
        }

        # Clamp?
        if( settings$get("clamp") ) full <- clamp_predictions(model, full)

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
        # Call garbage collector to save memory
        invisible(gc())

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

          # Make a prediction, but do in parallel so as to not overuse memory
          full$rowid <- 1:nrow(full)
          full_sub <- subset(full, stats::complete.cases(full))
          w_full_sub <- w_full[full_sub$rowid]
          assertthat::assert_that((nrow(full_sub) == length(w_full_sub)) || is.null(w_full_sub) )

          # Tile the problem
          splits <- cut(1:nrow(full_sub), nrow(full_sub) / (min(100, nrow(full_sub) / 10)) )

          # Now depending on parallization setting use foreach
          if(getOption("ibis.runparallel")){
            # Check that future is registered
            if(!foreach::getDoParRegistered()) ibis_future(cores = getOption("ibis.nthread"),
                                                            strategy = getOption("ibis.futurestrategy"))

            # Run the outgoing command
            # out <- foreach::foreach(s = unique(splits),
            #                         .combine = rbind,
            #                         .export = c("splits", "fit_breg", "full_sub",
            #                                     "w_full_sub", "fam", "params"),
            #                         .packages = c("matrixStats"),
            #                         .multicombine = TRUE,
            #                         .inorder = TRUE,
            #                         verbose = settings$get("verbose") ) %do% {
            out <- parallel::mclapply(unique(splits), function(s) {
              i <- which(splits == s)
              # -> external code in utils-boom
              pred_breg <- predict_boom(
                obj = fit_breg,
                newdata = full_sub[i,],
                w = w_full_sub[i],
                fam = fam,
                params = params
              )
              # Summarize the posterior
              preds <- as.data.frame(
                  cbind(
                  matrixStats::rowMeans2(pred_breg, na.rm = TRUE),
                  matrixStats::rowSds(pred_breg, na.rm = TRUE),
                  matrixStats::rowQuantiles(pred_breg, probs = c(.05,.5,.95), na.rm = TRUE),
                  apply(pred_breg, 1, mode)
                )
              )
              names(preds) <- c("mean", "sd", "q05", "q50", "q95", "mode")
              preds$cv <- preds$sd / preds$mean
              return(preds)
              })
            out <- do.call(rbind, out)
          } else {
            out <- data.frame()
            pb <- progress::progress_bar$new(total = length(levels(unique(splits))),
                                             format = "Creating model prediction (:spin) [:bar] :percent")
            for(s in unique(splits)){
              pb$tick()
              i <- which(splits == s)
              # -> external code in utils-boom
              pred_breg <- predict_boom(
                obj = fit_breg,
                newdata = full_sub[i,],
                w = w_full_sub[i],
                fam = fam,
                params = params
              )
              # Summarize the posterior
              preds <- cbind(
                matrixStats::rowMeans2(pred_breg, na.rm = TRUE),
                matrixStats::rowSds(pred_breg, na.rm = TRUE),
                matrixStats::rowQuantiles(pred_breg, probs = c(.05,.5,.95), na.rm = TRUE),
                apply(pred_breg, 1, mode)
              )  |> as.data.frame()
              names(preds) <- c("mean", "sd", "q05", "q50", "q95", "mode")
              preds$cv <- preds$sd / preds$mean
              out <- rbind(out, preds)
              rm(preds, pred_breg)
            }
          }
          assertthat::assert_that(is.data.frame(out), nrow(out)>0,
                                  msg = "Something went wrong withe prediction. Output empty!")
          # Fill output with summaries of the posterior
          stk <- terra::rast()
          for(v in colnames(out)){
            temp <- emptyraster(prediction)
            temp[full_sub$rowid] <- out[,v]
            names(temp) <- v
            suppressWarnings( stk <- c(stk, temp) )
          }
          prediction <- stk;rm(stk)
          prediction <- terra::mask(prediction, self$get_data("template"))
          try({rm(out, full, full_sub)},silent = TRUE)
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
          id = model$id,
          model = model,
          settings = settings,
          fits = list(
            "fit_best" = fit_breg,
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
            # Set type
            if(is.null(type)) type <- self$settings$get("type")
            type <- match.arg(type, c("link", "response"), several.ok = FALSE)
            settings$set("type", type)

            mod <- self$get_data('fit_best')
            model <- self$model
            df <- model$biodiversity[[1]]$predictors
            df <- subset(df, select = attr(mod$terms, "term.labels"))
            w <- model$biodiversity[[1]]$expect # Also get exposure variable

            # Match x.var to argument
            if(is.null(x.var)){
              x.var <- colnames(df)
            } else {
              x.var <- match.arg(x.var, colnames(df), several.ok = TRUE)
            }

            if(is.null(newdata)){
              # Calculate range of predictors
              if(any(model$predictors_types$type=="factor")){
                rr <- sapply(df[model$predictors_types$predictors[model$predictors_types$type=="numeric"]],
                             function(x) range(x, na.rm = TRUE)) |> as.data.frame()
              } else {
                rr <- sapply(df, function(x) range(x, na.rm = TRUE)) |> as.data.frame()
              }
              assertthat::assert_that(nrow(rr)>1, ncol(rr)>=1)

              df_partial <- list()
              if(!is.null(values)){ assertthat::assert_that(length(values) >= 1) }
              # Add all others as constant
              if(is.null(constant)){
                for(n in names(rr)) df_partial[[n]] <- rep( mean(df[[n]], na.rm = TRUE), variable_length )
              } else {
                for(n in names(rr)) df_partial[[n]] <- rep( constant, variable_length )
              }
            } else {
              df_partial <- newdata |> dplyr::select(dplyr::any_of(names(df)))
            }

            # create list to store results
            o <- vector(mode = "list", length = length(x.var))
            names(o) <- x.var

            # loop through x.var
            for(v in x.var) {

              df2 <- df_partial

              if(!is.null(values)){
                df2[[v]] <- values
              } else {
                df2[[v]] <- seq(rr[1, v], rr[2, v], length.out = variable_length)
              }
              df2 <- as.data.frame(df2)

              if(any(model$predictors_types$type=="factor")){
                lvl <- levels(model$predictors[[model$predictors_types$predictors[model$predictors_types$type=="factor"]]])
                df2[model$predictors_types$predictors[model$predictors_types$type=="factor"]] <- factor(lvl[1], levels = lvl)
              }

              # For Integrated model, take the last one
              fam <- model$biodiversity[[length(model$biodiversity)]]$family

              pred_breg <- predict_boom(
                obj = mod,
                newdata = df2,
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
              ) |> as.data.frame()

              names(pred_part) <- c("mean", "sd", "q05", "q50", "q95", "mode")
              pred_part$cv <- pred_part$sd / pred_part$mean
              pred_part <- cbind("partial_effect" = df2[[v]], pred_part)

              # Add variable name for consistency
              pred_part <- cbind("variable" = v, pred_part)

              o[[v]] <- pred_part

            }

            o <- do.call(what = rbind, args = c(o, make.row.names = FALSE))

            if(plot){
              # Make a plot
              g <- ggplot2::ggplot(data = o, ggplot2::aes(x = partial_effect)) +
                ggplot2::theme_classic(base_size = 18) +
                ggplot2::geom_ribbon(ggplot2::aes(ymin = q05, ymax = q95), fill = "grey90") +
                ggplot2::geom_line(ggplot2::aes(y = mean)) +
                ggplot2::facet_wrap(. ~ variable, scales = "free") +
                ggplot2::labs(x = "", y = "Partial effect")
              print(g)
            }
            # Return the data
            return(o)
          },
          # Spatial partial dependence plot
          spartial = function(self, x.var, constant = NULL, newdata = NULL,
                              plot = TRUE, type = NULL){
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

            # Check if newdata is defined, if yes use that one instead
            if(!is.null(newdata)){
              df <- newdata
              assertthat::assert_that(nrow(df) == nrow(model$biodiversity[[1]]$predictors))
            } else {
              # df <- model$biodiversity[[1]]$predictors
              df <- model$predictors
            }
            df <- subset(df, select = attr(mod$terms, "term.labels"))
            w <- model$biodiversity[[1]]$expect # Also get exposure variable

            # Match x.var to argument
            if(is.null(x.var)){
              x.var <- colnames(df)
            } else {
              x.var <- match.arg(x.var, names(df), several.ok = FALSE)
            }

            # Make spatial container for prediction
            template <- model_to_background(model)

            # Add all others as constant
            if(is.null(constant)){
              for(n in names(df)) if(n != x.var) df[[n]] <- suppressWarnings( mean(model$predictors[[n]], na.rm = TRUE) )
            } else {
              for(n in names(df)) if(n != x.var) df[[n]] <- constant
            }

            if(any(model$predictors_types$type=="factor")){
              lvl <- levels(model$predictors[[model$predictors_types$predictors[model$predictors_types$type=="factor"]]])
              df[[model$predictors_types$predictors[model$predictors_types$type=="factor"]]] <-
                factor(lvl[1], levels = lvl)
              # FIXME: Assigning the first level (usually reference) for now. But ideally find a way to skip factors from partial predictions
            }

            # For Integrated model, take the last one
            fam <- model$biodiversity[[length(model$biodiversity)]]$family

            pred_breg <- predict_boom(
              obj = mod,
              newdata = df,
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
            ) |> as.data.frame()
            names(pred_part) <- c("mean", "sd", "q05", "q50", "q95", "mode")
            pred_part$cv <- pred_part$sd / pred_part$mean

            # Now create spatial prediction
            template <- fill_rasters(pred_part, template)

            # Do plot and return result
            if(plot){
              terra::plot(template, col = ibis_colours$ohsu_palette,
                   main = paste0("Spartial effect of ", x.var, collapse = ","))
            }
            return(template)
          },
          # Get coefficients from breg
          get_coefficients = function(self){
            # Returns a vector of the coefficients with direction/importance
            obj <- self$get_data("fit_best")
            cofs <- posterior::summarise_draws(obj$beta)
            cofs <- subset(cofs, select = c("variable", "mean", "sd"))
            names(cofs) <- c("Feature", "Beta", "Sigma")
            # Remove intercept(s)
            int <- grep("Intercept",cofs$Feature,ignore.case = TRUE)
            if(length(int)>0) cofs <- cofs[-int,]
            return(cofs)
          },
          # Model convergence check
          has_converged = function(self){
            fit <- self$get_data("fit_best")
            if(is.Waiver(fit)) return(FALSE)
            return(TRUE)
          },
          # Residual function
          get_residuals = function(self){
            # Get best object
            obj <- self$get_data("fit_best")
            if(is.Waiver(obj)) return(obj)
            # Get residuals
            rd <- obj$deviance.residuals
            assertthat::assert_that(length(rd)>0)
            return(rd)
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
            df <- newdata
            df <- subset(df, select = attr(mod$terms, "term.labels"))

            # Clamp?
            if( settings$get("clamp") ) df <- clamp_predictions(model, df)

            if(!is.Waiver(settings$get('bias_variable'))){
              for(i in 1:length(settings$get('bias_variable'))){
                if(settings$get('bias_variable')[i] %notin% colnames(df)){
                  if(getOption('ibis.setupmessages')) myLog('[Estimation]','red','Did not find bias variable in prediction object!')
                  next()
                }
                df[,settings$get('bias_variable')[i]] <- settings$get('bias_value')[i]
              }
            }

            df$rowid <- 1:nrow(df)
            df_sub <- subset(df, stats::complete.cases(df))
            w <- model$biodiversity[[1]]$expect # Also get exposure variable

            # For Integrated model, take the last one
            fam <- model$biodiversity[[length(model$biodiversity)]]$family

            # Rather predict in steps than for the whole thing
            out <- data.frame()

            # Tile the problem
            splits <- cut(1:nrow(df_sub), nrow(df_sub) / (min(100, nrow(df_sub) / 10)) )

            pb <- progress::progress_bar$new(total = length(levels(unique(splits))),
                                             format = "Projecting on new data (:spin) [:bar] :percent")
            for(s in unique(splits)){
              pb$tick()
              i <- which(splits == s)
              # -> external code in utils-boom
              pred_breg <- predict_boom(
                obj = mod,
                newdata = df_sub[i,],
                w = unique(w)[2],
                fam = fam,
                params = settings$data
              )
              # Summarize the posterior
              preds <- cbind(
                matrixStats::rowMeans2(pred_breg, na.rm = TRUE),
                matrixStats::rowSds(pred_breg, na.rm = TRUE),
                matrixStats::rowQuantiles(pred_breg, probs = c(.05,.5,.95), na.rm = TRUE),
                apply(pred_breg, 1, mode)
              )  |> as.data.frame()
              names(preds) <- c("mean", "sd", "q05", "q50", "q95", "mode")
              preds$cv <- preds$sd / preds$mean
              out <- rbind(out, preds)
              rm(preds, pred_breg)
            }

            # Now create spatial prediction
            prediction <- try({emptyraster( self$model$predictors_object$get_data()[[1]] )},silent = TRUE) # Background
            if(inherits(prediction, "try-error")){
              prediction <- terra::rast(self$model$predictors[,c("x", "y")], crs = terra::crs(model$background),type = "xyz") |>
                emptyraster()
            }
            prediction[df_sub$rowid] <- out[,layer]
            names(prediction) <- layer

            return(prediction)
          }
        )
        return(out)
      }
    )
  ) # End of bdproto object
} # End of function
