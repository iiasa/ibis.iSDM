#' @include bdproto-engine.R utils-spatial.R bdproto-distributionmodel.R
NULL

#' Engine for use of Bayesian Additive Regression Trees (BART)
#'
#' @description The Bayesian regression approach to a sum of complementary trees
#' is to shrink the said fit of each tree through a regularization prior. BART
#' models provide non-linear highly flexible estimation and have been shown to
#' compare favourable among machine learning algorithms (Dorie et al. 2019).
#' Default prior preference is for trees to be small (few terminal nodes) and
#' shrinkage towards \code{0}.
#'
#' This package requires the \code{"dbarts"} R-package to be installed. Many
#' of the functionalities of this engine have been inspired by the
#' \code{"embarcadero"} R-package. Users are therefore advised to cite if they
#' make heavy use of BART.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param iter A [`numeric`] estimate of the number of trees to be used in the
#' sum-of-trees formulation (Default: \code{1000}).
#' @param nburn A [`numeric`] estimate of the burn in samples (Default: \code{250}).
#' @param chains A number of the number of chains to be used (Default: \code{4}).
#' @param type The mode used for creating posterior predictions. Either \code{"link"}
#' or \code{"response"} (Default: \code{"response"}).
#' @param ... Other options.
#'
#' @details Prior distributions can furthermore be set for:
#' * probability that a tree stops at a node of a given depth (Not yet implemented)
#' * probability that a given variable is chosen for a splitting rule
#' * probability of splitting that variable at a particular value (Not yet implemented)
#'
#' @returns An [Engine].
#'
#' @references
#' * Carlson, CJ. embarcadero: Species distribution modelling with Bayesian additive
#' regression trees in r. Methods Ecol Evol. 2020; 11: 850– 858. https://doi.org/10.1111/2041-210X.13389
#' * Dorie, V., Hill, J., Shalit, U., Scott, M., & Cervone, D. (2019). Automated
#' versus do-it-yourself methods for causal inference: Lessons learned from a data analysis competition. Statistical Science, 34(1), 43-68.
#' * Vincent Dorie (2020). dbarts: Discrete Bayesian Additive Regression Trees Sampler.
#' R package version 0.9-19. https://CRAN.R-project.org/package=dbarts
#'
#' @family engine
#'
#' @examples
#' \dontrun{
#' # Add BART as an engine
#' x <- distribution(background) |> engine_bart(iter = 100)
#' }
#'
#' @name engine_bart
NULL

#' @rdname engine_bart
#' @export
engine_bart <- function(x,
                        iter = 1000,
                        nburn = 250,
                        chains = 4,
                        type = "response",
                       ...) {

  # Check whether dbarts package is available
  check_package('dbarts')
  if(!("dbarts" %in% loadedNamespaces()) || ('dbarts' %notin% utils::sessionInfo()$otherPkgs) ) {
    try({requireNamespace('dbarts');attachNamespace("dbarts")},silent = TRUE)
  }

  # assert that arguments are valid
  assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                          inherits(x$background,'sf'),
                          is.numeric(nburn),
                          is.numeric(iter),
                          is.character(type),
                          is.numeric(chains)
                          )
  type <- match.arg(type, choices = c("link", "predictor", "response", "ppd"), several.ok = FALSE)
  if(type == "predictor") type <- "link"
  if(nburn > iter) nburn <- floor( iter / 4)

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

  # Set up dbarts control with some parameters, rest default
  dc <- dbarts::dbartsControl(keepTrees	= TRUE, # Keep trees
                              n.burn = nburn,
                              n.trees = iter,
                              n.chains = chains,
                              n.threads = ifelse( dbarts::guessNumCores() < getOption('ibis.nthread'),dbarts::guessNumCores(),getOption('ibis.nthread'))
  )
  # Other parameters
  # Set up the parameter list
  params <- list(
    type = type,
    ...
  )

  # Print a message in case there is already an engine object
  if(!is.Waiver(x$engine)) myLog('[Setup]','yellow','Replacing currently selected engine.')

  # Set engine in distribution object
  x$set_engine(
    bdproto(
      "BART-Engine",
      Engine,
      name = "<BART>",
      data = list(
        'template' = template,
        'dc' = dc,
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
                             iter = 1000,
                             nburn = 250,
                             chains = 4,
                             cores = dbarts::guessNumCores(),
                             verbose = TRUE,
                             ...
      ){
        # Set up boosting control
        dc <- dbarts::dbartsControl(verbose = verbose,
                                    n.burn = nburn,
                                    n.trees = iter,
                                    n.chains = chains,
                                    n.threads = cores,
                                    ...
        )
        # Overwrite existing
        self$data$dc <- dc
      },
      # Setup function
      setup = function(self, model, settings = NULL, ...){
        # Simple security checks
        assertthat::assert_that(
          assertthat::has_name(model, 'background'),
          assertthat::has_name(model, 'biodiversity'),
          inherits(settings,'Settings') || is.null(settings),
          nrow(model$predictors) == ncell(self$get_data('template')),
          length(model$biodiversity) == 1 # Only works with single likelihood. To be processed separately
        )
        # Messenger
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Engine setup.')

        # Add pseudo-absence points if necessary
        if('poipo' == model$biodiversity[[1]]$type && model$biodiversity[[1]]$family == 'poisson') {
          # Warning since PPMs are not really performing / correctly set up in bart
          if(getOption('ibis.setupmessages')) myLog('[Estimation]','yellow','Engine BART prone to overfit Poisson-distributed occurrence data.\nConsider non-linear xgboost as alternative!')

          # Get background layer
          bg <- self$get_data('template')
          assertthat::assert_that(!is.na(terra::global(bg, "min", na.rm = TRUE)[,1]))

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
                           weight = 1 # Set those to 1 so that absences become ratio of pres/abs
          )
          df$w <- w # Also add as column

          model$biodiversity[[1]]$predictors <- df
          model$biodiversity[[1]]$expect <- w * (1/model$biodiversity[[1]]$expect) # Multiply with provided weights
        } else {
          # If family is not poisson, assume factor distribution for response
          assertthat::assert_that(  length( unique(model$biodiversity[[1]]$observations[['observed']])) == 2)
          # calculating the case weights (equal weights)
          # the order of weights should be the same as presences and backgrounds in the training data
          prNum <- as.numeric(table(model$biodiversity[[1]]$observations[['observed']])["1"]) # number of presences
          bgNum <- as.numeric(table(model$biodiversity[[1]]$observations[['observed']])["0"]) # number of backgrounds
          w <- ifelse(model$biodiversity[[1]]$observations[['observed']] == 1, 1, prNum / bgNum)
          model$biodiversity[[1]]$expect <- w * model$biodiversity[[1]]$expect # Multiply with provided weights
          model$biodiversity[[1]]$observations[['observed']] <- factor(model$biodiversity[[1]]$observations[['observed']])
        }

        # Split up factors as this is done anyway during the fitting!
        # Check for factors and split them up
        train_cov <- model$biodiversity[[1]]$predictors[, c(model$biodiversity[[1]]$predictors_names)]
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
        pred_cov <- model$predictors[,c('x','y',model$biodiversity[[1]]$predictors_names)]
        if(any(model$predictors_types$type=='factor')){
          vf <- model$predictors_types$predictors[which(model$predictors_types$type == "factor")]
          # Get factors
          z <- explode_factor(pred_cov[[vf]], name = vf)
          # Remove variables from train_cov and append
          pred_cov[[vf]] <- NULL
          pred_cov <- cbind(pred_cov, z)
          pred_cov <- pred_cov[,c("x", "y", colnames(train_cov))]
          model$predictors <- pred_cov # Save new in model object
          model$predictors_types <- rbind(model$predictors_types, data.frame(predictors = colnames(z), type = "numeric"))
          model$biodiversity[[1]]$predictors_names <- colnames(train_cov)
          model$predictors_names <- colnames(pred_cov)
          assertthat::assert_that(all( colnames(train_cov) %in% colnames(pred_cov) ))
        }
        rm(train_cov, pred_cov)

        # Process and add priors if set
        params <- self$get_data("params")
        if(!is.Waiver(model$priors)){
          assertthat::assert_that(
            all( model$priors$varnames() %in% model$biodiversity[[1]]$predictors_names )
          )
          # Match position of variables with monotonic constrains
          mc <- rep(1 / length( model$biodiversity[[1]]$predictors_names ), length( model$biodiversity[[1]]$predictors_names ) )
          names(mc) <- model$biodiversity[[1]]$predictors_names
          for(v in model$priors$varnames()){
            mc[v] <- model$priors$get(v)
          }
          # Save the priors in the model parameters
          params[["priors"]] <- mc
        } else { params[["priors"]] <- new_waiver() }
        self$set_data("params", params)

        # Instead of invisible return the model object
        return( model )
      },
      # Training function
      train = function(self, model, settings, ...){
        assertthat::assert_that(
          inherits(settings,'Settings'),
          is.list(model),length(model) > 1,
          # Check that model id and setting id are identical
          settings$modelid == model$id
        )
        # Get name
        name <- model$biodiversity[[1]]$name

        # Messenger
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','green',paste0( 'Starting fitting: ', name))

        # Get output raster
        prediction <- self$get_data('template')

        # Get dbarts control and params
        dc <- self$get_data('dc')
        params <- self$get_data('params')

        # All other needed data for model fitting
        equation <- model$biodiversity[[1]]$equation
        data <- cbind(model$biodiversity[[1]]$predictors,
                      data.frame(observed = model$biodiversity[[1]]$observations[,'observed', drop = TRUE]) )
        # Subset to predictor names
        data <- subset(data, select = c('observed', model$biodiversity[[1]]$predictors_names) )
        if(model$biodiversity[[1]]$family=='binomial') data$observed <- factor(data$observed)
        w <- model$biodiversity[[1]]$expect # The expected weight
        full <- model$predictors # All predictors

        # Select predictors
        full <- subset(full, select = c('x','y', model$biodiversity[[1]]$predictors_names))
        full$cellid <- rownames(full) # Add rownames
        full <- subset(full, stats::complete.cases(full))

        # Clamp?
        if( settings$get("clamp") ) full <- clamp_predictions(model, full)

        assertthat::assert_that(
          is.null(w) || length(w) == nrow(data),
          is.formula(equation),
          all( model$biodiversity[[1]]$predictors_names %in% names(full) )
        )

        if(!is.Waiver(model$offset)){
          # Add offset to full prediction and load vector
          if(model$biodiversity[[1]]$family == "poisson"){
            # Offsets are only supported for binary dbarts models, but maybe there is an option
            if(getOption('ibis.setupmessages')) myLog('[Estimation]','yellow','Offsets are not supported for poisson models. Trying to modify weights.')
            w <- w + model$biodiversity[[1]]$offset[,"spatial_offset"]
            # Check and correct for issues
            if(any(w < 0, na.rm = TRUE)) {
              w <- scales::rescale(w, to = c(1e-6, 1))
            }
            if(anyNA(w)){
              w[is.na(w)] <- 1e-6
            }
          } else if(model$biodiversity[[1]]$family == "binomial"){
            # Set the created ranges and binaryOffset
            off <- model$biodiversity[[1]]$offset[,"spatial_offset"]
          }
        } else { off = 0.0 }

        # Specify splitprobs depending on whether priors have been set
        if( !is.Waiver(params[["priors"]]) ){
          splitprobs = params[["priors"]]
        } else { splitprobs = NULL }

        # --- #
        # Parameter tuning #
        if(settings$get('optim_hyperparam')){
          if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Starting hyperparameters search.')

          cv_bart <- dbarts::xbart(
            formula = equation, data = data,
            n.samples = round( nrow(data) * 0.1 ), # Draw posterior samples (10% of dataset)
            n.test = 5, # Number of folds
            method = "k-fold",
            n.reps = 4L, # For replications
            control = dc,
            offset = off,
            loss = ifelse(is.factor(data$observed), "log", "rmse"),
            n.trees = dc@n.trees,
            k = c(1, 2, 4), # Prior for node-mean SD
            power = c(1.5, 2), # Prior growth probability
            base = c(0.75, 0.8, 0.95), # Tree growth probability
            drop = TRUE, # Drop those with only one record
            n.threads = dc@n.threads,
            verbose = settings$get('verbose')
          )
          # An array of dimensions n.reps * length(n.trees) * length(k) * length(power) * length(base)
          # Convert to data.frame
          cv_bart <- as.data.frame.table(cv_bart)
          best <- which.min(cv_bart$Freq) # Get the setting with lowest loss/error
          k <- as.numeric( as.character(cv_bart$k[best]) )
          power <- as.numeric( as.character(cv_bart$power[best]) )
          base <- as.numeric( as.character(cv_bart$base[best]) )

        } else {
          # Pick default hyperparameters for bart
          # k = 2 implies that the maximum and minimum are each approximately
          # 2 standard deviations from the mean, or ≈ 95 percent prior probability
          # in the interval (ymin, ymax) when Y is continuous and symmetrically distributed.
          # Source: https://journals.sagepub.com/doi/10.1177/2378023119825886
          k <- 2.0; power = 2.0 ; base = 0.95
        }
        # --- #
        # Fit the model. Little hack to work correctly with binomials...
        if(is.factor(data$observed)){
          fit_bart <- dbarts::bart(y.train = data[,'observed'],
                                   x.train = data[,model$biodiversity[[1]]$predictors_names],
                                   # To make partial plots faster
                                   keeptrees = dc@keepTrees,
                                   keepevery = 10,
                                   # weights = w,
                                   binaryOffset = off,
                                   # Hyper parameters
                                   k = k, power = power, base =  base,
                                   splitprobs = splitprobs,
                                   ntree = dc@n.trees,
                                   nthread = dc@n.threads,
                                   nchain = dc@n.chains,
                                   nskip = dc@n.burn,
                                   verbose = settings$get('verbose')
          )
        } else {
          fit_bart <- dbarts::bart(y.train = data[,'observed'],
                                   x.train = data[,model$biodiversity[[1]]$predictors_names],
                                   # To make partial plots faster
                                   keeptrees = dc@keepTrees,
                                   keepevery = 10,
                                   weights = w,
                                   ntree = dc@n.trees,
                                   # Hyper parameters
                                   k = k, power = power, base =  base,
                                   splitprobs = splitprobs,
                                   nthread = dc@n.threads,
                                   nchain = dc@n.chains,
                                   nskip = dc@n.burn,
                                   verbose = settings$get('verbose')
          )
        }

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
          check_package("foreach")
          params <- self$get_data("params")

          full$rowid <- 1:nrow(full)

          # Tile the problem
          splits <- cut(1:nrow(full), nrow(full) / min(nrow(full) / 4, 5000) )

          # Get offset if existing
          if(is.Waiver(model$offset)) of <- NULL else of <- scales::rescale(model$offset[full$cellid, "spatial_offset"], to = c(1e-6, 1))

          # Make a prediction
          ms <- foreach::foreach(s = unique(splits),
                                 .inorder = TRUE,
                                 .combine = rbind,
                                 .errorhandling = "stop",
                                 .multicombine = TRUE,
                                 .export = c("splits", "fit_bart", "full", "model", "params", "of"),
                                 .packages = c("dbarts", "matrixStats")) %do% {
                                   i <- which(splits == s)

                                   pred_bart <- predict(object = fit_bart,
                                                                      newdata = full[i, model$biodiversity[[1]]$predictors_names],
                                                                      type = params$type,
                                                                      offset = of[i]
                                   )
                                   # Summarize quantiles and sd from posterior
                                   ms <- as.data.frame(
                                     cbind( apply(pred_bart, 2, function(x) mean(x, na.rm = TRUE)),
                                            matrixStats::colSds(pred_bart),
                                            matrixStats::colQuantiles(pred_bart, probs = c(.05,.5,.95)),
                                            apply(pred_bart, 2, mode)
                                     )
                                   )
                                   names(ms) <- c("mean","sd", "q05", "q50", "q95", "mode")
                                   ms$cv <- ms$sd / ms$mean
                                   rm(pred_bart)
                                   return( ms )
                                 } # End of processing
          assertthat::assert_that(nrow(ms)>0,
                                  nrow(ms) == nrow(full))

          # Add them through a loop since the cellid changed
          prediction <- terra::rast()
          for(post in names(ms)){
            prediction2 <- self$get_data('template')
            prediction2[as.numeric(full$cellid)] <- ms[[post]]; names(prediction2) <- post
            suppressWarnings( prediction <- c(prediction, prediction2) )
            rm(prediction2)
          }
          # plot(prediction$mean, col = ibis_colours$sdm_colour)
          try({rm(ms, full)},silent = TRUE)
        } else {
          # No prediction done
          prediction <- NULL
        }
        # Compute end of computation time
        settings$set('end.time', Sys.time())
        # Also append boosting control option to settings
        for(entry in methods::slotNames(dc)) settings$set(entry, methods::slot(dc,entry))
        for(entry in names(params)) settings$set(entry, params[[entry]])
        # Create output
        out <- bdproto(
          "BART-Model",
          DistributionModel,
          id = model$id,
          model = model,
          settings = settings,
          fits = list(
            "fit_best" = fit_bart,
            "params" = params,
            "fit_best_equation" = equation,
            "prediction" = prediction
          ),
          # Partial effects
          partial = function(self, x.var = NULL, constant = NULL, variable_length = 100,
                             values = NULL, newdata = NULL, plot = FALSE, type = NULL, ...){

            model <- self$get_data('fit_best')

            assertthat::assert_that(all(x.var %in% attr(model$fit$data@x,'term.labels')) || is.null(x.var),
                                    msg = 'Variable not in predicted model' )

            # Match x.var to argument
            if(is.null(x.var)){
              x.var <- attr(model$fit$data@x,'term.labels')
            } else {
              x.var <- match.arg(x.var, attr(model$fit$data@x,'term.labels'), several.ok = TRUE)
            }

            if(is.null(newdata)){
              bart_partial_effect(model, x.var = x.var,
                                  transform = self$settings$data$binary,
                                  variable_length = variable_length,
                                  values = values,
                                  equal = TRUE,
                                  plot = plot )
            } else {
              # Set the values to newdata
              bart_partial_effect(model, x.var = x.var,
                                  transform = self$settings$data$binary,
                                  values = newdata[[x.var]],
                                  plot = plot)
            }
          },
          # Spatial partial dependence plot option from embercardo
          spartial = function(self, x.var = NULL, newdata = NULL, equal = FALSE,
                              smooth = 1, transform = TRUE, type = NULL){
            fit <- self$get_data('fit_best')
            model <- self$model
            if(is.null(newdata)){
              predictors <- model$predictors_object$get_data()
            } else {
              predictors <- newdata
              assertthat::assert_that(all(x.var %in% colnames(predictors)),
                                      msg = 'Variable not in provided data!')
            }
            assertthat::assert_that(all(x.var %in% attr(fit$fit$data@x,'term.labels')),
                                    msg = 'Variable not in predicted model' )

            if( model$biodiversity[[1]]$family != 'binomial' && transform) warning('Check whether transform should not be set to False!')

            # Calculate
            p <- bart_partial_space(fit, predictors, x.var, equal, smooth, transform)

            terra::plot(p, col = ibis_colours$ohsu_palette,
                        main = paste0(x.var, collapse ='|'))
            # Also return spatial
            return(p)
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
            rd <- residuals(obj)
            if(length(rd)==0) rd <- new_waiver()
            return(rd)
          },
          # Coefficient function
          get_coefficients = function(self){
            # Returns a vector of the coefficients with direction/importance
            cofs <- self$summary()
            cofs$Sigma <- NA
            names(cofs) <- c("Feature", "Weights")
            return(cofs)
          },
          # Engine-specific projection function
          project = function(self, newdata, type = "response", layer = 'mean'){
            assertthat::assert_that(!missing(newdata),
                                    is.data.frame(newdata))

            # get model data
            model <- self$model

            # Define rowids as those with no missing data
            rownames(newdata) <- 1:nrow(newdata)
            newdata$rowid <- as.numeric( rownames(newdata) )
            newdata <- subset(newdata, stats::complete.cases(newdata))

            # Also get settings for bias values
            settings <- self$settings

            # Clamp?
            if( settings$get("clamp") ) newdata <- clamp_predictions(model, newdata)

            if(!is.Waiver(settings$get('bias_variable'))){
              for(i in 1:length(settings$get('bias_variable'))){
                if(settings$get('bias_variable')[i] %notin% names(newdata)){
                  if(getOption('ibis.setupmessages')) myLog('[Estimation]','red','Did not find bias variable in prediction object!')
                  next()
                }
                newdata[[settings$get('bias_variable')[i]]] <- settings$get('bias_value')[i]
              }
            }

            # Make a prediction
            suppressWarnings(
              pred_bart <- predict(object = self$get_data('fit_best'),
                                                 newdata = newdata,
                                                 type = type) |> t()
              )
            assertthat::assert_that(nrow(pred_bart) == nrow(newdata))
            # Fill output with summaries of the posterior
            prediction <- try({emptyraster( model$predictors_object$get_data()[[1]] )},silent = TRUE) # Background
            if(inherits(prediction, "try-error")){
              prediction <- terra::rast(self$model$predictors[,c("x", "y")], crs = terra::crs(model$background),type = "xyz") |>
                emptyraster()
            }

            if(layer == "mean"){
              prediction[newdata$rowid] <- matrixStats::rowMeans2(pred_bart)
            } else if(layer == "sd"){
              prediction[newdata$rowid] <- matrixStats::rowSds(pred_bart)
            } else if(layer == "q05"){
              prediction[newdata$rowid] <- matrixStats::rowQuantiles(pred_bart, probs = c(.05))
            } else if(layer == "q50" || layer == "median"){
              prediction[newdata$rowid] <- matrixStats::rowQuantiles(pred_bart, probs = c(.5))
            } else if(layer == "q95"){
              prediction[newdata$rowid] <- matrixStats::rowQuantiles(pred_bart, probs = c(.95))
            } else if(layer == "mode"){
              prediction[newdata$rowid] <- apply(pred_bart, 1, mode)
            } else if(layer == "cv"){
              prediction[newdata$rowid] <- matrixStats::rowSds(pred_bart) / matrixStats::rowMeans2(pred_bart)
            } else { message("Custom posterior summary not yet implemented.")}
            return(prediction)
          }
        )
        return(out)
      }
    )
  ) # End of bdproto object
} # End of function
