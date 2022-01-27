#' @include bdproto-engine.R utils-spatial.R bdproto-distributionmodel.R
NULL

#' Engine for use of Bayesian Additive Regression Trees (BART)
#'
#' @description Bayesian regression approach to a sum of complementary trees by shrinking
#' the fit of each tree using a regularization prior.
#' BART models have been shown to compare favourable among machine learning algorithms (Dorie et al. 2019)
#' Default prior preference is for trees to be small (few terminal nodes) and shrinkage towards 0.
#' Prior distributions can furthermore be set for:
#' - probability that a tree stops at a node of a given depth
#' - probability that a given variable is chosen for a splitting rule
#' - probability of splitting that variable at a particular value
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param nburn A [`numeric`] estimate of the burn in samples
#' @param ntree A [`numeric`] estimate of the number of trees to be used in the sum-of-trees formulation.
#' @param chains A number of the number of chains to be used (Default: 4)
#' @references  Carlson, CJ. embarcadero: Species distribution modelling with Bayesian additive regression trees in r. Methods Ecol Evol. 2020; 11: 850– 858. https://doi.org/10.1111/2041-210X.13389
#' @references  Vincent Dorie (2020). dbarts: Discrete Bayesian Additive Regression Trees Sampler. R package version 0.9-19. https://CRAN.R-project.org/package=dbarts
#' @name engine_bart
NULL
#' @rdname engine_bart
#' @export

engine_bart <- function(x,
                        nburn = 250,
                        ntree = 1000,
                        chains = 4,
                       ...) {

  # Check whether dbarts package is available
  check_package('dbarts')
  if(!("dbarts" %in% loadedNamespaces()) || ('dbarts' %notin% sessionInfo()$otherPkgs) ) {
    try({requireNamespace('dbarts');attachNamespace("dbarts")},silent = TRUE)
  }

  # assert that arguments are valid
  assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                          inherits(x$background,'sf'),
                          is.numeric(nburn),
                          is.numeric(ntree),
                          is.numeric(chains)
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

  # Set up dbarts control with some parameters, rest default
  dc <- dbarts::dbartsControl(keepTrees	= TRUE, # Keep trees
                              n.burn = nburn,
                              n.trees = ntree,
                              n.chains = chains,
                              n.threads = ifelse( dbarts::guessNumCores() < getOption('ibis.nthread'),dbarts::guessNumCores(),getOption('ibis.nthread'))
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
        'dc' = dc
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
                             nburn = 250,
                             ntree = 1000,
                             chains = 4,
                             cores = dbarts::guessNumCores(),
                             verbose = TRUE,
                             ...
      ){
        # Set up boosting control
        dc <- dbarts::dbartsControl(verbose = verbose,
                                    n.burn = nburn,
                                    n.trees = ntree,
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
        # Messager
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Engine setup.')

        # Add pseudo-absence points if necessary
        if('poipo' == model$biodiversity[[1]]$type && model$biodiversity[[1]]$family == 'poisson') {

          # Get background layer
          bg <- self$get_data('template')
          assertthat::assert_that(!is.na(cellStats(bg,min)))

          abs <- create_pseudoabsence(
            env = model$predictors,
            presence = model$biodiversity[[1]]$observations,
            bias = settings$get('bias_variable'),
            template = bg,
            npoints = ifelse(ncell(bg)<10000,ncell(bg),10000),
            replace = TRUE
          )
          abs$intercept <- 1 # Redundant for this engine
          # Combine absence and presence and save
          abs_observations <- abs[,c('x','y')]; abs_observations[['observed']] <- 0

          # Rasterize observed presences
          pres <- raster::rasterize(model$biodiversity[[1]]$predictors[,c('x','y')],
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
          envs <- get_ngbvalue(coords = obs[,c('x','y')],
                               env =  model$predictors[,c("x","y", model[['predictors_names']])],
                               longlat = raster::isLonLat(self$get_data("template")),
                               field_space = c('x','y')
          )
          envs$intercept <- 1
          # Overwrite observations
          model$biodiversity[[1]]$observations <- rbind(obs, abs_observations)

          # Format out
          df <- rbind(envs,
                      abs[,c('x','y','intercept', model$biodiversity[[1]]$predictors_names)]) %>%
            subset(., complete.cases(.) )

          # Preprocessing security checks
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
          df$w <- w # Also add as column

          model$biodiversity[[1]]$predictors <- df
          model$biodiversity[[1]]$expect <- w
        } else {
          # If family is not poisson, assume factor distribution for response
          assertthat::assert_that(  length( unique(model$biodiversity[[1]]$observations[['observed']])) == 2)
          # calculating the case weights (equal weights)
          # the order of weights should be the same as presences and backgrounds in the training data
          prNum <- as.numeric(table(model$biodiversity[[1]]$observations[['observed']])["1"]) # number of presences
          bgNum <- as.numeric(table(model$biodiversity[[1]]$observations[['observed']])["0"]) # number of backgrounds
          w <- ifelse(model$biodiversity[[1]]$observations[['observed']] == 1, 1, prNum / bgNum)
          model$biodiversity[[1]]$expect <- w

          model$biodiversity[[1]]$observations[['observed']] <- factor(model$biodiversity[[1]]$observations[['observed']])
        }

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
        # Messager
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Starting fitting...')

        # Get output raster
        prediction <- self$get_data('template')

        # Get dbarts control
        dc <- self$get_data('dc')

        # All other needed data for model fitting
        equation <- model$biodiversity[[1]]$equation
        data <- cbind(model$biodiversity[[1]]$predictors, data.frame(observed = model$biodiversity[[1]]$observations[,'observed']) )
        # Subset to predictor names
        data <- subset(data, select = c('observed', model$biodiversity[[1]]$predictors_names) )
        if(model$biodiversity[[1]]$family=='binomial') data$observed <- factor(data$observed)
        w <- model$biodiversity[[1]]$expect # The expected weight
        data$w <- w # Also add as predictor
        full <- model$predictors # All predictors

        # Select predictors
        full <- subset(full, select = c('x','y',model$biodiversity[[1]]$predictors_names))
        full$cellid <- rownames(full) # Add rownames
        full <- subset(full, complete.cases(full))

        assertthat::assert_that(
          is.null(w) || length(w) == nrow(data),
          is.formula(equation),
          all( model$biodiversity[[1]]$predictors_names %in% names(full) )
        )

        if('offset' %in% names(model$biodiversity[[1]]) ){
          # Add offset to full prediction and load vector
          # TBD
        } else { off = NULL }

        # --- #
        # Parameter tuning #
        if(settings$get('varsel') == "reg"){
          if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Starting hyperparameters search.')

          cv_bart <- dbarts::xbart(
            formula = equation, data = data,
            n.samples = round( nrow(data) * 0.1 ), # Draw posterior samples (10% of dataset)
            n.test = 5, # Number of folds
            method = "k-fold",
            n.reps = 4L, # For replications
            control = dc,
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
                                   keeptrees = dc@keepTrees,
                                   # Hyper parameters
                                   k = k, power = power, base =  base,
                                   ntree = dc@n.trees,
                                   nthread = dc@n.threads,
                                   nchain = dc@n.chains,
                                   nskip = dc@n.burn,
                                   verbose = settings$get('verbose')
          )
        } else {
          fit_bart <- dbarts::bart(equation,
                                   data,
                                   keeptrees = dc@keepTrees,
                                   weights = w,
                                   ntree = dc@n.trees,
                                   # Hyper parameters
                                   k = k, power = power, base =  base,
                                   nthread = dc@n.threads,
                                   nchain = dc@n.chains,
                                   nskip = dc@n.burn,
                                   verbose = settings$get('verbose')
          )
        }

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
            pred_bart <- dbarts:::predict.bart(object = fit_bart,
                                        newdata = full,
                                        type = 'response')
          )
          # Summarize quantiles and sd from posterior
          ms <- as.data.frame(
                 cbind( apply(pred_bart, 2, mean),
                        matrixStats::colSds(pred_bart),
                        matrixStats::colQuantiles(pred_bart, probs = c(.05,.5,.95)),
                        apply(pred_bart, 2, mode)
                       )
          )
          names(ms) <- c("mean","sd", "q05", "q50", "q95", "mode")
          ms$cv <- ms$mean / ms$sd

          # Add them
          prediction <- raster::stack()
          for(post in names(ms)){
            prediction2 <- self$get_data('template')
            prediction2[as.numeric(full$cellid)] <- ms[[post]]; names(prediction2) <- post
            prediction <- raster::addLayer(prediction, prediction2)
            rm(prediction2)
          }

        } else {
          # No prediction done
          prediction <- NULL
        }
        # Compute end of computation time
        settings$set('end.time', Sys.time())
        # Also append boosting control option to settings
        for(entry in slotNames(dc)) settings$set(entry, slot(dc,entry))

        # Create output
        out <- bdproto(
          "BART-Model",
          DistributionModel,
          id = new_id(),
          model = model,
          settings = settings,
          fits = list(
            "fit_best" = fit_bart,
            "prediction" = prediction
          ),
          # Partial effects
          partial = function(self, x.var = NULL, ...){
            model <- self$get_data('fit_best')
            assertthat::assert_that(x.var %in% attr(model$fit$data@x,'term.labels') || is.null(x.var),
                                    msg = 'Variable not in predicted model' )
            # Check if family is binomial, if so alter
            not_binomial = self$get_data('model')$biodiversity[[1]]$family != 'binomial'
            bart_partial_effect(model, x.var = x.var, transform = not_binomial, ... )
          },
          # Spatial partial dependence plot option from embercardo
          spartial = function(self, predictors, x.var = NULL, equal = FALSE, smooth = 1, transform = TRUE){
            model <- self$get_data('fit_best')
            assertthat::assert_that(x.var %in% attr(model$fit$data@x,'term.labels'),
                                    msg = 'Variable not in predicted model' )

            if( self$model$biodiversity[[1]]$family != 'binomial' && transform) warning('Check whether transform should not be set to False!')

            # Calculate
            p <- bart_partial_space(model, predictors, x.var, equal, smooth, transform)

            raster::plot(p, col = ibis_colours$viridis_plasma, main = paste0(x.var, collapse ='|'))
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
      }
    )
  ) # End of bdproto object
} # End of function
