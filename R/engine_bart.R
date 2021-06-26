#' @include bdproto-engine.R utils-spatial.R bdproto-distributionmodel.R
NULL

#' Engine for use of Bayesian Additive Regression Trees (BART)
#'
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
  if(!is.Waiver(x$engine)) message('Replacing currently selected engine.')

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
      # Dummy function for spatial latent effects (not yet implemented)
      calc_latent_spatial = function(self, type = NULL, priors = NULL){
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
      setup = function(self, model, ...){
        # Simple security checks
        assertthat::assert_that(
          assertthat::has_name(model, 'background'),
          assertthat::has_name(model, 'biodiversity'),
          nrow(model$predictors) == ncell(self$get_data('template')),
          length(model$biodiversity) == 1 # Only works with single likelihood. To be processed separately
        )
        # In case anything else needs to be specified, do it here
        if(model$biodiversity[[1]]$family=='binomial') assertthat::assert_that(  length( unique(model$biodiversity[[1]]$observations[['observed']])) == 2)
        invisible()
      },
      # Training function
      train = function(self, model, settings, ...){
        assertthat::assert_that(
          inherits(settings,'Settings'),
          is.list(model),length(model)>1,
          # Check that model id and setting id are identical
          settings$modelid == model$id
        )
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
        # Fit the model. Little hack to work correctly with binomials...
        if(is.factor(data$observed)){
          fit_bart <- dbarts::bart(y.train = data[,'observed'],
                                   x.train = data[,model$biodiversity[[1]]$predictors_names],
                                   keeptrees = dc@keepTrees,
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
                                   # weights = w,
                                   ntree = dc@n.trees,
                                   nthread = dc@n.threads,
                                   nchain = dc@n.chains,
                                   nskip = dc@n.burn,
                                   verbose = settings$get('verbose')
          )
        }

        # Predict spatially
        if(!settings$get('inference_only')){
          # Make a prediction
          suppressWarnings(
            pred_bart <- dbarts:::predict.bart(fit_bart,
                                        newdata = full,
                                        type = 'response')
          )
          # Fill output with summaries of the posterior
          prediction[as.numeric(full$cellid)] <- apply(pred_bart, 2, mean)
          names(prediction) <- 'mean'
          prediction <- raster::stack(prediction)

          # Summarize quantiles and sd from posterior
          ms <- as.data.frame(
                 cbind( matrixStats::colQuantiles(pred_bart, probs = c(.05,.5,.95)),
                        matrixStats::colSds(pred_bart)
                       )
          )
          names(ms) <- c('0.05ci','0.5ci','0.95ci','sd')
          # Add them
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
          partial = function(self, x.vars = NULL, ...){

            model <- self$get_data('fit_best')
            assertthat::assert_that(x.vars %in% attr(model$fit$data@x,'term.labels') || NULL,
                                    msg = 'Variable not in predicted model' )
            # Check if family is binomial, if so alter
            not_binomial = self$get_data('model')$biodiversity[[1]]$family != 'binomial'
            partial_effect(model, x.vars = NULL, transform = not_binomial, ... )
          },
          # Spatial partial dependence plot option from embercardo
          spartial = function(self, predictors, x.vars = NULL, equal = FALSE, smooth = 1, transform = TRUE){
            model <- self$get_data('fit_best')
            assertthat::assert_that(x.vars %in% attr(model$fit$data@x,'term.labels'),
                                    msg = 'Variable not in predicted model' )

            if( self$model$biodiversity[[1]]$family != 'binomial') warning('Check whether transform should not be set to False!')

            # Calculate
            p <- partial_space(model, predictors, x.vars, equal, smooth, transform)

            cols <- c("#000004FF","#1B0C42FF","#4B0C6BFF","#781C6DFF","#A52C60FF","#CF4446FF","#ED6925FF","#FB9A06FF","#F7D03CFF","#FCFFA4FF")
            plot(p, col = cols, main = paste0(x.vars, collapse ='|'))
            # Also return spatial
            return(p)
          },
          # Engine-specific projection function
          project = function(self, newdata, rowids){
            assertthat::assert_that(!missing(newdata),
                                    is.data.frame(newdata))

            # Make a prediction
            suppressWarnings(
              pred_bart <- dbarts:::predict.bart(object = self$get_data('fit_best'),
                                                 newdata = newdata,
                                                 type = 'response')
              )
            # Fill output with summaries of the posterior
            prediction <- emptyraster(self$fits$prediction) # Background
            prediction[] <- apply(pred_bart, 2, mean)
            names(prediction) <- 'mean'
            prediction <- raster::stack(prediction)

            # Summarize quantiles and sd from posterior
            ms <- as.data.frame(
              cbind( matrixStats::colQuantiles(pred_bart, probs = c(.05,.5,.95)),
                     matrixStats::colSds(pred_bart)
              )
            )
            names(ms) <- c('0.05ci','0.5ci','0.95ci','sd')
            # Add them
            for(post in names(ms)){
              prediction2 <- self$get_data('template')
              prediction2[as.numeric(full$cellid)] <- ms[[post]]; names(prediction2) <- post
              prediction <- raster::addLayer(prediction, prediction2)
              rm(prediction2)
            }

          }
        )
        return(out)
      }
    )
  ) # End of bdproto object
} # End of function
