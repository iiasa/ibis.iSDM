#' @include bdproto-engine.R utils-spatial.R bdproto-distributionmodel.R
NULL
#' Use of Gradient Descent Boosting for model estimation
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param fam Default is [`Poisson`]
#' @param boosting_iterations An [`integer`] giving the number of boosting iterations
#' @param learning_rate A bounded [`numeric`] value between 0 and 1 defining the shrinkage parameter.
#' @param empirical_risk method for empirical risk calculation ('inbag','oobag','none')
#' @param verbose Should progress be printed?
#' @param ... Other variables or control parameters
#' @references Hofner, B., Mayr, A., Robinzonov, N., & Schmid, M. (2014). Model-based boosting in R: a hands-on tutorial using the R package mboost. Computational statistics, 29(1-2), 3-35.
#' @references Hofner, B., Müller, J., Hothorn, T., 2011. Monotonicity-constrained species distribution models. Ecology 92, 1895–901.
#'
#' @name engine_gdb
NULL
#' @rdname engine_gdb
#' @export
engine_gdb <- function(x,
                       fam = 'Poisson',
                       boosting_iterations = 1000,
                       learning_rate = 0.1,
                       empirical_risk = 'inbag',
                       verbose = FALSE,
                        ...) {
  # Check whether mboost package is available
  check_package('mboost')
  if(!("mboost" %in% loadedNamespaces()) || ('mboost' %notin% sessionInfo()$otherPkgs) ) {
    try({requireNamespace('mboost');attachNamespace("mboost")},silent = TRUE)
    }

  # assert that arguments are valid
  assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                          inherits(x$background,'sf'),
                          is.numeric(boosting_iterations),
                          is.numeric(learning_rate),
                          is.character(empirical_risk),
                          assertthat::is.flag(verbose),
                          empirical_risk %in% c('inbag','oobag','none')
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

  # Set up boosting control
  bc <- mboost::boost_control(mstop = boosting_iterations,
                              nu = learning_rate,
                              risk = empirical_risk,
                              trace = verbose
                              )

  # Print a message in case there is already an engine object
  if(!is.Waiver(x$engine)) message('Replacing currently selected engine.')

  # Detect and format the family
  if( grep('pois',fam,ignore.case = TRUE) ) fam <- mboost::Poisson() else if( grep('bino',fam,ignore.case = TRUE) ) fam <- mboost::Binomial()

  assertthat::assert_that(inherits(fam,'boost_family'),msg = 'Family misspecified.')

  # Set engine in distribution object
  x$set_engine(
    bdproto(
      "GDB-Engine",
      Engine,
      name = "<GDB>",
      data = list(
        'template' = template,
        'bc' = bc,
        'fam' = fam
      ),
      # Function to respecify the control parameters
      set_control = function(self,
                             boosting_iterations = 20,
                             learning_rate = 0.1, # Set relatively low to not regularize too much
                             empirical_risk = 'inbag',
                             verbose = TRUE
                             ){
        # Set up boosting control
        bc <- mboost::boost_control(mstop = boosting_iterations,
                                    nu = learning_rate,
                                    risk = empirical_risk,
                                    trace = verbose
        )
        # Overwrite existing
        self$data$bc <- bc

      },
      # Dummy function for latent factors
      calc_latent_spatial = function(self,...){
        invisible()
      },
      # Get equation for spatial effect
      get_equation_latent_spatial = function(self, spatial_field = c('x','y'),
                                             df = 6, knots = list('x' = 10, 'y' = 10),...){
        return(
          paste0(
            'bspatial(',spatial_field[1],',',spatial_field[2],', center = TRUE, df = ',df,', knots = ',knots,')'
            )
          )
      },
      # Setup function
      setup = function(self, model, ...){
        # Simple security checks
        assertthat::assert_that(
         assertthat::has_name(model, 'background'),
         assertthat::has_name(model, 'data'),
         # Check that all predictors are present
         all( all.vars(model[['equation']][['poipo']]) %in% names(  model$data$poipo_values ) ),
         nrow(model$predictors) == ncell(self$get_data('template')),
         all( sapply(model$equation, is.formula) )
        )
      # Add in case anything needs to be further prepared here

      invisible()
      },
      # Training function
      train = function(self, model, inference_only = FALSE, ...){

        # Get output raster
        prediction <- self$get_data('template')
        # Get boosting control and family data
        bc <- self$get_data('bc')
        fam <- self$get_data('fam')

        full = model$predictors

        if('offset' %in% names(model)){
          # Add offset to full prediction and load vector
          n <- data.frame(model$offset[,3], log(model$offset[,3]) )
          names(n) <- c( names(model$offset)[3], paste0('offset(log(',names(model$offset)[3],'))') )

          full <- cbind(full, n)
          off <- model$data$poipo_values[, names(model$offset)[3] ]
        } else { off = NULL }

        # --- #
        # Fit the base model
        fit_gdb <- mboost::gamboost(
                          formula = model$equation$poipo,
                          data = model$data$poipo_values,
                          weights = model$data$poipo_expect,
                          family = fam,
                          offset = off, # Offset added already
                          control = bc
                          )
        # Variables included
        names(coef(fit_gdb))

        # 5 fold Cross validation to prevent overfitting
        # Andreas Mayr, Benjamin Hofner, and Matthias Schmid (2012). The importance of knowing when to stop - a sequential stopping rule for component-wise gradient boosting. Methods of Information in Medicine, 51, 178–186.
        grs <- seq(from = 10, to = 1000, by = 10)
        cvf <- mboost::cv(model.weights(fit_gdb),B = 5, type = "kfold")
        try({cvm <- mboost::cvrisk(fit_gdb,
                              folds = cvf, grid = grs,
                              papply = pbapply::pblapply )
        },silent = TRUE)

        # Check whether crossvalidation has run through successfully
        if(exists('cvm') && mstop(cvm)>0){
          # Set the model to the optimal mstop to limit overfitting
          fit_gdb[mstop(cvm)]
        }

        # Predict spatially
        if(!inference_only){
          # Make a prediction
          suppressWarnings(
            pred_gdb <- mboost::predict.mboost(object = fit_gdb, newdata = full,
                                               type = 'response', aggregate = 'sum')
          )

          # Fill output
          prediction[] <- pred_gdb[,1]
          names(prediction) <- 'mean'
        } else {
          prediction <- NULL
        }
        # Create output
        out <- bdproto(
          "GDB-Model",
          DistributionModel,
          model = model,
          fits = list(
            "fit_best" = fit_gdb,
            "fit_cv" = cvm,
            "fit_best_equation" = model$equation$poipo,
            "prediction" = prediction
          )
        )
        return(out)


      }
    )
  ) # End of bdproto object
} # End of function
