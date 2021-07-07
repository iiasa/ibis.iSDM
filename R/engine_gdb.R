#' @include bdproto-engine.R utils-spatial.R bdproto-distributionmodel.R
NULL
#' Use of Gradient Descent Boosting for model estimation
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param boosting_iterations An [`integer`] giving the number of boosting iterations
#' @param learning_rate A bounded [`numeric`] value between 0 and 1 defining the shrinkage parameter.
#' @param empirical_risk method for empirical risk calculation ('inbag','oobag','none')
#' @param ... Other variables or control parameters
#' @references Hofner, B., Mayr, A., Robinzonov, N., & Schmid, M. (2014). Model-based boosting in R: a hands-on tutorial using the R package mboost. Computational statistics, 29(1-2), 3-35.
#' @references Hofner, B., Müller, J., Hothorn, T., 2011. Monotonicity-constrained species distribution models. Ecology 92, 1895–901.
#'
#' @name engine_gdb
NULL
#' @rdname engine_gdb
#' @export
engine_gdb <- function(x,
                       boosting_iterations = 1000,
                       learning_rate = 0.1,
                       empirical_risk = 'inbag',
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
                              risk = empirical_risk
                              )

  # Print a message in case there is already an engine object
  if(!is.Waiver(x$engine)) message('Replacing currently selected engine.')

  # Set engine in distribution object
  x$set_engine(
    bdproto(
      "GDB-Engine",
      Engine,
      name = "<GDB>",
      data = list(
        'template' = template,
        'bc' = bc
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
                                             df = 6, knots = 10,...){
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
         assertthat::has_name(model, 'biodiversity'),
         # Check that all predictors are present
         all(all.vars(model$biodiversity[[1]]$equation)[-1] %in% names(model$biodiversity[[1]]$predictors)),
         nrow(model$predictors) == ncell(self$get_data('template')),
         length(model$biodiversity) == 1 # Only works with single likelihood. To be processed separately
        )
        # Add in case anything needs to be further prepared here

        # Detect and format the family
        fam <- model$biodiversity[[1]]$family
        if( fam == 'poisson' ) {fam <- mboost::Poisson()} else if( fam == 'binomial') fam <- mboost::Binomial(type = "glm", link = "cloglog")
        self$data[['family']] <- fam
        assertthat::assert_that(inherits(fam,'boost_family'),msg = 'Family misspecified.')

        invisible()
      },
      # Training function
      train = function(self, model, settings, ...){
        # Get output raster
        prediction <- self$get_data('template')
        # Get boosting control and family data
        bc <- self$get_data('bc')
        bc$trace <- settings$get('verbose') # Reset trace as specified
        fam <- self$get_data('family')

        # All other needed data for model fitting
        equation <- model$biodiversity[[1]]$equation
        data <- cbind(model$biodiversity[[1]]$predictors, data.frame(observed = model$biodiversity[[1]]$observations[,'observed']) )
        if(model$biodiversity[[1]]$family=='binomial') data$observed <- factor(data$observed)
        w <- model$biodiversity[[1]]$expect
        full <- model$predictors

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
          n <- data.frame(model$offset[as.numeric(full$cellid),3], log(model$offset[as.numeric(full$cellid),3]) )
          names(n) <- c( names(model$offset)[3], paste0('offset(log(',names(model$offset)[3],'))') )

          full <- cbind(full, n)
          off <- model$biodiversity[[1]]$offset[, names(model$offset)[3] ]
        } else { off = NULL }

        # --- #
        # Fit the base model
        fit_gdb <- mboost::gamboost(
                          formula = equation,
                          data = data,
                          weights = w,
                          family = fam,
                          # offset = off, # Offset added already
                          control = bc
                          )
        # Variables included
        names(coef(fit_gdb))

        # 5 fold Cross validation to prevent overfitting
        # Andreas Mayr, Benjamin Hofner, and Matthias Schmid (2012). The importance of knowing when to stop - a sequential stopping rule for component-wise gradient boosting. Methods of Information in Medicine, 51, 178–186.
        grs <- seq(from = 10, to = max( bc$mstop *5), by = 10)
        cvf <- mboost::cv(model.weights(fit_gdb),B = 5, type = "kfold")
        try({cvm <- mboost::cvrisk(fit_gdb,
                              folds = cvf, grid = grs,
                              papply = pbapply::pblapply )
        },silent = TRUE)
        # TODO: parallize better? parallel::mclapply(mc.cores = getOption('ibis.nthread'))

        # Check whether crossvalidation has run through successfully
        if(exists('cvm') && mstop(cvm) > 0){
          # Set the model to the optimal mstop to limit overfitting
          fit_gdb[mstop(cvm)]
        } else {
          cvm <- new_waiver()
        }

        # Predict spatially
        if(!settings$get('inference_only')){
          # Make a prediction
          suppressWarnings(
            pred_gdb <- mboost::predict.mboost(object = fit_gdb, newdata = full,
                                               type = 'response', aggregate = 'sum')
          )
          # Fill output
          prediction[as.numeric(full$cellid)] <- pred_gdb[,1]
          names(prediction) <- 'mean'

        } else {
          prediction <- NULL
        }

        # Compute end of computation time
        settings$set('end.time', Sys.time())
        # Also append boosting control option to settings
        for(entry in names(bc)) settings$set(entry, bc[[entry]])

        # Create output
        out <- bdproto(
          "GDB-Model",
          DistributionModel,
          id = new_id(),
          model = model,
          settings = settings,
          fits = list(
            "fit_best" = fit_gdb,
            "fit_cv" = cvm,
            "fit_best_equation" = equation,
            "prediction" = prediction
          ),
          # Project function
          project = function(self, newdata){
            assertthat::assert_that('fit_best' %in% names(self$fits),
                                    is.data.frame(newdata) || is.matrix(newdata),
                                    assertthat::has_name(newdata,c('x','y'))
                                    )
            # Get model
            mod <- self$get_data('fit_best')
            # Check that all variables are in provided data.frame
            assertthat::assert_that(all( as.character(mboost::extract(mod,'variable.names')) %in% names(newdata) ))

            # Make empty template
            temp <- raster::rasterFromXYZ(newdata[,c('x','y')])
            # Predict
            y <- suppressWarnings(
              mboost::predict.mboost(object = mod,newdata = newdata,
                                     type = 'response', aggregate = 'sum')
            )
            assertthat::assert_that(raster::ncell(temp)==nrow(y))
            temp[] <- y[,1]
            return(temp)
          },
          # Spatial partial effect plots
          spartial = function(self, what, plot = TRUE){
            assertthat::assert_that('fit_best' %in% names(self$fits),
                                    is.character(what))
            # Get model and make empty template
            mod <- self$get_data('fit_best')
            # Also check that what is present in coefficients of model
            assertthat::assert_that(what %in%  as.character( mboost::extract(mod,'variable.names') ))

            # Make template of target variable(s)
            temp <- raster::rasterFromXYZ(cbind(self$model$predictors$x,self$model$predictors$y))
            # Get target variables and predict
            target <- self$model$predictors[,c('x','y',what)]
            target$rowid <- as.numeric( rownames(target) )
            assertthat::assert_that(nrow(target)==ncell(temp))

            y <- suppressWarnings(
              mboost::predict.mboost(mod,newdata = target, which = what)
            )
            temp[target$rowid[which(!is.na(target[[what]]))]] <- y[,ncol(y)]
            names(temp) <- paste0('partial__',what)

            if(plot){
              # Plot both partial spatial partial
              par.ori <- par(no.readonly = TRUE)
              par(mfrow = c(1,3))
              raster::plot(temp,main = expression(f[partial]), col =
                             c("#2C194C","#284577","#4B76A0","#8CA7C3","#D0DCE6","#D4E6D6","#98C39B","#5C9F61","#3E7229","#434C01")
              )
              mboost::plot.mboost(mod,which = what)
              par(par.ori)
            }

            return(temp)
          },
          # Spatial latent effect
          plot_spatial = function(self, plot = TRUE){
            assertthat::assert_that('fit_best' %in% names(self$fits) )
            # Get model and make empty template
            mod <- self$get_data('fit_best')
            # Also check that what is present in coefficients of model
            vars <- as.character( mboost::extract(mod,'bnames') )
            assertthat::assert_that(length(grep('bspatial',vars))>0,
                                    msg = 'No spatial effect found in model!')

            # Make template of target variable(s)
            temp <- raster::rasterFromXYZ(cbind(self$model$predictors$x,self$model$predictors$y))
            # Get target variables and predict
            target <- self$model$predictors[,c('x','y')]

            y <- suppressWarnings(
              mboost::predict.mboost(mod,newdata = target, which = c('x','y'))
            )
            assertthat::assert_that(nrow(target)==nrow(y))
            temp[] <- y[,2]
            names(temp) <- paste0('partial__','space')
            # Mask with background
            temp <- raster::mask(temp, self$model$background )

            if(plot){
              # Plot both partial spatial partial
              raster::plot(temp, main = expression(f[partial]), col =
                             c("#2C194C","#284577","#4B76A0","#8CA7C3","#D0DCE6","#D4E6D6","#98C39B","#5C9F61","#3E7229","#434C01")
              )
            }

            return(temp)

          }
        )
        return(out)
      }
    )
  ) # End of bdproto object
} # End of function
