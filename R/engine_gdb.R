#' @include bdproto-engine.R utils-spatial.R bdproto-distributionmodel.R
NULL
#' Use of Gradient Descent Boosting for model estimation
#'
#' @description
#' Gradient descent boosting is an efficient way to optimize any loss function
#' of a generalized linear or additive model (such as the GAMs available through the [mgcv] R-package).
#' It furthermore automatically regularizes the fit, thus the resulting model only contains the
#' covariates whose baselearners have some influence on the response.
#' Depending on the type of the [add_biodiversity] data, either poisson process models or
#' logistic regressions are estimated. If the \code{only_linear} term in [train] is set to \code{FALSE},
#' splines are added to the estimation, thus providing a non-linear additive inference.
#'
#' @details:
#' This package requires the [mboost] R-package to be installed.
#' It is in philosophy somewhat related to the [engine_xgboost] and [XGBoost] R-package,
#' however providing some additional desirable features that make estimation quicker and
#' particularly useful for spatial projections. Such as for instance the ability to specifically add
#' spatial baselearners via [add_latent] or the specification of monotonically constrained priors
#' via [GDBPrior].
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param boosting_iterations An [`integer`] giving the number of boosting iterations
#' @param learning_rate A bounded [`numeric`] value between \code{0} and \code{1} defining the shrinkage parameter.
#' @param empirical_risk method for empirical risk calculation.
#' Available options are \code{'inbag'}, \code{'oobag'} and \code{'none'}. (Default: \code{'inbag'}).
#' @param type The mode used for creating posterior predictions. Either making \code{"link"}, \code{"response"} or \code{"class"} (Default: \code{"response"}).
#' @param ... Other variables or control parameters
#' @references
#' * Hofner, B., Mayr, A., Robinzonov, N., & Schmid, M. (2014). Model-based boosting in R: a hands-on tutorial using the R package mboost. Computational statistics, 29(1-2), 3-35.
#' * Hofner, B., Müller, J., Hothorn, T., 2011. Monotonicity-constrained species distribution models. Ecology 92, 1895–901.
#'
#' @family engine
#' @name engine_gdb
NULL
#' @rdname engine_gdb
#' @export
engine_gdb <- function(x,
                       boosting_iterations = 1000,
                       learning_rate = 0.1,
                       empirical_risk = 'inbag',
                       type = "response",
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
                          is.character(type),
                          empirical_risk %in% c('inbag','oobag','none')
                          )
  # Match type
  type <- match.arg(type, choices = c("link", "response", "class"), several.ok = FALSE)
  # Get background
  background <- x$background

  # Create a background raster
  if(is.Waiver(x$predictors)){
    # Create from background
    template <- raster::raster(
      ext = raster::extent(background),
      crs = raster::projection(background),
      res = c(diff( (sf::st_bbox(background)[c(1,3)]) ) / 100, # Simplified assumption for resolution
              diff( (sf::st_bbox(background)[c(1,3)]) ) / 100
                    )
                      )
  } else {
    # If predictor existing, use them
    template <- emptyraster(x$predictors$get_data() )
  }

  # Burn in the background
  template <- raster::rasterize(background, template, field = 0)

  # Set up boosting control
  bc <- mboost::boost_control(mstop = boosting_iterations,
                              nu = learning_rate,
                              risk = empirical_risk
                              )

  # Set up the parameter list
  params <- list(
    type = type,
    ...
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
        'bc' = bc,
        'params' = params
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
                                             df = 6, knots = 4,...){
        return(
          paste0(
            'bspatial(',spatial_field[1],',',spatial_field[2],', center = TRUE, df = ',df,', knots = ',knots,')',
            ' + ',
            'bols(',spatial_field[1],')', '+', 'bols(',spatial_field[2],')', '+', 'bols(',spatial_field[1],',',spatial_field[2],')'
            )
          )
      },
      # Setup function
      setup = function(self, model, settings = NULL, ...){
        # Simple security checks
        assertthat::assert_that(
         assertthat::has_name(model, 'background'),
         assertthat::has_name(model, 'biodiversity'),
         inherits(settings,'Settings') || is.null(settings),
         # Check that all predictors are present
         nrow(model$predictors) == raster::ncell(self$get_data('template')),
         length(model$biodiversity) == 1 # Only works with single likelihood. To be processed separately
        )
        # Add in case anything needs to be further prepared here
        # Messenger
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Engine setup.')

        # Add pseudo-absence points if necessary
        # Include nearest predictor values for each
        if('poipo' == model$biodiversity[[1]]$type && model$biodiversity[[1]]$family == 'poisson') {

          # Get background layer
          bg <- self$get_data("template")
          assertthat::assert_that(is.Raster(bg), !is.na(raster::cellStats(bg,min)))

          # # Aggregate observations
          # MJ: Removed for now. For this to work, the non-duplicates need to matched with the envs
          # model$biodiversity[[1]]$observations <- aggregate_observations2grid(
          #   df = model$biodiversity[[1]]$observations,
          #   template = bg,
          #   field_occurrence = "observed")

          # Add pseudo-absence points
          presabs <- add_pseudoabsence(df = model$biodiversity[[1]]$observations,
                                   field_occurrence = 'observed',
                                   template = bg,
                                   settings = model$biodiversity[[1]]$pseudoabsence_settings)
          if(inherits(presabs, 'sf')) presabs <- presabs %>% sf::st_drop_geometry()
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
          if(length(any_missing)>0) presabs <- presabs[-any_missing,] # This works as they are in the same order
          df <- subset(df, complete.cases(df))
          assertthat::assert_that(nrow(presabs) == nrow(df))

          # Check that factors have been correctly set if any
          if(any(model$predictors_types$type=="factor")){
            df[,model$predictors_types$predictors[model$predictors_types$type=="factor"]] <- factor(df[,model$predictors_types$predictors[model$predictors_types$type=="factor"]])
          }

          # Overwrite observation data
          model$biodiversity[[1]]$observations <- presabs

          # Preprocessing security checks
          assertthat::assert_that( all( model$biodiversity[[1]]$observations[['observed']] >= 0 ),
                                   any(!is.na(presabs[['observed']])),
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
            #                     longlat = raster::isLonLat(bg),
            #                     field_space = c('x','y')
            #                     )
            model$biodiversity[[1]]$offset <- ofs
          }
          # Define expectation as very small vector following Renner et al.
          w <- ppm_weights(df = df,
                           pa = model$biodiversity[[1]]$observations[['observed']],
                           bg = bg,
                           weight = 1e-6
          )
          df$w <- w # Also add as column

          model$biodiversity[[1]]$predictors <- df
          model$biodiversity[[1]]$expect <- w
        } else if(model$biodiversity[[1]]$family != 'poisson'){
          # calculating the case weights (equal weights)
          # the order of weights should be the same as presences and backgrounds in the training data
          prNum <- as.numeric(table(model$biodiversity[[1]]$observations[['observed']])["1"]) # number of presences
          bgNum <- as.numeric(table(model$biodiversity[[1]]$observations[['observed']])["0"]) # number of backgrounds
          w <- ifelse(model$biodiversity[[1]]$observations[['observed']] == 1, 1, prNum / bgNum)
          # Weights for IWLR
          # w <- (10^6)^(1 - model$biodiversity[[1]]$observations[['observed']])

          # If family is not poisson, assume factor distribution for response
          model$biodiversity[[1]]$observations[['observed']] <- factor(model$biodiversity[[1]]$observations[['observed']])

          model$biodiversity[[1]]$expect <- w
        }

        # ---- #
        # Detect and format the family
        fam <- model$biodiversity[[1]]$family
        fam <- switch (fam,
          "poisson" = mboost::Poisson(),
          "binomial" = mboost::Binomial(type = "glm", link = "cloglog"),
          "gaussian" = Gaussian(),
          "hurdle" = Hurdle(nuirange = c(0,100))
        )
        self$data[['family']] <- fam
        assertthat::assert_that(inherits(fam,'boost_family'),msg = 'Family misspecified.')

        # Instead of invisible return the model object
        return( model )
      },
      # Training function
      train = function(self, model, settings, ...){
        # Messager
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Starting fitting...')

        # Get output raster
        prediction <- self$get_data('template')
        # Get boosting control and family data
        bc <- self$get_data('bc')
        bc$trace <- settings$get('verbose') # Reset trace as specified
        params <- self$get_data('params')
        fam <- self$get_data('family')

        # All other needed data for model fitting
        equation <- model$biodiversity[[1]]$equation
        data <- cbind(model$biodiversity[[1]]$predictors, data.frame(observed = model$biodiversity[[1]]$observations[,'observed']) )
        w <- model$biodiversity[[1]]$expect
        full <- model$predictors

        # Select predictors
        full <- subset(full, select = c('x','y',model$biodiversity[[1]]$predictors_names))
        full$cellid <- rownames(full) # Add row.names
        full <- subset(full, complete.cases(full))

        assertthat::assert_that(
          is.null(w) || length(w) == nrow(data),
          is.formula(equation),
          all( model$biodiversity[[1]]$predictors_names %in% names(full) )
        )

        if(!is.Waiver(model$offset)){
          # Add offset to full prediction and load vector
          n <- data.frame(model$offset[as.numeric(full$cellid), "spatial_offset"], model$offset[as.numeric(full$cellid), "spatial_offset"] )
          names(n) <- c( "spatial_offset", paste0('offset(',"spatial_offset",')') )
          full <- cbind(full, n)
          # And for biodiversity object
          n <- cbind(model$biodiversity[[1]]$offset[,"spatial_offset"],
                     model$biodiversity[[1]]$offset[,"spatial_offset"]) |> as.data.frame()
          names(n) <- c( "spatial_offset", paste0('offset(',"spatial_offset",')') )
          data <- cbind(data, n)
        }

        # --- #
        # Fit the base model
        fit_gdb <- try({
          mboost::gamboost(
            formula = equation,
            data = data,
            # weights = w,
            family = fam,
            offset = w, # Add exposure as offset
            control = bc
          )
        },silent = FALSE)
        # If error, decrease step size by a factor of 10 and try again.
        if(inherits(fit_gdb, "try-error") || length(names(coef(fit_gdb)))< 2){
          if(getOption('ibis.setupmessages')) myLog('[Estimation]','red','Reducing learning rate by 1/100.')
          bc$nu <- bc$nu * 0.01
          fit_gdb <- try({
            mboost::gamboost(
              formula = equation,
              data = data,
              # weights = w,
              family = fam,
              offset = w, # Add exposure as offset
              control = bc
            )
          },silent = FALSE)
          if(inherits(fit_gdb, "try-error")) myLog('[Estimation]','red','Fitting failed. Check model and alter parameters!')
        }

        if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Starting cross.validation.')
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
            pred_gdb <- mboost::predict.mboost(object = fit_gdb, newdata = full,
                                               type = self$get_data('params')$type, aggregate = 'sum')
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
            "params" = params,
            "fit_best_equation" = equation,
            "prediction" = prediction
          ),
          # Project function
          project = function(self, newdata, type = NULL){
            assertthat::assert_that('fit_best' %in% names(self$fits),
                                    is.data.frame(newdata) || is.matrix(newdata),
                                    assertthat::has_name(newdata,c('x','y'))
                                    )
            # Get model
            mod <- self$get_data('fit_best')
            assertthat::assert_that(inherits(mod,'mboost'),msg = 'No model found!')
            if(is.null(type)) type <- self$get_data('params')$type
            # Check that all variables are in provided data.frame
            assertthat::assert_that(all( as.character(mboost::extract(mod,'variable.names')) %in% names(newdata) ))

            # Add rowid
            newdata$rowid <- 1:nrow(newdata)
            # Subset to non-missing data
            newdata <- subset(newdata, complete.cases(newdata)) # FIXME: Potentially limit to relevant variables only
            # Make empty template
            temp <- raster::rasterFromXYZ(newdata[,c('x','y')])
            # Predict
            y <- suppressWarnings(
              mboost::predict.mboost(object = mod,newdata = newdata,
                                     type = type, aggregate = 'sum')
            )
            temp[as.numeric(newdata$rowid)] <- y[,1]
            names(temp) <- "mean" # Rename to mean
            return(temp)
          },
          # Partial effect
          partial = function(self, x.var, constant = NULL, variable_length = 100, plot = FALSE){
            # Assert that variable(s) are in fitted model
            assertthat::assert_that( is.character(x.var),inherits(self$get_data('fit_best'), 'mboost') )
            # Unlike the effects function, build specific predictor for target variable(s) only
            variables <- mboost::extract(self$get_data('fit_best'),'variable.names')
            assertthat::assert_that( all( x.var %in% variables), msg = 'x.var variable not found in model!' )
            # Special treatment for factors
            if(any(model$predictors_types$type=="factor")){
              if(x.var %in% model$predictors_types$predictors[model$predictors_types$type=="factor"]){
                variable_range <- levels(self$model$predictors[,x.var])
              } else {
                variable_range <- range(self$model$predictors[[x.var]],na.rm = TRUE)
              }
            } else {
              variable_range <- range(self$model$predictors[[x.var]],na.rm = TRUE)
            }

            # Create dummy data.frame
            if(is.character(variable_range)){
              # For factors, just add them
              dummy <- as.data.frame(matrix(nrow = length(variable_range)))
              dummy[,x.var] <- factor(variable_range)
            } else {
              dummy <- as.data.frame(matrix(nrow = variable_length))
              dummy[,x.var] <- seq(variable_range[1],variable_range[2],length.out = variable_length)
            }
            # For the others
            if(is.null(constant)){
              if(any(self$model$predictors_types$type=='factor')){
                # Numeric names
                nn <- self$model$predictors_types$predictors[which(self$model$predictors_types$type=='numeric')]
                constant <- apply(self$model$predictors[,nn], 2, function(x) mean(x, na.rm=T))
                dummy <- cbind(dummy,t(constant))
                # For each factor duplicate the entire matrix and add factor levels
                # nf <- self$model$predictors_types$predictors[which(self$model$predictors_types$type=='factor')]
                # for(n in nf){
                #   var <- data.frame( factor( rep(levels(self$model$predictors[,nf]), nrow(dummy)) ) )
                #   names(var) <- n
                #   dummy <- cbind(dummy,var);rm(var)
                # }
              } else {
                # Calculate mean
                constant <- apply(self$model$predictors, 2, function(x) mean(x, na.rm=T))
                dummy <- cbind(dummy,t(constant))
              }
            } else {
              dummy[,variables] <- constant
            }

            # Now predict with model
            pp <- mboost::predict.mboost(object = self$get_data('fit_best'), newdata = dummy,
                                         which = x.var,
                                         type = self$get_data('params')$type, aggregate = 'sum')
            # Combine with
            out <- data.frame(partial = pp[,grep(x.var, colnames(pp))] ); out[[x.var]] <- dummy[[x.var]]

            # If plot, make plot, otherwise
            if(plot){
              mboost::plot.mboost(self$get_data('fit_best'), which = x.var, newdata = dummy)
            }
            return(out)
          },
          # Spatial partial effect plots
          spartial = function(self, x.var, constant = NULL, plot = TRUE,...){
            assertthat::assert_that('fit_best' %in% names(self$fits),
                                    is.character(x.var), length(x.var) == 1)
            # Get model and make empty template
            mod <- self$get_data('fit_best')
            model <- self$model
            # Also check that what is present in coefficients of model
            variables <- as.character( mboost::extract(mod,'variable.names') )
            assertthat::assert_that(x.var %in% variables )

            # Make template of target variable(s)
            temp <- raster::rasterFromXYZ(cbind(model$predictors$x,model$predictors$y),
                                          crs = raster::projection(model$background))
            # Get target variables and predict
            target <- model$predictors
            # Set all variables other the target variable to constant
            if(is.null(constant)){
              # Calculate mean
              nn <- self$model$predictors_types$predictors[which(self$model$predictors_types$type=='numeric')]
              constant <- apply(target[,nn], 2, function(x) mean(x, na.rm=T))
              for(v in variables[ variables %notin% x.var]){
                if(v %notin% names(target) ) next()
                target[!is.na(target[v]),v] <- constant[v]
              }
            } else {
              target[!is.na(target[,x.var]), variables] <- constant
            }
            target$rowid <- as.numeric( rownames(target) )
            assertthat::assert_that(nrow(target)==ncell(temp))

            pp <- suppressWarnings(
              mboost::predict.mboost(mod, newdata = target, which = x.var,
                                     type = 'link', aggregate = 'sum')
            )
            # If both linear and smooth effects are in model
            if(length(target$rowid[which(!is.na(target[[x.var]]))] ) == length(pp[,ncol(pp)])){
              temp[ target$rowid[which(!is.na(target[[x.var]]))] ] <- pp[,ncol(pp)]
            } else { temp[] <- pp[, ncol(pp) ]}
            names(temp) <- paste0('partial__',x.var)

            if(plot){
              # Plot both partial spatial partial
              par.ori <- par(no.readonly = TRUE)
              par(mfrow = c(1,3))
              raster::plot(temp, main = expression(f[partial]), col = ibis_colours$divg_bluegreen)
              mboost::plot.mboost(mod,which = x.var)
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
              raster::plot(temp, main = expression(f[partial]), col = ibis_colours$divg_bluegreen )
            }

            return(temp)

          }
        )
        return(out)
      }
    )
  ) # End of bdproto object
} # End of function
