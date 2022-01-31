#' @include utils.R bdproto.R identifier.R
NULL

#' @export
if (!methods::isClass("DistributionModel")) methods::setOldClass("DistributionModel")
NULL

#' Prototype for the trained Model object
#'
#' All trained Models should inherit the options here
#'
#' @name DistributionModel-class
#' @aliases DistributionModel
#' @family bdproto
#' @keywords bdproto
NULL

#' @export
DistributionModel <- bdproto(
  "DistributionModel",
  id = character(), # An id for any trained model
  model = list(),
  fits = list(), # List of fits with data
  # Print message with summary of model
  print = function(self) {
    # TODO: Have a lot more information in here and to be prettified

    # Check whether prediction exists and number of layers
    has_prediction <- "prediction" %in% self$show_rasters()
    # Check whether threshold has been calculated
    has_threshold <- grep('threshold',self$show_rasters(),value = TRUE)[1]

    # FIXME: Have engine-specific code moved to engine
    if(inherits(self, 'INLA-Model') || inherits(self, 'INLABRU-Model') ){
      if( length( self$fits ) != 0 ){
        # Get strongest effects
        ms <- subset(tidy_inla_summary(self$get_data('fit_best')),
                     select = c('variable', 'mean'))
        ms <- ms[order(ms$mean,decreasing = TRUE),] # Sort

        message(paste0(
          'Trained ',class(self)[1],' (',self$show(),')',
          '\n  \033[2mStrongest summary effects:\033[22m',
          '\n     \033[34mPositive:\033[39m ', name_atomic(ms$variable[ms$mean>0]),
          '\n     \033[31mNegative:\033[39m ', name_atomic(ms$variable[ms$mean<0]),
          ifelse(has_prediction,
                 paste0("\n  Prediction fitted: ",text_green("yes")),
                 ""),
          ifelse(!is.na(has_threshold),
                 paste0("\n  Threshold created: ",text_green("yes")),
                 "")
        ))
      }
    } else if( inherits(self, 'GDB-Model') ) {

        # Get Variable importance
        vi <- mboost::varimp(
          self$get_data('fit_best')
        )
        vi <- sort( vi[which(vi>0)],decreasing = TRUE )

        message(paste0(
          'Trained ',class(self)[1],' (',self$show(),')',
          '\n  \033[2mStrongest effects:\033[22m',
          '\n     ', name_atomic(names(vi)),
          ifelse(has_prediction,
                 paste0("\n  Prediction fitted: ",text_green("yes")),
                 ""),
          ifelse(!is.na(has_threshold),
                 paste0("\n  Threshold created: ",text_green("yes")),
                 "")
        ))
    } else if( inherits(self, 'BART-Model') ) {
      # Calculate variable importance from the posterior trees
      vi <- varimp.bart(self$get_data('fit_best'))

      message(paste0(
        'Trained ',class(self)[1],' (',self$show(),')',
        '\n  \033[2mStrongest effects:\033[22m',
        '\n     ', name_atomic(vi$names),
        ifelse(has_prediction,
               paste0("\n  Prediction fitted: ",text_green("yes")),
               ""),
        ifelse(!is.na(has_threshold),
               paste0("\n  Threshold created: ",text_green("yes")),
               "")
      ))
    } else if( inherits(self, 'STAN-Model') ) {
      # Calculate variable importance from the posterior trees
      summary(self$get_data('fit_best'))$summary |> as.data.frame() |>
        tibble::rownames_to_column(var = "parameter") |> tibble::as_tibble(rownames = NULL)

      # message(paste0(
      #   'Trained ',class(self)[1],' (',self$show(),')',
      #   '\n  \033[2mStrongest effects:\033[22m',
      #   '\n     ', name_atomic(vi$names)
      # ))
    } else if( inherits(self, 'XGBOOST-Model') ) {
      vi <- xgboost::xgb.importance(model = self$get_data('fit_best'))

      message(paste0(
        'Trained ',class(self)[1],' (',self$show(),')',
        '\n  \033[2mStrongest effects:\033[22m',
        '\n     ', name_atomic(vi$Feature),
        ifelse(has_prediction,
               paste0("\n  Prediction fitted: ",text_green("yes")),
               ""),
        ifelse(!is.na(has_threshold),
               paste0("\n  Threshold created: ",text_green("yes")),
               "")
      ))
    } else if( inherits(self, 'BREG-Model') ) {
      obj <- self$get_data('fit_best')
      # Summarize the beta coefficients from the posterior
      ms <- posterior::summarise_draws(obj$beta) |>
        subset(select = c('variable', 'mean'))
      # Reorder
      ms <- ms[order(ms$mean,decreasing = TRUE),] # Sort

      message(paste0(
        'Trained ',class(self)[1],' (',self$show(),')',
        '\n  \033[2mStrongest summary effects:\033[22m',
        '\n     \033[34mPositive:\033[39m ', name_atomic(ms$variable[ms$mean>0]),
        '\n     \033[31mNegative:\033[39m ', name_atomic(ms$variable[ms$mean<0]),
        ifelse(has_prediction,
               paste0("\n  Prediction fitted: ",text_green("yes")),
               ""),
        ifelse(!is.na(has_threshold),
               paste0("\n  Threshold created: ",text_green("yes")),
               "")
      ))


    } else {
      message(paste0(
        'Trained distribution model (',self$show(),')',
        text_red('\n     No fitted model found!'),
        ifelse(has_prediction,
               paste0("\n  Prediction fitted: ",text_green("yes")),
               ""),
        ifelse(!is.na(has_threshold),
               paste0("\n  Threshold created: ",text_green("yes")),
               "")
      ))
    }
  },
  # Show the name of the Model
  show = function(self) {
    self$model$runname
  },
  # Plot the prediction
  plot = function(self, what = 'mean',...){
    if( length( self$fits ) != 0 && !is.null( self$fits$prediction ) ){
      pred <- self$get_data('prediction')
      assertthat::assert_that(is.Raster(pred))
      # Match arguement
      what <- match.arg(what, names(pred), several.ok = FALSE)
      assertthat::assert_that( what %in% names(pred),msg = paste0('Prediction type not found. Available: ', paste0(names(pred),collapse = '|')))
      raster::plot(pred[[what]],
           main = paste0(self$model$runname, ' prediction (',what,')'),
           box = FALSE,
           axes = TRUE,
           colNA = NA, col = ibis_colours[['sdm_colour']],...)

    } else {
      message(
        paste0('No model predictions found.')
      )
    }
  },
  # Show model run time if settings exist
  show_duration = function(self){
    if(!is.Waiver(self$settings)) self$settings$duration()
  },
  # Get effects or importance tables from model
  summary = function(self, x = 'fit_best'){
    # Distinguishing between model types
    if(inherits(self, 'GDB-Model')){
      mboost:::summary.mboost(self$get_data(x))
    } else if(inherits(self, 'INLA-Model') || inherits(self, 'INLABRU-Model')){
      tidy_inla_summary(self$get_data(x))
    } else if(inherits(self, 'BART-Model')){
      # Number of times each variable is used by a tree split
      # Tends to become less informative with higher numbers of splits
      varimp.bart(self$get_data(x)) %>% tibble::remove_rownames()
    } else if(inherits(self, 'STAN-Model')){
      summary(self$get_data(x))
    } else if(inherits(self, 'BREG-Model')){
      posterior::summarise_draws(self$get_data(x)$beta)
    } else if(inherits(self, "XGBOOST-Model")){
      xgboost::xgb.importance(model = self$get_data('fit_best'))
      # xgboost::xgb.ggplot.importance(o)
    }
  },
  # Dummy partial response calculation. To be overwritten per engine
  partial = function(self){
    new_waiver()
  },
  # Dummy spartial response calculation. To be overwritten per engine
  spartial = function(self){
    new_waiver()
  },
  # Generic plotting function for effect plots
  effects = function(self, x = 'fit_best', what = 'fixed', ...){
    assertthat::assert_that(is.character(what))
    if(inherits(self, 'GDB-Model')){
      # How many effects
      n <- length( coef( self$get_data(x) ))
      # Use the base plotting
      par.ori <- par(no.readonly = TRUE)
      par(mfrow = c(ceiling(n/3),3))

      mboost:::plot.mboost(x = self$get_data(x),
                           type = 'b',cex.axis=1.5, cex.lab=1.5)

      par(par.ori)#dev.off()
    } else if(inherits(self, 'INLA-Model')) {
      plot_inla_marginals(self$get_data(x),what = 'fixed')
    } else if(inherits(self, 'STAN-Model')) {
      # Get true beta parameters
      ra <- grep("beta", names(self$get_data(x)),value = TRUE) # Get range
      rstan::stan_plot(self$fits$fit_best,pars = ra)
    } else if(inherits(self, 'INLABRU-Model')) {
      # Use inlabru effect plot
      ggplot2::ggplot() +
        inlabru:::gg(self$get_data(x)$summary.fixed, bar = TRUE)
    } else if(inherits(self, 'BART-Model')){
      message('Calculating partial dependence plots')
      self$partial(self$get_data(x), x.vars = what, ...)
    } else if(inherits(self, 'BREG-Model')){
      obj <- self$get_data(x)
      if(what == "fixed") what <- "coefficients"
      what <- match.arg(what, choices = c("coefficients", "scaled.coefficients","residuals",
                                           "size", "fit", "help", "inclusion"), several.ok = FALSE)
      if( length( grep("poisson", obj$call) ) > 0 ){
        BoomSpikeSlab::plot.poisson.spike(obj, y = what)
      } else if( length( grep("binomial", obj$call) ) > 0 ){
        BoomSpikeSlab::plot.logit.spike(obj, y = what)
      } else {
        BoomSpikeSlab::plot.lm.spike(obj, y = what)
      }
    } else {
      self$partial(self$get_data(x), x.vars = NULL)
    }
  },
  # Get specific fit from this Model
  get_data = function(self, x) {
    if (!x %in% names(self$fits))
      return(new_waiver())
    return(self$fits[[x]])
  },
  # Set fit for this Model
  set_data = function(self, x, value) {
    # Get biodiversity dataset collection
    ff <- self$fits
    # Set the object
    ff[[x]] <- value
    bdproto(NULL, self, fits = ff )
  },
  # List all rasters in object
  show_rasters = function(self){
    rn <- names(self$fits)
    rn <- rn[ which( sapply(rn, function(x) is.Raster(self$get_data(x)) ) ) ]
    return(rn)
  },
  # Save object
  save = function(self, fname, type = 'gtif', dt = 'FLT4S'){
    assertthat::assert_that(
      is.character(fname),
      type %in% c('gtif','gtiff','tif','nc','ncdf'),
      'fits' %in% self$ls(),
      dt %in% c('LOG1S','INT1S','INT1U','INT2S','INT2U','INT4S','INT4U','FLT4S','FLT8S')
    )
    type <- tolower(type)

    # Get raster file in fitted object
    cl <- sapply(self$fits, class)
    ras <- self$fits[[grep('raster', cl,ignore.case = T)]]

    # Check that no-data value is not present in ras
    assertthat::assert_that(any(!cellStats(ras,min) <= -9999),msg = 'No data value -9999 is potentially in prediction!')

    if(file.exists(fname)) warning('Overwritting existing file...')
    if(type %in% c('gtif','gtiff','tif')){
      # Save as geotiff
      writeGeoTiff(ras, fname = fname, dt = dt)
    } else if(type %in% c('nc','ncdf')) {
      # Save as netcdf
      # TODO: Potentially change the unit descriptions
      writeNetCDF(ras, fname = fname, varName = 'iSDM prediction', varUnit = "",varLong = "")
    }
    invisible()
  }
)
