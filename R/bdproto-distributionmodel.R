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
    # FIXME: Have engine-specific code moved to engine
    if(inherits(self, 'INLA-Model')){
      if( length( self$fits ) != 0 ){
        # Get strongest effects
        ms <- subset(tidy_inla_summary(self$get_data('fit_best')),
                     select = c('variable', 'mean'))
        ms <- ms[order(ms$mean,decreasing = TRUE),] # Sort

        # TODO: See if this can be plotted and extracted and shown here
        # Provides the distance value (in the unit of the point coordinates) above
        # which spatial dependencies become negligible
        # model0.res<-inla.spde2.result(model0, 'spatial.field', spde, do.transf=TRUE)
        # model0.res$summary.log.range.nominal

        message(paste0(
          'Trained ',class(self)[1],' (',self$show(),')',
          '\n  \033[2mStrongest summary effects:\033[22m',
          '\n     \033[34mPositive:\033[39m ', name_atomic(ms$variable[ms$mean>0]),
          '\n     \033[31mNegative:\033[39m ', name_atomic(ms$variable[ms$mean<0])
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
          '\n     ', name_atomic(names(vi))
        ))
    } else if( inherits(self, 'BART-Model') ) {
      # Calculate variable importance from the posterior trees
      vi <- varimp.bart(self$get_data('fit_best'))

      message(paste0(
        'Trained ',class(self)[1],' (',self$show(),')',
        '\n  \033[2mStrongest effects:\033[22m',
        '\n     ', name_atomic(vi$names)
      ))
    } else if( inherits(self, 'STAN-Model') ) {
      # Calculate variable importance from the posterior trees
      summary(self$get_data('fit_best'))

      # message(paste0(
      #   'Trained ',class(self)[1],' (',self$show(),')',
      #   '\n  \033[2mStrongest effects:\033[22m',
      #   '\n     ', name_atomic(vi$names)
      # ))
    } else {
      message(paste0(
        'Trained distribution model (',self$show(),')',
        '\n     No fitted model found!'
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
      assertthat::assert_that(
        inherits(pred,'Raster')
      )
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
    } else if(inherits(self, 'INLA-Model')){
      tidy_inla_summary(self$get_data(x))
    } else if(inherits(self, 'BART-Model')){
      # Number of times each variable is used by a tree split
      # Tends to become less informative with higher numbers of splits
      varimp.bart(self$get_data('fit_best')) %>% tibble::remove_rownames()
    } else if(inherits(self, 'INLA-Model')){
      summary(self$get_data(x))
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
      plot( self$get_data(x) )
    } else if(inherits(self, 'BART-Model')){
      message('Calculating partial dependence plots')
      self$partial(self$get_data(x), x.vars = what, ...)
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
