#' @include utils.R bdproto.R
NULL

#' @export
if (!methods::isClass("BiodiversityScenario")) methods::setOldClass("BiodiversityScenario")
NULL

#' Prototype for a biodiversity scenario from a trained model
#'
#' @name BiodiversityScenario-class
#'
#' @aliases BiodiversityScenario
NULL

#' @export
BiodiversityScenario <- bdproto(
  "BiodiversityScenario",
  modelobject = new_waiver(), # The id of the model
  modelid = new_waiver(),
  predictors = new_waiver(),
  constraints = new_waiver(),
  scenarios = new_waiver(),
  # Print message with summary of model
  print = function(self) {
    # Check that model exists
    fit <- self$get_model()
    timeperiod <- self$get_timeperiod()
    # Get set predictors and time period
    pn = ifelse(is.Waiver(self$get_predictor_names()),'None',name_atomic(self$get_predictor_names(), "predictors"))
    tp = ifelse(is.Waiver(timeperiod),'None',
                paste0(
                  paste0( timeperiod,collapse = ' -- '),
                  ' (',round(as.numeric(difftime(self$get_timeperiod()[2],self$get_timeperiod()[1],unit="weeks"))/52.25,1),' years)'
                )
    )
    # Constrains
    cs = self$get_constrains()

    message(paste0('Spatial-temporal scenario:',
                   '\n  Used model: ',ifelse(is.Waiver(fit), text_red('None'), class(fit)[1] ),
                   "\n --------- ",
                   "\n  Predictors:     ", pn,
                   "\n  Time period:    ", tp,
                   "\n --------- ",
                   "\n  Scenarios fitted: ", ifelse(is.Waiver(self$scenarios),text_yellow('None'),'Yes')
      )
    )
  },
  # Verify that set Model exist
  verify = function(self){
    assertthat::validate_that( !is.Waiver(self$modelobject),
                               !is.Waiver(self$modelid),
                               exists(self$modelobject) )
    # Get Model object and check that ID is correct
    x <- get(self$modelobject)
    assertthat::validate_that(x$id == self$modelid)
    invisible()
  },
  # Generic projection function | Model type specific (based on model type)
  # Intended to do the basic data preparation and aggregations of predictions
  project = function(self){
    # Get Model object
    fit <- self$get_model()
    # Get predictors
    new_preds <- self$get_predictors()
    pred_names <- self$get_predictor_names()

    # Now convert to data.frame
    new_preds <- new_preds$get_data(df = TRUE)
    # convert time dimension to Posix
    new_preds$Time <- as.POSIXct( new_preds$Time )
    # Convert all units classes to numeric to avoid problems
    new_preds[,pred_names] <- apply(new_preds[,pred_names], 2, function(x) as.numeric(x) )

    # --- #
    # Now for each unique element, loop and project in order
    for(times in sort(unique(new_preds$Time))){
      myLog('Projecting time stamp', times)
      nd <- subset(new_preds, Time == times)
      out <- fit$project(newdata = nd)
    }

    # Check whether projection function exists in model
    if('project' %notin% names(fit) ){
      warning('Model engine has no projection method yet implemented!')
      return(NULL)
    }

    if(inherits(fit, 'GDB-Model')) {
    } else if(inherits(fit, 'INLA-Model')) {
    } else if(inherits(fit, 'BART-Model')) {

    }

    if(!is.Waiver(fit) && !is.Waiver(new_preds)){
      proj <- predict(fit, newdata = new_preds)
      return(proj)
    } else { NULL }
  },
  # Get Model
  get_model = function(self){
    if(is.Waiver(self$modelobject)) return( new_waiver() )
      else return( get(self$modelobject) )
  },
  # Get Model predictors
  get_predictor_names = function(self) {
    if(is.Waiver(self$predictors)) return(self$predictors)
    if(inherits(self$predictors, "PredictorDataset")) {
      self$predictors$get_names()
    } else {
      stop("Feature data is of an unrecognized class")
    }
  },
  # Get time period of projection
  get_timeperiod = function(self){
    if(is.Waiver(self$predictors)) return(self$predictors)
    if(inherits(self$predictors, "PredictorDataset")) {
      return(
        c( min(self$predictors$timeperiod), max(self$predictors$timeperiod) )
      )
    }
  },
  # Get constrains for model
  get_constrains = function(self){
    # TODO:
    return(new_waiver())
  },
  # Show the name of the Model
  show = function(self) {
    self$modelobject
  },
  # Set Predictors
  set_predictors = function(self, x){
    assertthat::assert_that(inherits(x, "PredictorDataset"))
    bdproto(NULL, self, predictors = x)
  },
  # Get Predictors
  get_predictors = function(self){
    return(self$predictors)
  },
  # Get scenario predictions
  get_scenarios = function(self, what = NULL){
    if(is.null(self$scenarios)) {
      return(self$scenarios)
    } else {
      if(what %in% names(self$scenarios)) return( self$scenarios[[what]] )
    }
  },
  # Plot the prediction
  plot = function(self, what = 'mean',...){
    raster::plot( self$scenarios[[what]] )
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
    cl <- sapply(self$scenarios, class)
    ras <- self$scenarios[[grep('raster', cl,ignore.case = T)]]

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
