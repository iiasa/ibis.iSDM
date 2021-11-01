#' @include utils.R
NULL

#' Threshold a continuous prediction to a categorical layer
#'
#' @param obj A trained [`DistributionModel`] or alternatively a [`Raster`] object
#' @param method A specifc method for thresholding. One of 'fixed', 'mtp', 'percentile'
#' @param value A [`numeric`] value for thresholding if method is fixed (Default: NULL)
#' @param poi A [`sf`] object containing observational data used for model training
#' @param return_threshold Should threshold value be returned instead (Default: FALSE)
#' @details
#' The following options are currently implemented:
#' * 'fixed' = applies a single pre-determined threshold
#' * 'mtp' = minimum training presence find and sets the lowest predicted suitability for any occurrence point
#' * 'percentile' = For a percentile threshold
#' * 'TSS' = Determines the optimal TSS (True Skill Statistic). Requires the ['modEvA'] package
#' * 'kappa' = Determines the optimal kappa value (Kappa). Requires the ['modEvA'] package
#' * 'F1score' = Determines the optimal F1score (also known as Sorensen similarity). Requires the ['modEvA'] package
#' * 'F1score' = Determines the optimal sensitivity of presence records. Requires the ['modEvA'] package
#' * 'Sensitivity' = Determines the optimal sensitivity of presence records. Requires the ['modEvA'] package
#' * 'Specificity' = Determines the optimal sensitivity of presence records. Requires the ['modEvA'] package
#' @name threshold
#' @examples
#' \dontrun{
#' print('test')
#' }
#' @export
NULL

#' @name threshold
#' @rdname threshold
#' @exportMethod threshold
#' @export
methods::setGeneric(
  "threshold",
  signature = methods::signature("obj", "method", "value"),
  function(obj, method = 'fixed', value = NULL, poi = NULL, return_threshold = FALSE, ...) standardGeneric("threshold"))

#' Generic threshold with supplied DistributionModel object
#' @name threshold
#' @rdname threshold
#' @usage \S4method{threshold}{ANY}(obj)
methods::setMethod(
  "threshold",
  methods::signature(obj = "ANY"),
  function(obj, method = 'mtp', value = NULL, return_threshold = FALSE, ...) {
    assertthat::assert_that(inherits(obj,c('GDB-Model','BART-Model','INLA-Model','STAN-Model')),
                            is.character(method),
                            is.null(value) || is.numeric(value)
    )
    # Get raster
    ras <- obj$get_data('prediction')

    # Check that the object actually contains a prediction
    assertthat::assert_that(
      is.Raster(ras),
      !is.Waiver(ras),
      msg = 'No fitted prediction in object!'
    )
    # Matching for correct method
    method <- match.arg(method, c('fixed','mtp','percentile',
                                           'TSS','kappa','F1score','Sensitivity','Specificity'), several.ok = FALSE)

    # Get all point data in distribution model
    # FIXME: Adapt to more/ combined datasets
    poi <- sf::st_as_sf( obj$model$biodiversity[[1]]$observations, coords = c('x','y') )

    # If TSS or kappa is chosen, check whether there is poipa data among the sources
    if(!any(sapply(mod1$model$biodiversity, function(x) 0 %in% x$observations[,'observed'])) & method %in% c('TSS','kappa','F1score','Sensitivity','Specificity')){
      if(getOption('ibis.setupmessages')) myLog('[Threshold]','red','Threshold method needs absence-data. Generating some now...')
      bg <- raster::rasterize(obj$model$background, emptyraster(obj$get_data('prediction')))
      abs <- create_pseudoabsence(
        env = obj$model$predictors,
        presence = poi,
        bias = obj$settings$get('bias_variable'),
        template = bg,
        npoints = ifelse(ncell(bg)<10000,ncell(bg),10000),
        replace = TRUE
      )
      abs <- subset(abs, select = c('x','y'));abs$observed <- 0
      poi <- rbind(poi, st_as_sf(abs, coords = c('x','y')))
    }

    # Now self call threshold
    out <- threshold(ras, method = method, value = value, poi = poi, ...)
    assertthat::assert_that(is.Raster(out))
    # Add result to new obj
    new_obj <- obj
    if(inherits(out,'RasterLayer')){
      new_obj <- new_obj$set_data(names(out), out)
    } else if(inherits(out,'RasterStack')) {
      # When stack loop through and add
      for(n in names(out)){ new_obj <- new_obj$set_data(n, out[[n]]) }
    }
    # Return altered object
    return(new_obj)
  }
)

#' @noRd
#' @keywords noexport
.stackthreshold <- function(obj, method = 'fixed', value = NULL,
                            poi = NULL, return_threshold = FALSE) {
  assertthat::assert_that(is.Raster(obj),
                          is.character(method),
                          inherits(poi,'sf'),
                          is.null(value) || is.numeric(value)
  )
  # Apply threshold on each entry
  if(return_threshold){
    # Return the threshold directly
    out <- vector()
    for(i in names(obj)) out <- c(out, threshold(obj[[i]], method = method,
                                                                value = value, poi = poi, return_threshold = return_threshold) )
    names(out) <- names(obj)
  } else {
    # Return the raster instead
    out <- raster::stack()
    for(i in names(obj)) out <- raster::addLayer(out, threshold(obj[[i]], method = method,
                                                                value = value, poi = poi, return_threshold = return_threshold) )
  }
  return(out)
}

#' @name threshold
#' @rdname threshold
#' @inheritParams threshold
#' @usage \S4method{threshold}{RasterBrick}(obj)
methods::setMethod("threshold",methods::signature(obj = "RasterBrick"),.stackthreshold)
#' @usage \S4method{threshold}{RasterStack}(obj)
methods::setMethod("threshold",methods::signature(obj = "RasterStack"),.stackthreshold)

#' @name threshold
#' @rdname threshold
#' @usage \S4method{threshold}{RasterLayer}(obj)
methods::setMethod(
  "threshold",
  methods::signature(obj = "RasterLayer"),
  function(obj, method = 'fixed', value = NULL, poi = NULL, return_threshold = FALSE, plot = FALSE) {
    assertthat::assert_that(is.Raster(obj),
                            inherits(obj,'RasterLayer'),
                            is.character(method),
                            is.null(value) || is.numeric(value)
    )
    # If poi is set, try to convert sf
    if(!is.null(poi)) try({poi <- sf::st_as_sf(poi)})
    assertthat::assert_that(is.null(poi) || inherits(poi,'sf'))

    # Match to correct spelling mistakes
    method <- match.arg(method, c('fixed','mtp','percentile',
                                           'TSS','kappa','F1score','Sensitivity','Specificity'), several.ok = FALSE)

    # Check that raster has at least a mean prediction in name
    if(!is.null(poi)) assertthat::assert_that(unique(sf::st_geometry_type(poi)) %in% c('POINT','MULTIPOINT'))
    assertthat::assert_that(hasName(poi, 'observed'))
    poi_pres <- subset(poi, observed > 0) # Remove any eventual absence data

    # Get the raster layer
    raster_thresh <- obj

    # Specify by type:
    if(method == "fixed"){
      # Fixed threshold. Confirm to be set
      assertthat::assert_that(is.numeric(value),msg = 'Fixed value is missing!')
      tr <- value
    } else if(method == "mtp"){
      # minimum training presence
      pointVals <- raster::extract(raster_thresh, poi_pres) # Extract point only estimates
      # Minimum threshold
      tr <- min( na.omit(pointVals) )

    } else if(method == "percentile"){
      pointVals <- raster::extract(raster_thresh, poi_pres) # Extract point only estimates
      # percentile training threshold
      if(is.null(value)) value <- 0.1 # If value is not set, use 10%
      if(length(pointVals) < 10) {
        perc <- floor(length(pointVals) * (1-value))
      } else {
        perc <- ceiling(length(pointVals) * (1-value))
      }
      tr <- rev(sort(pointVals))[perc] # Percentile threshold

    } else {
      # Optimized threshold statistics using the modEvA package
      # FIXME: Could think of porting these functions but too much effort for now. Rather have users install the package here
      check_package("modEvA")
      assertthat::assert_that('modEvA' %in% installed.packages()[,1])
      # Assure that point data is correctly specified
      assertthat::assert_that(inherits(poi, 'sf'), hasName(poi, 'observed'))
      poi$observed <- ifelse(poi$observed>1,1,poi$observed) # Ensure that observed is <=1
      assertthat::assert_that(all( unique(poi$observed) %in% c(0,1) ))

      # Re-extract point vals but with the full dataset
      pointVals <- raster::extract(raster_thresh, poi)
      assertthat::assert_that(length(pointVals)>2)
      # Calculate the optimal thresholds
      suppressWarnings(
        opt <- modEvA::optiThresh(obs = poi$observed, pred = pointVals,
                                  measures = c("TSS","kappa","F1score","Misclass","Omission","Commission",
                                               "Sensitivity","Specificity"),
                                  optimize = "each", plot = plot)
      )
      if(method %in% opt$optimals.each$measure){
        tr <- opt$optimals.each$threshold[which(opt$optimals.each$measure==method)]
      } else {
        # Returning a collection of them as vector
        tr <- opt$optimals.each$threshold; names(tr) <- opt$optimals.each$measure
      }
    }
    # Security check
    assertthat::assert_that(is.numeric(tr) || is.vector(tr))

    # -- Threshold -- #
    if(return_threshold){
      names(tr) <- method
      return(tr)
    } else {
      # Finally threshold the raster
      raster_thresh[raster_thresh < tr[1]] <- 0
      raster_thresh[raster_thresh >= tr[1]] <- 1
      names(raster_thresh) <- paste0('threshold_',names(obj),'_',method)
      raster_thresh <- raster::asFactor(raster_thresh)
      # Assign attributes
      base::attr(raster_thresh, 'method') <- method
      base::attr(raster_thresh, 'threshold') <- tr
    }
    # Return result
    return(raster_thresh)
  }
)

#### For scenarios ####

#' Thresholds in scenario estimation
#'
#' @description For [`BiodiversityScenario`] objects store a threshold attribute in
#' the scenario object
#' @name threshold
#' @inheritParams threshold
#' @rdname threshold
#' @usage \S4method{threshold}{BiodiversityScenario}(obj)
methods::setMethod(
  "threshold",
  methods::signature(obj = "BiodiversityScenario"),
  function(obj, ...) {
    # Assert that predicted raster is present
    assertthat::assert_that( is.Raster(obj$get_model()$get_data('prediction')) )
    # Check that a threshold layer is available and get the methods and data from it
    assertthat::assert_that( grep('threshold', obj$get_model()$show_rasters())>0 ,
                             msg = 'Call \' threshold \' for prediction first!')
    # Get threshold layer
    ras_tr <- obj$get_model()$get_data( grep('threshold', obj$get_model()$show_rasters(),value = TRUE) )
    tr <- attr(ras_tr, 'threshold')
    names(tr) <- attr(ras_tr, 'method')
    # Otherwise sample?
    # tr <- threshold(  obj = obj$get_model()$get_data('prediction'),
    #                   method = method,
    #                   value = value,
    #                   poi = poi,
    #                   return_threshold = return_threshold
    #                 )
    bdproto(NULL, obj, threshold = tr)
  }
)
