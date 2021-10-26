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
#' 'fixed' applies a single determined threshold
#' 'mtp' minimum training presence find and sets the lowest predicted suitability for any occurrence point
#' 'percentile' For a percentile threshold
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
    # Match to correct spelling mistakes
    method <- match.arg(tolower(method), c('fixed','mtp','percentile'), several.ok = FALSE)
    # Check that provided method is supported
    assertthat::assert_that(
      method %in% c('fixed','mtp','percentile'),
      msg = 'Method not yet supported.'
    )

    # Get all point data in distribution model
    # FIXME: Adapt to more/ combined datasets
    poi <- sf::st_as_sf( obj$model$biodiversity[[1]]$observations, coords = c('x','y') )
    poi <- subset(poi, observed > 0) # Remove any eventual absence data

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
  function(obj, method = 'fixed', value = NULL, poi = NULL, return_threshold = FALSE) {
    assertthat::assert_that(is.Raster(obj),
                            inherits(obj,'RasterLayer'),
                            is.character(method),
                            is.null(poi) || inherits(poi,'sf'),
                            is.null(value) || is.numeric(value)
    )
    # Match to correct spelling mistakes
    method <- match.arg(tolower(method), c('fixed','mtp','percentile'), several.ok = FALSE)

    # Check that raster has at least a mean prediction in name
    if(!is.null(poi)) assertthat::assert_that(unique(sf::st_geometry_type(poi)) %in% c('POINT','MULTIPOINT'))
    assertthat::assert_that(hasName(poi, 'observed'))
    poi_pres <- subset(poi, observed > 0) # Remove any eventual absence data

    # Get the raster layer
    raster_thresh <- obj

    # If defined by type
    if(method != 'fixed'){
      if (!is.null(poi_pres)) {
        pointVals <- raster::extract(raster_thresh, poi_pres)
        # minimum training presence
        if (method == "mtp") {
          tr <- min( na.omit(pointVals) )
        } else
        # percentile training threshold
        if (method == "percentile") {
          if(is.null(value)) value <- 0.1 # If value is not set, use 10%
          if(length(pointVals) < 10) {
            perc <- floor(length(pointVals) * (1-value))
          } else {
            perc <- ceiling(length(pointVals) * (1-value))
          }
          tr <- rev(sort(pointVals))[perc]
        } else
          # Optimized True Skill Statistic
          if(method == 'optiTSS'){
            # TODO: Implement optimal TSS selection
            # Make sure that poi includes [0,1]
            # https://github.com/r-forge/modeva/blob/master/pkg/R/optiThresh.R
            opt <- optiThresh(obs = dat[, myspecies], pred = dat[, "BART_P"], measures = "TSS",
                              optimize = "each", interval = 1e-04)
          }
      } else {
       # No point data defined. Raise error
       stop('Training points needed for this function to work.')
      }
    } else {
      # Fixed threshold. Confirm to be set
      assertthat::assert_that(is.numeric(value))
      tr <- value
    }
    # -- Threshold -- #
    # Security check
    assertthat::assert_that(is.numeric(tr))
    if(return_threshold){
      return(tr)
    } else {
      # Finally threshold the raster
      raster_thresh[raster_thresh < tr] <- 0
      raster_thresh[raster_thresh >= tr] <- 1
      names(raster_thresh) <- paste0('threshold_',names(obj),'_',method)
      raster_thresh <- raster::asFactor(raster_thresh)
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
  function(obj, method = 'mtp', value = NULL, poi = NULL, return_threshold = TRUE) {
    # Assert that predicted raster is present
    assertthat::assert_that( is.Raster(obj$get_model()$get_data('prediction'))  )
    tr <- threshold(  obj = obj$get_model()$get_data('prediction'),
                      method = method,
                      value = value,
                      poi = poi,
                      return_threshold = return_threshold
                    )
    bdproto(NULL, obj, threshold = tr)
  }
)
