#' @include utils.R
NULL

#' Threshold a continuous prediction to a categorical layer
#'
#' @param obj A trained [`DistributionModel`] or alternatively a [`Raster`] object
#' @param method A specifc method for thresholding. One of 'fixed', 'mtp', 'percentile'
#' @param value A [`numeric`] value for thresholding if method is fixed (Default: NULL)
#' @param poi A [`sf`] object containing observational data used for model training
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
  function(obj, method = 'fixed', value = NULL,...) standardGeneric("threshold"))

#' Threshold with supplied DistributionModel object
#' @name threshold
#' @rdname threshold
#' @usage \S4method{threshold}{ANY, character}(obj, method)
methods::setMethod(
  "threshold",
  methods::signature(obj = "ANY", method = "character"),
  function(obj, method = 'mtp', value = NULL,...) {
    assertthat::assert_that("DistributionModel" %in% class(obj),
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
    poi <- sf::st_as_sf( obj$model$biodiversity[[1]]$observations)

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
                            layers = 'mean', poi = NULL) {
  assertthat::assert_that(is.Raster(obj),
                          is.character(method),
                          is.character(layers),
                          is.null(poi) || inherits(poi,'sf'),
                          is.null(value) || is.numeric(value)
  )
  # Get specified layers in obj
  wl <- which(names(obj) %in% layers)
  assertthat::assert_that(length(wl)>0,msg = 'Specified layers not found in object.')
  obj <- obj[[wl]]
  # Now apply threshold on each entry
  out <- raster::stack()
  for(i in names(obj)[wl]) out <- raster::addLayer(out, threshold(obj[[i]], method = method, value = value, poi = poi) )
  return(out)
}
#' @name threshold
#' @rdname threshold
#' @usage \S4method{threshold}{RasterBrick, character}(obj, method)
methods::setMethod("threshold",methods::signature(obj = "RasterBrick", method = "character"),.stackthreshold)
#' @usage \S4method{threshold}{RasterStack, character}(obj, method)
methods::setMethod("threshold",methods::signature(obj = "RasterStack", method = "character"),.stackthreshold)

#' @name threshold
#' @rdname threshold
#' @usage \S4method{threshold}{RasterLayer, character}(obj, method)
methods::setMethod(
  "threshold",
  methods::signature(obj = "RasterLayer", method = "character"),
  function(obj, method = 'fixed', value = NULL, poi = NULL) {
    assertthat::assert_that(is.Raster(obj),
                            inherits(obj,'RasterLayer'),
                            is.character(method),
                            is.null(poi) || inherits(poi,'sf'),
                            is.null(value) || is.numeric(value)
    )

    # Check that raster has at least a mean prediction in name
    if(!is.null(poi)) assertthat::assert_that(unique(sf::st_geometry_type(poi)) %in% c('POINT','MULTIPOINT'))

    # Get the raster layer
    raster_thresh <- obj

    # dismo::evaluate
    # dismo::threshold
    # out <- ras$mean >= quantile(ras$mean)[4]
    # out <- raster::cut(ras$mean, quantile(ras$mean))

    # If defined by type
    if(method != 'fixed'){
      if (!is.null(poi)) {
        pointVals <- raster::extract(raster_thresh, poi)
        # minimum training presence
        if (method == "mtp") {
          tr <- min( na.omit(pointVals) )
        }
        # percentile training threshold
        if (method == "percentile") {
          if(is.null(value)) value <- 0.1 # If value is not set, use 10%
          if(length(pointVals) < 10) {
            perc <- floor(length(pointVals) * (1-value))
          } else {
            perc <- ceiling(length(pointVals) * (1-value))
          }
          tr <- rev(sort(pointVals))[perc]
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
    # Finally threshold the raster
    raster_thresh[raster_thresh < tr] <- 0
    raster_thresh[raster_thresh >= tr] <- 1
    names(raster_thresh) <- paste0('threshold_',names(obj),'_',method)
    raster_thresh <- raster::asFactor(raster_thresh)

    # Return result
    return(raster_thresh)
  }
)
