#' @include utils.R bdproto-biodiversityscenario.R
NULL

#' Create a new scenario based on trained model parameters
#'
#' @param fit A trained [`BiodiversityDistribution`] model
#' @param limits A [`raster`] or [`sf`] object that limits the projection surface when
#' intersected with the prediction data (Default: NULL).
#' @aliases scenario
#' @exportMethod scenario
#' @name scenario
#'
#' @examples
#' \dontrun{
#' print('test')
#' }
#' @export
methods::setGeneric("scenario",
                    signature = methods::signature("fit"),
                    function(fit, limits = NULL) standardGeneric("scenario"))

#' @name scenario
#' @usage \S4method{scenario}{ANY}(fit)
#' @rdname scenario
methods::setMethod(
  "scenario",
  methods::signature(fit = "ANY"),
  function(fit, limits = NULL) {
    # Check that arguments are valid
    assertthat::assert_that(!missing(fit) || inherits(fit,'DistributionModel'),
                            inherits(limits,'Raster') || inherits(limits, 'sf') || inherits(limits, 'Spatial') || is.null(limits),
                            msg = 'No trained model supplied!')

    # Get model object name and id
    modelobject <- deparse(substitute(fit))
    modelid <- fit$id

    # Convert limits if provided
    if(!is.null(limits)){
      # Convert to polygon if raster
      if(inherits(limits,'Raster')){
        if(is.null(levels(limits))) stop('Provided limit raster needs to be ratified (categorical)!')
        # Remove 0 from ratified raster assuming this is no-data
        limits[limits == 0] <- NA
        limits <- sf::st_as_sf( raster::rasterToPolygons(limits, n = 16, dissolve = TRUE) )
      }
      # Ensure that limits has the same projection as background
      if(sf::st_crs(limits) != sf::st_crs(fit$model$background)) limits <- sf::st_transform(limits, fit$model$background)
      # Ensure that limits is intersecting the background
      if(suppressMessages(length( sf::st_intersects(limits, fit$model$background)))==0) { limits <- NULL; warning('Provided limits do not intersect the background!') }

      # Get fir column and rename
      limits <- limits[,1]; names(limits) <- c('limit','geometry')
    }
    if(is.null(limits)){
      # Convert to waiver if NULL
      limits <- new_waiver()
    }
    # Create BiodiversityScenario object
    bdproto(NULL, BiodiversityScenario,
            modelobject = modelobject,
            modelid = modelid,
            limits = limits
    )
  })
