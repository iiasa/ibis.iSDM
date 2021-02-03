#' @include utils.R bdproto-biodiversitydistribution.R
NULL

#' Setup biodiversity distribution modelling procedure
#'
#' @param background Specification of the modelling background. Must be a
#' [`raster`], [`sf`] or [`extent`] object
#' @param formula Model formula. Has to be [`character`] or [`formula`] object (optional)
#' @param ... not used.
#'
#' @details TBD. Say something about PPMs, INLA and co
#'
#' @return [`BiodiversityDistribution-class`] object containing
#'   data for building a biodiversity distribution modelling problem.
#'
#' @seealso [bdproto]
#'
#' @aliases distribution
#'
#' @exportMethod distribution
#'
#' @name distribution
#'
#' @examples
#' \dontrun{
#' print('test')
#' }
#' @export
methods::setGeneric("distribution",
                    signature = methods::signature("background"),
                    function(background, ...) standardGeneric("distribution"))

#' @name distribution
#' @import rgeos
#' @usage \S4method{distribution}{raster}(background, ...)
#' @rdname pdistribution
methods::setMethod(
  "distribution",
  methods::signature(background = "Raster"),
  function(background, ...) {
    # Check that arguments are valid
    assertthat::assert_that( inherits(background,'Raster')  )

    # Convert raster to dissolved polygons to get a study boundary
    newbg <- sf::st_as_sf(
      raster::rasterToPolygons(background, dissolve = TRUE)
    )

    # Rerun the distribution call with the object
    distribution(newbg, ...)
  })

#' @name distribution
#' @usage \S4method{distribution}{sf}(background, ...)
#' @rdname distribution
methods::setMethod(
  "distribution",
  methods::signature(background = "sf"),
  function(background, ...) {
    # Check that arguments are valid
    assertthat::assert_that(
      inherits(background,'sf'),
      unique(st_geometry_type(background)) %in% c('MULTIPOLYGON','POLYGON')
    )

    # Create BiodiversityDistribution object
    bdproto(NULL, BiodiversityDistribution,
            background = background,
            ...
            )
  })
