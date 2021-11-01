#' @include utils.R bdproto-biodiversitydistribution.R
NULL

#' Setup biodiversity distribution modelling procedure
#'
#' @param background Specification of the modelling background. Must be a
#' [`raster`], [`sf`] or [`extent`] object
#' @param limits A [`raster`] or [`sf`] object that limits the prediction surface when
#' intersected with input data (Default: NULL).
#'
#' @details TBD. Say something about PPMs, INLA and co.
#'
#' Option with constrain -> use reference
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
                    function(background, limits = NULL) standardGeneric("distribution"))

#' @name distribution
#' @usage \S4method{distribution}{raster}(background)
#' @rdname distribution
methods::setMethod(
  "distribution",
  methods::signature(background = "Raster"),
  function(background, limits = NULL) {
    assertthat::assert_that(!missing(background) || !exists('background'),
                            inherits(limits,'Raster') || inherits(limits, 'sf') || inherits(limits, 'Spatial') || is.null(limits),
                            msg = 'No background file supplied!')
    # Check that arguments are valid
    assertthat::assert_that( inherits(background,'Raster')  )

    # Convert raster to dissolved polygons to get a study boundary
    newbg <- sf::st_as_sf(
      raster::rasterToPolygons(background, dissolve = TRUE)
    )

    # Rerun the distribution call with the object
    distribution(newbg, limits)
  })

#' @name distribution
#' @usage \S4method{distribution}{sf}(background)
#' @rdname distribution
methods::setMethod(
  "distribution",
  methods::signature(background = "sf"),
  function(background, limits = NULL) {
    # Check that arguments are valid
    assertthat::assert_that(!missing(background) || !exists('background'),
                            inherits(limits,'Raster') || inherits(limits, 'sf') || inherits(limits, 'Spatial') || is.null(limits),
                            msg = 'No background file supplied!')
    assertthat::assert_that(
      inherits(background,'sf'),
      unique(st_geometry_type(background)) %in% c('MULTIPOLYGON','POLYGON')
    )

    # Messager
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Creating distribution object...')

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
      if(st_crs(limits) != st_crs(background)) limits <- st_transform(limits, background)
      # Ensure that limits is intersecting the background
      if(suppressMessages(length( st_intersects(limits, background)))==0) { limits <- NULL; warning('Provided limits do not intersect the background!') }

      # Get fir column and rename
      limits <- limits[,1]; names(limits) <- c('limit','geometry')
    }

    # Convert to waiver if NULL
    if(is.null(limits)) limits <- new_waiver()

    # Create BiodiversityDistribution object
    bdproto(NULL, BiodiversityDistribution,
            background = background,
            limits = limits
            )
  })
