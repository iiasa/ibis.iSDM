#' @include utils.R bdproto-biodiversitydistribution.R
NULL

#' Create distribution modelling procedure
#'
#' @description
#' This function creates an object that contains all the data, parameters and settings
#' for building an (integrated) species distribution model.
#' Key functions to add data are [add_biodiversity], [add_predictors],
#' [add_latent], [engine], [add_priors] and [add_offset]. This creates a
#' prototype [`BiodiversityDistribution`] object with its own functions.
#' After setting input data and parameters, model predictions can then be created
#' via the [train] function and predictions be created.
#' See **Details** for further functions available to modify or summarize the created object.
#'
#' **Note that any model requires at minimum a single added [biodiversity] dataset
#' as well as a specified [engine].**
#'
#' @param background Specification of the modelling background. Must be a
#' [`raster`], [`sf`] or [`extent`] object.
#' @param limits A [`raster`] or [`sf`] object that limits the prediction surface when
#' intersected with input data (Default: \code{NULL}).
#'
#' @details
#' This function creates a [`BiodiversityDistribution-class`] object that in itself contains
#' other functions and stores parameters and (pre-)processed data.
#' A full list of functions available can be queried via \code{names(object)}.
#' Some of the functions are not intended to be manipulated directly,
#'  but rather through convenience functions (e.g. [`object$set_predictors()`]).
#' Similarly other objects are stored in the [`BiodiversityDistribution-class`] object that
#' have their own functions as well and can be queried (e.g. [`names(object)`]). For a list of
#' functions see the reference documentation. By default,
#' if some datasets are not set, then a [`Waiver`] object is returned instead.
#'
#' The following objects can be stored:
#' * \code{object$biodiversity} A [`BiodiversityDatasetCollection`] object with the added biodiversity data.
#' * \code{object$engine} An [`engine`] object (e.g. [engine_inlabru]) with function depended on the added engine.
#' * \code{object$predictors} A [`PredictorDataset`] object with all set predictions.
#' * \code{object$priors} A [`PriorList`] object with all specified priors.
#' * \code{object$log} A [`Log`] object that captures.
#'
#' Useful high-level functions to address those objects are for instance:
#' * \code{object$show()} A generic summary of the [`BiodiversityDistribution-class`] object contents. Can also be called via [print].
#' * \code{object$get_biodiversity_equations()} Lists the equations used for each biodiversity dataset with given id. Defaults to all predictors.
#' * \code{object$get_biodiversity_types()} Lists the type of each specified biodiversity dataset with given id.
#' * \code{object$get_extent()} Outputs the [raster::extent] of the modelling region.
#' * \code{object$show_background_info()} Returns a [`list`] with the [raster::extent] and the [sp::proj4string].
#' * \code{object$get_extent_dimensions()} Outputs the [raster::extent] dimension by calling the [`extent_dimensions()`] function.
#' * \code{object$get_predictor_names()} Returns a [character] vector with the names of all added predictors.
#' * \code{object$get_prior_variables()} Returns a description of [`priors`] added.
#'
#' There are other functions as well but those are better accessed through their respective wrapper functions.
#' @returns [`BiodiversityDistribution-class`] object containing data for building a biodiversity distribution modelling problem.
#'
#' @seealso [`bdproto`] on the general definition of [`proto`] objects and in particular [`bdproto-biodiversitydistribution`].
#'
#' @references
#' * Fletcher, R.J., Hefley, T.J., Robertson, E.P., Zuckerberg, B., McCleery, R.A., Dorazio, R.M., (2019) A practical guide for combining data to model species distributions. Ecology 100, e02710. [https://doi.org/10.1002/ecy.2710](https://doi.org/10.1002/ecy.2710)
#' @aliases distribution
#' @exportMethod distribution
#' @name distribution
#'
#' @examples
#' \dontrun{
#'  # Load background raster
#'  background <- raster::raster(system.file("inst/extdata/europegrid_50km.tif",package = "ibis.iSDM"))
#'  Define model
#'  x <- distribution(background)
#'  x
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
