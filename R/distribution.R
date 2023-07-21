#' @include utils.R bdproto-biodiversitydistribution.R
NULL

#' Create distribution modelling procedure
#'
#' @description
#' This function creates an object that contains all the data, parameters and settings
#' for building an (integrated) species distribution model.
#' Key functions to add data are [`add_biodiversity_poipo`] and the like, [`add_predictors`],
#' [`add_latent_spatial`], [`engine_glmnet`] or similar, [`add_priors`] and [`add_offset`].
#' It creates a prototype [`BiodiversityDistribution`] object with its own functions.
#' After setting input data and parameters, model predictions can then be created
#' via the [train] function and predictions be created.
#'
#' Additionally, it is possible to specify a \code{"limit"} to any predictions conducted on
#' the background. This can be for instance a buffered layer by a certain dispersal distance (Cooper and Soberon, 2018)
#' or a categorical layer representing biomes or soil conditions.
#' Another option is to create a constraint by constructing a minimum convex polygon (MCP) using
#' the supplied biodiversity data. This option can be enabled by setting
#' \code{"limits_method"} to \code{"mcp"}. It is also possible to provide a small buffer
#' to constructed MCP that way.
#' See the frequently asked question (FAQ) section on the homepage for more information.
#'
#' See **Details** for a description of the internal functions available
#' to modify or summarize data within the created object.
#'
#' **Note that any model requires at minimum a single added biodiversity dataset
#' as well as a specified engine.**
#'
#' @param background Specification of the modelling background. Must be a
#' [`SpatRaster`] or [`sf`] object.
#' @param limits A [`SpatRaster`] or [`sf`] object that limits the prediction surface when
#' intersected with input data (Default: \code{NULL}).
#' @param limits_method A [`character`] of the method used for hard limiting a projection.
#' Available options are \code{"none"} (Default), \code{"zones"} or \code{"mcp"}.
#' @param mcp_buffer A [`numeric`] distance to buffer the mcp (Default \code{0}). Only used if
#'  \code{"mcp"} is used.
#' @param limits_clip [`logical`] Should the limits clip all predictors before fitting
#' a model (\code{TRUE}) or just the prediction (\code{FALSE}, default).
#'
#' @details
#' This function creates a [`BiodiversityDistribution-class`] object that in itself contains
#' other functions and stores parameters and (pre-)processed data.
#' A full list of functions available can be queried via \code{"names(object)"}.
#' Some of the functions are not intended to be manipulated directly,
#'  but rather through convenience functions (e.g. \code{"object$set_predictors()"}).
#' Similarly other objects are stored in the [`BiodiversityDistribution-class`] object that
#' have their own functions as well and can be queried (e.g. \code{"names(object)"}). For a list of
#' functions see the reference documentation. By default,
#' if some datasets are not set, then a \code{"Waiver"} object is returned instead.
#'
#' The following objects can be stored:
#' * \code{object$biodiversity} A [`BiodiversityDatasetCollection`] object with the added biodiversity data.
#' * \code{object$engine} An \code{"engine"} object (e.g. [`engine_inlabru()`]) with function depended on the added engine.
#' * \code{object$predictors} A [`PredictorDataset`] object with all set predictions.
#' * \code{object$priors} A [`PriorList`] object with all specified priors.
#' * \code{object$log} A [`Log`] object that captures.
#'
#' Useful high-level functions to address those objects are for instance:
#' * \code{object$show()} A generic summary of the [`BiodiversityDistribution-class`] object contents. Can also be called via [print].
#' * \code{object$get_biodiversity_equations()} Lists the equations used for each biodiversity dataset with given id. Defaults to all predictors.
#' * \code{object$get_biodiversity_types()} Lists the type of each specified biodiversity dataset with given id.
#' * \code{object$get_extent()} Outputs the [terra::ext] of the modelling region.
#' * \code{object$show_background_info()} Returns a [`list`] with the [terra::ext] and the [terra::crs].
#' * \code{object$get_extent_dimensions()} Outputs the [terra::ext] dimension by calling the \code{"extent_dimensions()"} function.
#' * \code{object$get_predictor_names()} Returns a [character] vector with the names of all added predictors.
#' * \code{object$get_prior_variables()} Returns a description of [`priors`] added.
#'
#' There are other functions as well but those are better accessed through their respective wrapper functions.
#'
#' @returns [`BiodiversityDistribution-class`] object containing data for building a biodiversity distribution modelling problem.
#'
#' @seealso \code{"bdproto"} on the general definition of [`proto`] objects and in particular [`BiodiversityDistribution`].
#'
#' @references
#' * Fletcher, R.J., Hefley, T.J., Robertson, E.P., Zuckerberg, B., McCleery, R.A., Dorazio, R.M., (2019) A practical guide for combining data to model species distributions. Ecology 100, e02710. https://doi.org/10.1002/ecy.2710
#' * Cooper, Jacob C., and Jorge Soberón. "Creating individual accessible area hypotheses improves stacked species distribution model performance." Global Ecology and Biogeography 27, no. 1 (2018): 156-165.
#' @aliases distribution
#' @exportMethod distribution
#' @name distribution
#'
#' @examples
#' \dontrun{
#'  # Load background raster
#'  background <- terra::rast(system.file("inst/extdata/europegrid_50km.tif",package = "ibis.iSDM"))
#'  # Define model
#'  x <- distribution(background)
#'  x
#'  # Show names of the functions within the object
#'  names(x)
#'
#' }
#' @export
methods::setGeneric("distribution",
                    signature = methods::signature("background"),
                    function(background, limits = NULL, limits_method = "none", mcp_buffer = 0,limits_clip = FALSE) standardGeneric("distribution"))

#' @name distribution
#' @usage \S4method{distribution}{SpatRaster,ANY,character,numeric,logical}(background,limits,limits_method,mcp_buffer,limits_clip)
#' @rdname distribution
methods::setMethod(
  "distribution",
  methods::signature(background = "SpatRaster"),
  function(background, limits = NULL, limits_method = "none", mcp_buffer = 0,limits_clip = FALSE) {
    assertthat::assert_that(!missing(background) || !exists('background'),
                            inherits(limits,'SpatRaster') || inherits(limits, 'sf') || inherits(limits, 'Spatial') || is.null(limits),
                            is.character(limits_method),
                            is.numeric(mcp_buffer),
                            is.logical(limits_clip),
                            msg = 'No background file supplied or limits misspecified!')
    # Check that arguments are valid
    assertthat::assert_that( inherits(background,'SpatRaster')  )

    # Convert raster to dissolved polygons to get a study boundary
    newbg <- sf::st_as_sf(
      terra::as.polygons(background, dissolve = TRUE)
    )

    # Rerun the distribution call with the object
    distribution(newbg, limits, limits_method, mcp_buffer, limits_clip)
  })

#' @name distribution
#' @usage \S4method{distribution}{sf,ANY,character,numeric,logical}(background,limits,limits_method,mcp_buffer,limits_clip)
#' @rdname distribution
methods::setMethod(
  "distribution",
  methods::signature(background = "sf"),
  function(background, limits = NULL, limits_method = "none", mcp_buffer = 0, limits_clip = FALSE) {
    # Check that arguments are valid
    assertthat::assert_that(!missing(background) || !exists('background'),
                            inherits(limits,'SpatRaster') || inherits(limits, 'sf') || inherits(limits, 'Spatial') || is.null(limits),
                            is.character(limits_method),
                            is.numeric(mcp_buffer),
                            is.logical(limits_clip),
                            msg = 'No background file supplied!')
    assertthat::assert_that(
      inherits(background,'sf'),
      unique(st_geometry_type(background)) %in% c('MULTIPOLYGON','POLYGON')
    )

    # Small checks on alternative limit functionalities
    limits_method <- match.arg(limits_method, c("none","zones", "mcp"), several.ok = FALSE)
    assertthat::assert_that(mcp_buffer>=0, msg = "Buffered mcp distance has to be positive!")

    # Messenger
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Creating distribution object...')

    # Convert limits if provided
    if(!is.null(limits)){
      # Set methods in case a layer was supplied
      limits_method <- "zones"
      # Convert to polygon if raster
      if(inherits(limits,'SpatRaster')){
        assertthat::assert_that(terra::is.factor(limits),
                                msg = 'Provided limit raster needs to be ratified (categorical)!')
        limits <- sf::st_as_sf( terra::as.polygons(limits, dissolve = TRUE) ) |> sf::st_cast("MULTIPOLYGON")
      }
      assertthat::assert_that(inherits(limits, "sf"))

      # Ensure that limits has the same projection as background
      if(sf::st_crs(limits) != sf::st_crs(background)) limits <- sf::st_transform(limits, background)
      # Ensure that limits is intersecting the background
      if(is.Raster(background)){
        if(suppressMessages(length( sf::st_intersects(limits, terra::as.polygons(background) |> sf::st_as_sf()) )) == 0 ) { limits <- NULL; warning('Provided limits do not intersect the background!') }
      } else {
        if(suppressMessages(length( sf::st_intersects(limits, background |> sf::st_as_sf()) )) == 0 ) { limits <- NULL; warning('Provided limits do not intersect the background!') }
      }

      # Get fir column and rename
      limits <- limits[,1]; names(limits) <- c('limit','geometry')
      limits <- list(layer = limits, "limits_method" = "zones",
                     "mcp_buffer" = mcp_buffer, "limits_clip" = limits_clip)
    } else if(limits_method == "mcp"){
      # Specify the option to calculate a mcp based on the added data.
      # This is done directly in train.
      limits <- list("layer" = NULL, "limits_method" = "mcp",
                     "mcp_buffer" = mcp_buffer, "limits_clip" = limits_clip)

    } else {
      # Convert to waiver if NULL
      if(is.null(limits)) limits <- new_waiver()
    }

    # Create BiodiversityDistribution object
    bdproto(NULL, BiodiversityDistribution,
            background = background,
            limits = limits
            )
  })
