#' @include class-biodiversitydistribution.R
NULL

#' Create distribution modelling procedure
#'
#' @description This function creates an object that contains all the data,
#' parameters and settings for building an (integrated) species distribution
#' model. Key functions to add data are [`add_biodiversity_poipo`] and the like,
#' [`add_predictors`], [`add_latent_spatial`], [`engine_glmnet`] or similar,
#' [`add_priors`] and [`add_offset`]. It creates a prototype
#' [`BiodiversityDistribution`] object with its own functions. After setting
#' input data and parameters, model predictions can then be created via the
#' [train] function and predictions be created.
#'
#' Additionally, it is possible to specify a \code{"limit"} to any predictions
#' conducted on the background. This can be for instance a buffered layer by a
#' certain dispersal distance (Cooper and Soberon, 2018) or a categorical layer
#' representing biomes or soil conditions. Another option is to create a
#' constraint by constructing a minimum convex polygon (MCP) using the supplied
#' biodiversity data. This option can be enabled by setting
#' \code{"limits_method"} to \code{"mcp"}. It is also possible to provide a
#' small buffer to constructed MCP that way. See the frequently asked question
#' (FAQ) section on the homepage for more information.
#'
#' See **Details** for a description of the internal functions available to
#' modify or summarize data within the created object.
#'
#' **Note that any model requires at minimum a single added biodiversity dataset
#' as well as a specified engine.**
#'
#' @param background Specification of the modelling background. Must be a
#' [`SpatRaster`] or [`sf`] object.
#' @param limits A [`SpatRaster`], [`sf`] or [`stars`] object that limits the prediction
#' surface when intersected with input data (Default: \code{NULL}). In case of a
#' [`stars`] object the first factorized time entry is taken.
#' @param limits_method A [`character`] of the method used for hard limiting a
#' projection. Available options are \code{"none"} (Default), \code{"zones"}
#' or \code{"mcp"}. See also [`add_limits_extrapolation()`].
#' @param mcp_buffer A [`numeric`] distance to buffer the mcp (Default \code{0}).
#' Only used if \code{"mcp"} is used.
#' @param limits_clip [`logical`] Should the limits clip all predictors before
#' fitting a model (\code{TRUE}) or just the prediction (\code{FALSE}, default).
#'
#' @details This function creates a [`BiodiversityDistribution-class`] object
#' that in itself contains other functions and stores parameters and
#' (pre-)processed data. A full list of functions available can be queried via
#' \code{"names(object)"}. Some of the functions are not intended to be
#' manipulated directly, but rather through convenience functions (e.g.
#' \code{"object$set_predictors()"}). Similarly other objects are stored in the
#' [`BiodiversityDistribution-class`] object that have their own functions as
#' well and can be queried (e.g. \code{"names(object)"}). For a list of
#' functions see the reference documentation. By default, if some datasets are
#' not set, then a \code{"Waiver"} object is returned instead.
#'
#' The following objects can be stored:
#' * \code{object$biodiversity} A [`BiodiversityDatasetCollection`] object with
#' the added biodiversity data.
#' * \code{object$engine} An \code{"engine"} object (e.g. [`engine_inlabru()`])
#' with function depended on the added engine.
#' * \code{object$predictors} A [`PredictorDataset`] object with all set predictions.
#' * \code{object$priors} A [`PriorList`] object with all specified priors.
#' * \code{object$log} A [`Log`] object that captures.
#'
#' Useful high-level functions to address those objects are for instance:
#' * \code{object$show()} A generic summary of the [`BiodiversityDistribution-class`]
#' object contents. Can also be called via [print].
#' * \code{object$get_biodiversity_equations()} Lists the equations used for each
#' biodiversity dataset with given id. Defaults to all predictors.
#' * \code{object$get_biodiversity_types()} Lists the type of each specified
#' biodiversity dataset with given id.
#' * \code{object$get_extent()} Outputs the [terra::ext] of the modelling region.
#' * \code{object$show_background_info()} Returns a [`list`] with the [terra::ext]
#' and the [terra::crs].
#' * \code{object$get_extent_dimensions()} Outputs the [terra::ext] dimension by
#' calling the \code{"extent_dimensions()"} function.
#' * \code{object$get_predictor_names()} Returns a [character] vector with the
#' names of all added predictors.
#' * \code{object$get_prior_variables()} Returns a description of [`priors`] added.
#'
#' There are other functions as well but those are better accessed through their
#' respective wrapper functions.
#'
#' @returns [`BiodiversityDistribution-class`] object containing data for building
#' a biodiversity distribution modelling problem.
#'
#' @seealso [`BiodiversityDistribution`] and other classes.
#'
#' @references
#' * Fletcher, R.J., Hefley, T.J., Robertson, E.P., Zuckerberg, B., McCleery, R.A.,
#' Dorazio, R.M., (2019) A practical guide for combining data to model species
#' distributions. Ecology 100, e02710. https://doi.org/10.1002/ecy.2710
#' * Cooper, Jacob C., and Jorge Soberón. "Creating individual accessible area
#' hypotheses improves stacked species distribution model performance." Global
#' Ecology and Biogeography 27, no. 1 (2018): 156-165.
#'
#' @examples
#' # Load background raster
#' background <- terra::rast(system.file("extdata/europegrid_50km.tif",package = "ibis.iSDM"))
#' # Define model
#' x <- distribution(background)
#' x
#'
#' @name distribution
NULL

#' @rdname distribution
#' @export
methods::setGeneric("distribution",
                    signature = methods::signature("background"),
                    function(background, limits = NULL, limits_method = "none",
                             mcp_buffer = 0,limits_clip = FALSE) standardGeneric("distribution"))

#' @rdname distribution
methods::setMethod(
  "distribution",
  methods::signature(background = "SpatRaster"),
  function(background, limits = NULL, limits_method = "none", mcp_buffer = 0,limits_clip = FALSE) {
    assertthat::assert_that(!missing(background) || !exists('background'),
                            inherits(limits,'SpatRaster') || inherits(limits, 'sf') ||
                              inherits(limits, 'Spatial') || is.null(limits) || inherits(limits, "stars"),
                            is.character(limits_method),
                            is.numeric(mcp_buffer),
                            is.logical(limits_clip),
                            msg = 'No background file supplied or limits misspecified!')
    # Check that arguments are valid
    assertthat::assert_that( inherits(background,'SpatRaster')  )

    # Check that provided background has a valid crs
    if(is.na(terra::crs(background,describe=TRUE)[['code']])){
      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Setup]','red','Provide a background with a valid projection!')
    }

    # Make an error check that units have a single units
    vals <- unique(background)[[1]]
    if(length(vals)>1){
      # Set values of 0 to NA
      background[background == 0] <- NA
      vals <- unique(background)[[1]]
      assertthat::assert_that(length(vals)==1,
                              msg = "Values in background layer can only be unique!")
    }

    # Convert raster to dissolved polygons to get a study boundary
    newbg <- terra_to_sf(background)

    # Check crs is set to be sure
    if(is.na(sf::st_crs(newbg))){
      newbg <- sf::st_set_crs(newbg, sf::st_crs(background))
    }

    # Rerun the distribution call with the object
    distribution(newbg, limits, limits_method, mcp_buffer, limits_clip)
  })

#' @rdname distribution
methods::setMethod(
  "distribution",
  methods::signature(background = "sf"),
  function(background, limits = NULL, limits_method = "none", mcp_buffer = 0, limits_clip = FALSE) {
    # Check that arguments are valid
    assertthat::assert_that(!missing(background) || !exists('background'),
                            inherits(limits,'SpatRaster') || inherits(limits, 'sf') ||
                              inherits(limits, 'Spatial') || is.null(limits) || inherits(limits,"stars"),
                            is.character(limits_method),
                            is.numeric(mcp_buffer),
                            is.logical(limits_clip),
                            msg = 'No background file supplied!')
    assertthat::assert_that(
      inherits(background,'sf'),
      unique(sf::st_geometry_type(background)) %in% c('MULTIPOLYGON','POLYGON')
    )

    # Check that provided background has a valid crs
    assertthat::assert_that(!is.na(sf::st_crs(background)),
                            msg = "Please provide a background with valid projection!")

    # Small checks on alternative limit functionalities
    limits_method <- match.arg(limits_method, c("none","zones", "mcp"), several.ok = FALSE)
    assertthat::assert_that(mcp_buffer>=0, msg = "Buffered mcp distance has to be positive!")

    # Messenger
    if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Setup]','green','Creating distribution object...')

    # Convert limits if provided
    if(!is.null(limits)){
      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Deprecated]','yellow','Use add_limits_extrapolation() to add limits!')
      # Set methods in case a layer was supplied
      limits_method <- "zones"
      # If stars, convert to raster
      if(inherits(limits, 'stars')){
        assertthat::assert_that( length(limits)==1,
                                 msg = "More than 1 attribute found for multi-temporal limits!")
        # Convert to SpatRaster and take the first entry (first time slot)
        limits <- stars_to_raster(limits)
        limits <- Reduce('c', limits) # Combine all
        # Assume it was categorical (stars does not recognize categorical factors)
        layer <- terra::as.factor(layer)
      }
      # Convert to polygon if raster
      if(inherits(limits,'SpatRaster')){
        # If there are more one layer, take the first
        assertthat::assert_that(all(terra::is.factor(limits)),
                                msg = 'Provided limit raster needs to be ratified (categorical)!')
        limits <- terra_to_sf(limits)
      }
      assertthat::assert_that(inherits(limits, "sf"))
      # If there is no time dimension, assign a dummy
      if(!utils::hasName(limits,"time")) limits$time <- "No set time"

      # Ensure that limits has the same projection as background
      if(sf::st_crs(limits) != sf::st_crs(background)) limits <- sf::st_transform(limits, background)
      # Ensure that limits is intersecting the background
      if(is.Raster(background)){
        if(suppressMessages(length( sf::st_intersects(limits, terra::as.polygons(background) |> sf::st_as_sf()) )) == 0 ) { limits <- NULL; cli::cli_alert_warning('Provided limits do not intersect the background!') }
      } else {
        if(suppressMessages(length( sf::st_intersects(limits, background |> sf::st_as_sf()) )) == 0 ) { limits <- NULL; cli::cli_alert_warning('Provided limits do not intersect the background!') }
      }

      # Rename geometry just to be sure
      limits <- rename_geometry(limits, "geometry")
      names(limits)[1] <- "limit" # The first entry is the layer name
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

    # Create a dummy biodiversity dataset at initalization
    bio <- BiodiversityDatasetCollection$new()

    # Create BiodiversityDistribution object
    bd <- BiodiversityDistribution$new(
      background = background,
      limits = limits,
      biodiversity = bio
    )
    return(bd)
  })
