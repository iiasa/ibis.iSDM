#' @include utils.R bdproto-biodiversityscenario.R
NULL

#' Create a new scenario based on trained model parameters
#'
#' @description This function creates a new [BiodiversityScenario-class] object
#' that contains the projections of a model.
#' @note
#' If a limit has been defined already during [train()], for example by adding
#' an extrapolation limit [add_control_extrapolation()], this zonal layer can be
#' reused for the projections. **Note: This effectively fixes the projections to certain areas.**
#'
#' @param fit A [`BiodiversityDistribution`] object containing a trained model.
#' @param limits A [`SpatRaster`] or [`sf`] object that limits the projection
#'   surface when intersected with the prediction data (Default: \code{NULL}).
#'   This can for instance be set as an expert-delineated constrain to limit
#'   spatial projections.
#' @param reuse_limits A [`logical`] on whether to reuse limits if found in the
#' trained [`BiodiversityDistribution`] object (Default: \code{FALSE}). See also notes!
#'
#' @param copy_model A [`logical`] of whether the model object is to be copied
#'   to the scenario object. Note that setting this option to \code{TRUE} can
#'   increase the required amount of memory (Default: \code{FALSE}).
#' @aliases scenario
#' @name scenario
#'
#' @examples
#' \dontrun{
#'   scenario(fit, limits = island_area)
#' }
#' @export
methods::setGeneric("scenario",
                    signature = methods::signature("fit"),
                    function(fit, limits = NULL, reuse_limits = FALSE, copy_model = FALSE) standardGeneric("scenario"))

#' @name scenario
#' @rdname scenario
methods::setMethod(
  "scenario",
  methods::signature(fit = "ANY"),
  function(fit, limits = NULL, reuse_limits = FALSE, copy_model = FALSE) {
    # Check that arguments are valid
    assertthat::assert_that(!missing(fit) || inherits(fit,'DistributionModel'),
                            inherits(limits, 'SpatRaster') || inherits(limits, 'sf') || inherits(limits, 'Spatial') || is.null(limits),
                            is.logical(reuse_limits),
                            is.logical(copy_model),
                            msg = 'No trained model supplied!')

    # Get model object name and id
    if(!copy_model){
      modelobject <- deparse(substitute(fit))
    } else {
      if(getOption('ibis.setupmessages')) myLog('[Scenario]','yellow', "Saving model directly in scenario object!")
      modelobject <- fit
    }
    modelid <- fit$id

    # Convert limits if provided
    if(!is.null(limits)){
      # Convert to polygon if raster
      if(inherits(limits,'SpatRaster')){
        if(terra::is.factor(limits)) stop('Provided limit raster needs to be ratified (categorical)!')
        # Remove 0 from ratified raster assuming this is no-data
        limits[limits == 0] <- NA
        limits <- sf::st_as_sf( terra::as.polygons(limits, trunc = TRUE, dissolve = TRUE) )
      }
      # Ensure that limits has the same projection as background
      if(sf::st_crs(limits) != sf::st_crs(fit$model$background)) limits <- sf::st_transform(limits, fit$model$background)
      # Ensure that limits is intersecting the background
      if(suppressMessages(length( sf::st_intersects(limits, fit$model$background)))==0) { limits <- NULL; warning('Provided limits do not intersect the background!') }

      # Get fir column and rename
      limits <- limits[,1]; names(limits) <- c('limit','geometry')
    }

    # Also check if limits are to be reused if found
    if(reuse_limits){
      # Check if limits have been found.
      if(fit$has_limits()){
        if(getOption('ibis.setupmessages')) myLog('[Scenario]','yellow', "Found existing extrapolation limits and will reuse them!")
        settings <- fit$settings
        limits <- settings$get('limits')[['layer']]
      }
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
