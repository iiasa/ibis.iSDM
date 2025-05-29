#' @include class-biodiversityscenario.R
NULL

#' Create a new scenario based on trained model parameters
#'
#' @description This function creates a new [BiodiversityScenario-class] object
#' that contains the projections of a model.
#'
#' @param fit A [`BiodiversityDistribution`] object containing a trained model.
#' @param limits A [`SpatRaster`] or [`sf`] object that limits the projection
#' surface when intersected with the prediction data (Default: \code{NULL}).
#' This can for instance be set as an expert-delineated constrain to limit spatial
#' projections.
#' @param reuse_limits A [`logical`] on whether to reuse limits if found in the
#' trained [`BiodiversityDistribution`] object (Default: \code{FALSE}). See also notes!
#' @param copy_model A [`logical`] of whether the model object is to be copied to
#' the scenario object. Note that setting this option to \code{TRUE} can increase
#' the required amount of memory (Default: \code{FALSE}).
#'
#' @note
#' If a limit has been defined already during [train()], for example by adding
#' an extrapolation limit [add_limits_extrapolation()], this zonal layer can be
#' reused for the projections. **Note: This effectively fixes the projections to certain areas.**
#'
#' @examples
#' \dontrun{
#'   scenario(fit, limits = island_area)
#' }
#'
#' @name scenario
NULL

#' @rdname scenario
#' @export
methods::setGeneric("scenario",
                    signature = methods::signature("fit"),
                    function(fit, limits = NULL, reuse_limits = FALSE, copy_model = FALSE) standardGeneric("scenario"))

#' @rdname scenario
methods::setMethod(
  "scenario",
  methods::signature(fit = "ANY"),
  function(fit, limits = NULL, reuse_limits = FALSE, copy_model = FALSE) {
    # Check that arguments are valid
    assertthat::assert_that(!missing(fit) || inherits(fit,'DistributionModel'),
                            inherits(limits, 'SpatRaster') || inherits(limits, 'sf') || inherits(limits, 'Spatial') ||
                              is.null(limits) || inherits(limits, "stars"),
                            is.logical(reuse_limits),
                            is.logical(copy_model),
                            msg = 'No trained model supplied!')

    # Get model object name and id
    if(!copy_model){
      modelobject <- deparse(substitute(fit))
    } else {
      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Scenario]','yellow', "Saving model directly in scenario object!")
      modelobject <- fit$clone(deep = TRUE)
    }
    modelid <- fit$id

    # Convert limits if provided
    if(!is.null(limits)){
      # Convert to polygon if raster
      if(is.Raster(limits)){
        # Remove 0 from ratified raster assuming this is no-data
        limits[limits == 0] <- NA
        limits <- terra_to_sf(limits)
      } else if(inherits(limits, "stars")){
        limits <- stars_to_sf(limits)
      }
      # Ensure that limits has the same projection as background
      if(sf::st_crs(limits) != sf::st_crs(fit$model$background)) limits <- sf::st_transform(limits, fit$model$background)
      # Ensure that limits is intersecting the background
      if(suppressWarnings( suppressMessages(length( sf::st_intersects(limits, fit$model$background)))==0)) {
        limits <- NULL; cli::cli_alert_warning('Provided limits do not intersect the background!')
      }

      # Rename geometry just to be sure
      if(inherits(limits, 'sf')){
        limits <- rename_geometry(limits, "geometry")
      }
    }

    # Also check if limits are to be reused if found
    if(reuse_limits){
      # Check if limits have been found.
      if(fit$has_limits()){
        if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Scenario]','yellow', "Found existing extrapolation limits and will reuse them!")
        settings <- fit$settings
        limits <- settings$get('limits')[['layer']]
      }
    }

    if(is.null(limits)){
      # Convert to waiver if NULL
      limits <- new_waiver()
    }
    # Create BiodiversityScenario object
    sc <- BiodiversityScenario$new()
    sc$modelobject <- modelobject
    sc$modelid <- modelid
    sc$limits <- limits

    return(sc)

  })
