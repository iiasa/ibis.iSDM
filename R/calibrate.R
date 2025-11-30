#' Calibrate a fitted engine again given a set of parameters and data
#'
#' @description
#' This function (re-)calibrates a model in an engine to find the one that best
#' matches the data. Data in this context can be either a parameter range
#' or a new of observational records against which the model should be calibrated
#' against. See details for more information.
#'
#' @note
#' This function largely conducts semi-automatic calibration and does not allow
#' manual fine-tuning of parameters. For more advanced tuning, it is recommended
#' to extract the model from the engine and conduct manual calibration.
#'
#' @details
#' The function `calibrate()` is used to re-calibrate a fitted model in an engine.
#' This is often also known as "tuning" a model. The function takes a fitted model object
#' and a set of parameters, and uses the data to find the best set of parameters.
#'
#' @param mod A [`BiodiversityDistribution`] object.
#' @param point A [`sf`] object containing observational data used for model training.
#' @param field_occurrence A [`character`] location of biodiversity point records (Default: \code{'observed'}).
#' @param project A [`logical`] statement on whether to make a new projection (Default: \code{FALSE}).
#' @param verbose A [`logical`] on whether be chatty. Default option from ibis-wide options.
#' @param ... Any other parameters not directly referenced.
#' @returns Add a new model to a [`BiodiversityDistribution`] object.
#'
#' @references
#' * Sillero, N., Arenas-Castro, S., Enriquez‐Urzelai, U., Vale, C. G., Sousa-Guedes, D., Martínez-Freiría, F., ... & Barbosa, A. M. (2021). Want to model a species niche? A step-by-step guideline on correlative ecological niche modelling. Ecological Modelling, 456, 109671.
#' @keywords tuning
#'
#' @examples
#' \dontrun{
#'  mod <- distribution(background) |>
#'  add_predictors(mypredictors) |>
#'  add_biodiversity_poipa(species) |>
#'  engine_glm()
#'
#'  # Now calibrate again
#'  mod_better <- calibrate(mod)
#' }
#'
#' @name calibrate
NULL

#' @rdname calibrate
#' @export
methods::setGeneric(
  "calibrate",
  signature = methods::signature("mod"),
  function(mod, point = NULL, field_occurrence = "observed",
           project = FALSE, verbose = getOption('ibis.setupmessages', default = TRUE), ...) standardGeneric("calibrate"))

#' @rdname calibrate
methods::setMethod(
  "calibrate",
  methods::signature(mod = "DistributionModel"),
  function(mod, point = NULL, field_occurrence = "observed",
           project = FALSE, verbose = getOption('ibis.setupmessages', default = TRUE), ...) {
    assertthat::assert_that(inherits(mod, "DistributionModel"),
                            is.null(point) || nrow(point)>0,
                            is.character(field_occurrence),
                            is.logical(project),
                            is.logical(verbose)
    )

    # Check that the model is fitted
    assertthat::assert_that(
      mod$has_converged(),
      !is.null(mod$get_data('fit_best')),
      msg = "The model must be fitted before calibration."
    )

    # If data is provided, check
    if(!is.null(point)) {
      assertthat::assert_that(
        inherits(point, "data.frame") || inherits(point, 'sf'),
        nrow(point) > 0,
        msg = "The testing point data must be a non-empty data.frame or sf object."
      )
      assertthat::assert_that(
        all(field_occurrence %in% colnames(point)),
        msg = paste0("The testing point data must contain the column: ", field_occurrence)
      )

      cli::cli_alert_warning("Calibration with point data not yet implemented.")
      point <- NULL
    }

    # --- #
    # Start the calibration
    if (verbose) cli::cli_alert_info("Starting calibration...")

    # Calibrate a model
    modc <- mod$clone(deep = TRUE)
    fit <- try({
      modc$calibrate(newdata = point)
    }, silent = TRUE)
    if(inherits(fit, "try-error")){
      cli::cli_abort("Calibration failed: {fit}. Defaulting to original model.")
      return(mod)
    }
    if (verbose) cli::cli_alert_success("Calibration completed.")

    # Replace the fitted model
    modc$set_data('fit_best', fit)

    # Recreating any projections if found
    if( modc$has_prediction() ){
      # Remake a projection
      if (verbose) cli::cli_alert_info("Recreating model projection...")
      ras <- modc$project(newdata = modc$model$predictors)
      modc$set_data('prediction', ras)
    }

    return(modc)
    # --- #
  }
)
