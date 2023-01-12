#' @include utils.R
NULL

#' Summarises a trained model or predictor object
#'
#' @description
#' This helper function summarizes a given object, including
#' [DistributionModel] and [PredictorDataset].
#' @param x Any prepared object.
#' @param ... not used.
#'
#' @return None.
#' @seealso [base::summary()].
#' @keywords summary
#' @name summary
NULL

#' @rdname summary
#' @method summary distribution
#' @keywords summary
#' @export
summary.distribution <- function(x) x$summary()

#' @rdname summary
#' @method summary DistributionModel
#' @keywords summary
#' @export
summary.DistributionModel <- function(x) x$summary()

#' @rdname summary
#' @method summary PredictorDataset
#' @keywords summary
#' @export
summary.PredictorDataset <- function(x) x$summary()

#' @rdname summary
#' @method summary BiodiversityScenario
#' @keywords summary
#' @export
summary.BiodiversityScenario <- function(x) x$summary()

#' Obtains the coefficients of a trained model
#'
#' @description
#' Similar as [`summary`], this helper function obtains the coefficients from
#' a given [DistributionModel] object.
#' **Note:**
#' For models trained with machine-learning approaches (e.g. [`engine_bart`] etc) this function
#' will return variable importance estimates rather than linear coefficients.
#' Similar can be said for trained non-linear models.
#' @param object Any prepared object.
#' @param ... not used.
#' @seealso [stats::coef()].
#' @keywords coef
#' @name coef
NULL

#' @rdname coef
#' @method coef DistributionModel
#' @keywords coef
#' @export
coef.DistributionModel <- function(x) x$get_coefficients()
