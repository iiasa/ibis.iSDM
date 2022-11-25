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

#' @rdname summary
#' @method summary PriorList
#' @keywords summary
#' @export
summary.PriorList <- function(x) x$summary()
