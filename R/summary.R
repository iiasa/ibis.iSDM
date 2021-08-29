#' @include utils.R
NULL

#' Summarises a trained model or predictor object
#'
#' @param x Any prepared object.
#' @param ... not used.
#'
#' @return None.
#' @seealso [base::summary()].
#'
#' @name summary
NULL

#' @rdname summary
#' @method summary distribution
#' @export
summary.distribution <- function(x) x$summary()

#' @rdname summary
#' @method summary DistributionModel
#' @export
summary.DistributionModel <- function(x) x$summary()

#' @rdname summary
#' @method summary PredictorDataset
#' @export
summary.PredictorDataset <- function(x) x$summary()
