#' @include utils.R
NULL

#' Summarises a trained model or predictor object
#'
#' @param x Any prepared object.
#' @param ... not used.
#'
#' @return None.
#' @seealso [base::summary()].
#' @keywords misc
#' @name summary
NULL

#' @rdname summary
#' @method summary distribution
#' @keywords misc
#' @export
summary.distribution <- function(x) x$summary()

#' @rdname summary
#' @method summary DistributionModel
#' @keywords misc
#' @export
summary.DistributionModel <- function(x) x$summary()

#' @rdname summary
#' @method summary PredictorDataset
#' @keywords misc
#' @export
summary.PredictorDataset <- function(x) x$summary()
