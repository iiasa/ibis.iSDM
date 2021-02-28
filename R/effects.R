#' @include utils.R
NULL

#' effects
#'
#' Create a partial response or effect plot of a trained model
#'
#' @param x Any prepared object.
#' @param ... not used.
#'
#' @return None.
#' @name effects
NULL

#' @rdname effects
#' @method effects DistributionModel
#' @noRd
effects.DistributionModel <- function(x, ...) x$effects()
