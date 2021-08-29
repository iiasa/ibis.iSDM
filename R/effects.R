#' @include utils.R
NULL

#' Plot effects of trained model
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
#' @export
effects.DistributionModel <- function(x, ...) x$effects()
