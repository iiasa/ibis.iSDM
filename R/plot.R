#' @include utils.R
NULL

#' Plot wrappers
#'
#' Plots information from a given object where a plotting object is available
#'
#' @param x Any object.
#' @param ... not used.
#'
#' @return None.
#' @name plot
NULL

#' @rdname plot
#' @method plot DistributionModel
#'
#' @export
#'
plot.DistributionModel <- function(x, ...) x$plot(...)
