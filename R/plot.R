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
#' @keywords misc
#' @name plot
NULL

#' @rdname plot
#' @method plot DistributionModel
#' @keywords misc
#' @export
plot.DistributionModel <- function(x, ...) x$plot(...)

#' @rdname plot
#' @method plot BiodiversityDatasetCollection
#' @keywords misc
#' @export
plot.BiodiversityDatasetCollection <- function(x, ...) x$plot(...)

#' @rdname plot
#' @method plot PredictorDataset
#' @keywords misc
#' @export
plot.PredictorDataset <- function(x, ...) x$plot(...)

#' @rdname plot
#' @method plot Engine
#' @keywords misc
#' @export
plot.Engine <- function(x,...) x$plot(...)

#' @rdname plot
#' @method plot BiodiversityScenario
#' @keywords misc
#' @export
plot.BiodiversityScenario <- function(x,...) x$plot(...)
