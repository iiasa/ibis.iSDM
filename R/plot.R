#' @include utils.R
NULL

#' Plot wrappers
#'
#' Plots information from a given object where a plotting object is available
#'
#' @param x Any object belonging to [DistributionModel], [BiodiversityDatasetCollection], [PredictorDataset] or [BiodiversityScenario].
#' @param ... Further arguments passed on to \code{x$plot}
#'
#' @details
#' The plotted outputs vary depending on what object is being plotted.
#' For example for fitted [DistributionModel] the output is usually the fitted spatial prediction (Default: 'mean')
#' @return Graphical output
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
