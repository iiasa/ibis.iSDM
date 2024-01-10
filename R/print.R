#' Print
#'
#' Display information about any object created through the \pkg{ibis.iSDM}
#' R-package.
#'
#' @param x Any object created through the package.
#' @param ... not used.
#'
#' @return Object specific.
#'
#' @seealso [base::print()].
#' @keywords misc
#'
#' @examples
#' \dontrun{
#' # Where mod is fitted object
#' mod
#' print(mod)
#' }
#'
#' @name print
NULL

#' @rdname print
#' @method print distribution
#' @export
print.distribution <- function(x, ...) x$print()

#' @rdname print
#' @method print BiodiversityDistribution
#' @export
print.BiodiversityDistribution <- function(x, ...) x$print()

#' @rdname print
#' @method print BiodiversityDatasetCollection
#' @export
print.BiodiversityDatasetCollection <- function(x, ...) x$print()

#' @rdname print
#' @method print BiodiversityDataset
#' @export
print.BiodiversityDataset <- function(x, ...) x$print()

#' @rdname print
#' @method print PredictorDataset
#' @export
print.PredictorDataset <- function(x, ...) x$print()

#' @rdname print
#' @method print DistributionModel
#' @export
print.DistributionModel <- function(x, ...) x$print()

#' @rdname print
#' @method print BiodiversityScenario
#' @export
print.BiodiversityScenario <- function(x, ...) x$print()

#' @rdname print
#' @method print Prior
#' @export
print.Prior <- function(x, ...) x$print()

#' @rdname print
#' @method print PriorList
#' @export
print.PriorList <- function(x, ...) x$print()

#' @rdname print
#' @method print Engine
#' @export
print.Engine <- function(x, ...) x$print()

#' @rdname print
#' @method print Settings
#' @export
print.Settings <- function(x, ...) x$print()

#' @rdname print
#' @method print Log
#' @export
print.Log <- function(x, ...) x$print()

#' @rdname print
#' @method print Id
#' @export
print.Id <- function(x, ...) message("id: ", x)

#' @rdname print
methods::setMethod("print", "Id", function(x, ...) print.Id(x))

#' @rdname print
#' @keywords misc
methods::setMethod("print", "tbl_df", function(x, ...) base::print(x, ...))
