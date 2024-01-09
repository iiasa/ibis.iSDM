#' @include utils.R
NULL

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
#' @export
print.distribution <- function(x, ...) x$print()

#' @rdname print
#' @export
print.BiodiversityDistribution <- function(x, ...) x$print()

#' @rdname print
#' @export
print.BiodiversityDatasetCollection <- function(x, ...) x$print()

#' @rdname print
#' @export
print.BiodiversityDataset <- function(x, ...) x$print()

#' @rdname print
#' @export
print.PredictorDataset <- function(x, ...) x$print()

#' @rdname print
#' @export
print.DistributionModel <- function(x, ...) x$print()

#' @rdname print
#' @export
print.BiodiversityScenario <- function(x, ...) x$print()

#' @rdname print
#' @export
print.Prior <- function(x, ...) x$print()

#' @rdname print
#' @export
print.PriorList <- function(x, ...) x$print()

#' @rdname print
#' @export
print.Engine <- function(x, ...) x$print()

#' @rdname print
#' @export
print.Settings <- function(x, ...) x$print()

#' @rdname print
#' @export
print.Log <- function(x, ...) x$print()

#' @rdname print
#' @export
print.Id <- function(x, ...) message("id: ", x)

#' @rdname print
methods::setMethod("print", "Id", function(x, ...) print.Id(x))

#' @rdname print
#' @keywords misc
methods::setMethod("print", "tbl_df", function(x, ...) base::print(x, ...))
