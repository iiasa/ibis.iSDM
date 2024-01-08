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
#' @seealso [base::print()].
#'
#' @name print
#' @aliases print, Id-method print, tbl_df-method
#' @keywords misc
#' @examples
#' \dontrun{
#' # Where mod is fitted object
#' mod
#' print(mod)
#' }
NULL

#' @rdname print
#' @method print distribution
#' @keywords misc
#' @export
print.distribution <- function(x, ...) x$print()

#' @rdname print
#' @method print BiodiversityDistribution
#' @keywords misc
#' @export
print.BiodiversityDistribution <- function(x, ...) x$print()

#' @rdname print
#' @method print BiodiversityDatasetCollection
#' @keywords misc
#' @export
print.BiodiversityDatasetCollection <- function(x, ...) x$print()

#' @rdname print
#' @method print BiodiversityDataset
#' @keywords misc
#' @export
print.BiodiversityDataset <- function(x, ...) x$print()

#' @rdname print
#' @method print PredictorDataset
#' @keywords misc
#' @export
print.PredictorDataset <- function(x, ...) x$print()

#' @rdname print
#' @method print DistributionModel
#' @keywords misc
#' @export
print.DistributionModel <- function(x, ...) x$print()

#' @rdname print
#' @method print BiodiversityScenario
#' @keywords misc
#' @export
print.BiodiversityScenario <- function(x, ...) x$print()

#' @rdname print
#' @method print Prior
#' @keywords misc
#' @export
print.Prior <- function(x, ...) x$print()

#' @rdname print
#' @method print PriorList
#' @keywords misc
#' @export
print.PriorList <- function(x, ...) x$print()

#' @rdname print
#' @method print Engine
#' @keywords misc
#' @export
print.Engine <- function(x, ...) x$print()

#' @rdname print
#' @method print Settings
#' @keywords misc
#' @export
print.Settings <- function(x, ...) x$print()

#' @rdname print
#' @method print Log
#' @keywords misc
#' @export
print.Log <- function(x, ...) x$print()

#' @rdname print
#' @method print Id
#' @keywords misc
#' @export
print.Id <- function(x, ...) message("id: ", x)

#' @name print
#' @rdname print
#' @keywords misc
methods::setMethod("print", "Id", function(x, ...) print.Id(x))

#' @name print
#' @rdname print
#' @keywords misc
methods::setMethod("print", "tbl_df", function(x, ...) base::print(x, ...))
