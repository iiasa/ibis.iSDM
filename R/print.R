#' @include utils.R
NULL

#' Print
#'
#' Display information about an object.
#'
#' @param x Any object.
#' @param ... not used.
#'
#' @return None.
#' @seealso [base::print()].
#'
#' @name print
#' @aliases print,Id-method print,tbl_df-method
#'
#' @examples
#' a <- 1:4
#' print(a)
NULL

#' @rdname print
#' @method print distribution
#'
#' @export
#'
print.distribution <- function(x, ...) x$print()

#' @rdname print
#' @method print BiodiversityDistribution
#'
#' @export
#'
print.BiodiversityDistribution <- function(x, ...) x$print()

#' @rdname print
#' @method print BiodiversityDatasetCollection
#'
#' @export
#'
print.BiodiversityDatasetCollection <- function(x, ...) x$print()

#' @rdname print
#' @method print BiodiversityDataset
#'
#' @export
#'
print.BiodiversityDataset <- function(x, ...) x$print()

#' @rdname print
#' @method print PredictorDataset
#'
#' @export
#'
print.PredictorDataset <- function(x, ...) x$print()

#' @rdname print
#' @method print DistributionModel
#'
#' @export
#'
print.DistributionModel <- function(x, ...) x$print()

#' @rdname print
#' @method print BiodiversityScenario
#'
#' @export
#'
print.BiodiversityScenario <- function(x, ...) x$print()

#' @rdname print
#' @method print Prior
#'
#' @export
#'
print.Prior <- function(x, ...) x$print()

#' @rdname print
#' @method print PriorList
#'
#' @export
#'
print.PriorList <- function(x, ...) x$print()

#' @rdname print
#' @method print Engine
#'
#' @export
#'
print.Engine <- function(x, ...) x$print()

#' @rdname print
#' @method print Settings
#'
#' @export
#'
print.Settings <- function(x, ...) x$print()


#' @rdname print
#' @method print Log
#'
#' @export
#'
print.Log <- function(x, ...) x$print()

#' @rdname print
#'
#' @method print Id
#'
#' @export
#'
print.Id <- function(x, ...) message("id: ", x)

#' @name print
#'
#' @rdname print
#'
#' @usage \S4method{print}{Id}(x)
#'
methods::setMethod("print", "Id", function(x, ...) print.Id(x))

#' @name print
#'
#' @rdname print
#'
#' @usage \S4method{print}{tbl_df}(x)
#'
methods::setMethod("print", "tbl_df", function(x, ...) base::print(x, ...))
