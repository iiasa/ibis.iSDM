#' @include utils.R
NULL

#' @export
if (!methods::isClass("Id")) methods::setOldClass("Id")
NULL

#' Identifier
#'
#' Generate a new unique identifier.
#'
#' @details Identifiers are made using the [uuid::UUIDgenerate()].
#'
#' @return `Id` object.
#'
#' @seealso [uuid::UUIDgenerate()].
#'
#' @examples
#' # create new id
#' i <- new_id()
#'
#' # print id
#' print(i)
#'
#' # convert to character
#' as.character(i)
#'
#' # check if it is an Id object
#' is.Id(i)
#'
#' @aliases Id
#'
#' @export
new_id <- function() {
  x <- uuid::UUIDgenerate()
  class(x) <- c("Id", class(x))
  x
}

#' As Id
#' @rdname as
#' @export
as.Id <- function(x, ...) UseMethod("as.Id")

#' As Id character
#' @rdname as
#' @export
as.Id.character <- function(x, ...) {
  class(x) <- c("Id", class(x))
  x
}

#' Is Id
#' @rdname is
#' @export
is.Id <- function(x) inherits(x, "Id")

#' Is waiver
#' @rdname is
#' @export
is.Waiver <- function(x) inherits(x, "Waiver")

#' Check whether a formula is valid
#'
#' @param x A [`character`] object
#' @return Boolean evaluation
#' @rdname is
#' @export

is.formula <- function(x) inherits(x, "formula")
