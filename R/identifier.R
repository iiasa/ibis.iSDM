if (!methods::isClass("Id")) methods::setOldClass("Id")

#' Identifier
#'
#' Generate a new unique identifier.
#'
#' @details Identifiers are made using the [uuid::UUIDgenerate()].
#'
#' @return \code{"Id"} object.
#'
#' @seealso [uuid::UUIDgenerate()].
#' @keywords misc
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
#' @export
new_id <- function() {
  x <- uuid::UUIDgenerate()
  class(x) <- c("Id", class(x))
  x
}

#' As Id
#'
#' @param x A [`character`] to be converted as id.
#' @param ... Other arguements
#'
#' @keywords misc
#'
#' @name as.Id
NULL

#' @rdname as.Id
#' @export
as.Id <- function(x, ...) UseMethod("as.Id")

#' @rdname as.Id
#' @export
as.Id.character <- function(x, ...) {
  class(x) <- c("Id", class(x))
  x
}

#' Check whether a provided object is truly of a specific type
#'
#' @param x A provided Id object
#'
#' @return Boolean evaluation with [logical] output.
#'
#' @keywords misc
#'
#' @export
is.Id <- function(x) inherits(x, "Id")

#' Is the provided object of type waiver?
#'
#' @param x A provided \code{Waiver} object
#'
#' @return Boolean evaluation with [logical] output.
#'
#' @keywords misc
#'
#' @export
is.Waiver <- function(x) inherits(x, "Waiver")

#' Check whether a formula is valid
#'
#' @param x A [`character`] object
#'
#' @return Boolean evaluation with [logical] output.
#'
#' @keywords misc
#'
#' @export
is.formula <- function(x) inherits(x, "formula")

#' Tests if an input is a SpatRaster object.
#'
#' @param x an R Object.
#'
#' @return Boolean evaluation with [logical] output.
#'
#' @keywords misc
#'
#' @export
is.Raster <- function(x)
{
  return((class(x)[1]=="SpatRaster" || class(x)[1]=="SpatRasterDataset" || class(x)[1]=="SpatRasterCollection"))
}

#' Tests if an input is a stars object.
#'
#' @param x an R Object.
#'
#' @return Boolean evaluation with [logical] output.
#'
#' @keywords misc
#'
#' @export
is.stars <- function(x)
{
  return(inherits(x, "stars") || is.list(x))
}
