if (!methods::isClass("Waiver")) methods::setOldClass("Waiver")

#' Waiver
#'
#' @description Create a `waiver` object.
#'
#' @details This object is used to represent that the user has not manually
#' specified a setting, and so defaults should be used. By explicitly using a
#' `new_waiver()`, this means that `NULL` objects can be a valid setting. The
#' use of a "waiver" object was inspired by the `ggplot2` and `prioritizr` package.
#'
#' @return Object of class `Waiver`.
#'
#' @keywords misc
#'
#' @examples
#' # create new waiver object
#' w <- new_waiver()
#'
#' # print object
#' print(w)
#'
#' # is it a waiver object?
#' is.Waiver(w)
#'
#' @export
new_waiver <- function() structure(list(), class = "Waiver")
