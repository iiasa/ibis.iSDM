#' @include utils.R
NULL

#' Pipe operator
#'
#' This package uses the pipe operator (`\%>\%`) to express nested code
#' as a series of imperative procedures.
#'
#' @param lhs,rhs An object and a function.
#' @seealso [magrittr::%>%()], [tee()].
#' @return An object.
#' @examples
#' # set seed for reproducibility
#' set.seed(500)
#'
#' # generate 100 random numbers and calculate the mean
#' mean(runif(100))
#'
#' # reset the seed
#' set.seed(500)
#'
#' # repeat the previous procedure but use the pipe operator instead of nesting
#' # function calls inside each other.
#' runif(100) %>% mean()
#'
#' @name %>%
#' @rdname pipe
#' @aliases pipe
#' @importFrom magrittr %>%
#' @export
NULL
