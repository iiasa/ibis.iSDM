#' @include utils.R
NULL

#' Plot effects of trained model
#'
#' Create a partial response or effect plot of a trained model
#'
#' @description This functions is handy wrapper that calls the default plotting
#' functions for the model of a specific [engine]. Equivalent to
#' calling \code{effects} of a fitted [distribution] function.
#' @param x Any fitted [distribution] object.
#' @param ... not used.
#'
#' @examples
#' \dontrun{
#' # Where mod is an estimated distribution model
#' mod$effects()
#' }
#' @return None.
#' @keywords misc
#' @references
#' @name effects
NULL

#' @rdname effects
#' @method effects DistributionModel
#' @keywords misc
#' @export
effects.DistributionModel <- function(x, ...) x$effects()
