#' Plot effects of trained model
#'
#' @description This functions is handy wrapper that calls the default plotting
#' functions for the model of a specific engine. Equivalent to calling \code{effects}
#' of a fitted [distribution] function.
#'
#' @param object Any fitted [distribution] object.
#' @param ... Not used.
#'
#' @note For some models, where default coefficients plots are not available,
#' this function will attempt to generate [partial] dependency plots instead.
#'
#' @return None.
#'
#' @keywords partial
#'
#' @examples
#' \dontrun{
#' # Where mod is an estimated distribution model
#' mod$effects()
#' }
#'
#' @name effects
NULL

#' @rdname effects
#' @export
effects.DistributionModel <- function(object, ...) object$effects()
