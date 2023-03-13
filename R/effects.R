#' @include utils.R
NULL

#' Plot effects of trained model
#'
#' @description This functions is handy wrapper that calls the default plotting
#' functions for the model of a specific [engine]. Equivalent to
#' calling \code{effects} of a fitted [distribution] function.
#' @note
#' For some models, where default coefficients plots are not available,
#' this function will attempt to generate [partial] dependency plots instead.
#' @param x Any fitted [distribution] object.
#'
#' @examples
#' \dontrun{
#' # Where mod is an estimated distribution model
#' mod$effects()
#' }
#' @return None.
#' @keywords partial
#' @importFrom stats effects
#' @name effects
NULL

#' @rdname effects
#' @method effects DistributionModel
#' @keywords partial
#' @export
effects.DistributionModel <- function(x) x$effects()
