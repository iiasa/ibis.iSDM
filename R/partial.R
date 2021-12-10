#' @include utils.R
NULL

#' Obtain partial effects of trained model
#'
#' Create a partial response or effect plot of a trained model
#'
#' @param mod A trained `DistributionModel` model
#' @param x.var A [`character`] indicating the variable for which a partial effect is to be calculated
#' @param constant A [`numeric`] constant to be inserted for all other variables. Default calculates a mean per variable
#' @param variable_length The interpolation depth (nr. of points) to be used (Default: 100)
#' @param plot A [`logical`] indication of whether the result is to be plotted ?
#' @param ... Other engine specific parameters
#' @return None
#' @keywords misc
#' @name partial
methods::setGeneric(
  "partial",
  signature = methods::signature("mod","x.var"),
  function(mod, x.var, constant = NULL, variable_length = 100, plot = FALSE, ...) standardGeneric("partial"))

#' @name partial
#' @rdname partial
#' @usage \S4method{partial}{ANY,character}(mod, x.var)
methods::setMethod(
  "partial",
  methods::signature(mod = "ANY", x.var = "character"),
  function(mod, x.var, constant = NULL, variable_length = 100, plot = FALSE, ...) {
    assertthat::assert_that(!missing(x.var),msg = 'Specify a variable name in the model!')
    assertthat::assert_that(inherits(mod, "DistributionModel"),
                            is.character(x.var),
                            is.null(constant) || is.numeric(constant),
                            is.numeric(variable_length),
                            is.logical(plot)
    )
    # Work around to call partial response directly
    if(inherits(mod,'DistributionModel')){
      partial.DistributionModel(mod, x.var, constant, variable_length, plot, ...)
    } else {
      stop('Partial response calculation not supported!')
    }
  }
)

#' @rdname partial
#' @method partial DistributionModel
#' @keywords misc
#' @export
partial.DistributionModel <- function(x, ...) x$partial(...)

#' Obtain spatial partial effects of trained model
#'
#' Similar as \code{partial} calculates partial response of a trained model for a given variable
#' However the result is a raster showing the spatial magnitude of the partial response
#'
#' @param mod A trained `DistributionModel` model
#' @param x.var A [`character`] indicating the variable for which a partial effect is to be calculated
#' @param constant A [`numeric`] constant to be inserted for all other variables. Default calculates a mean per variable
#' @param plot A [`logical`] indication of whether the result is to be plotted ?
#' @param ... Other engine specific parameters
#' @keywords misc
#' @name spartial
methods::setGeneric(
  "spartial",
  signature = methods::signature("mod","x.var"),
  function(mod, x.var, constant = NULL, variable_length = 100, plot = FALSE, ...) standardGeneric("spartial"))

#' @name spartial
#' @rdname spartial
#' @usage \S4method{spartial}{ANY,character}(mod, x.var)
methods::setMethod(
  "spartial",
  methods::signature(mod = "ANY", x.var = "character"),
  function(mod, x.var, constant = NULL, plot = FALSE, ...) {
    assertthat::assert_that(!missing(x.var),msg = 'Specify a variable name in the model!')
    assertthat::assert_that(inherits(mod, "DistributionModel"),
                            is.character(x.var),
                            is.null(constant) || is.numeric(constant),
                            is.logical(plot)
    )
    # Work around to call partial response directly
    if(inherits(mod,'DistributionModel')){
      spartial.DistributionModel(mod, x.var, constant, plot, ...)
    } else {
      stop('Spatial partial response calculation not supported!')
    }
  }
)

#' @rdname spartial
#' @method spartial DistributionModel
#' @keywords misc
#' @export
spartial.DistributionModel <- function(x, ...) x$spartial(...)

