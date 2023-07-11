#' @include utils.R
NULL

#' Obtain partial effects of trained model
#'
#' Create a partial response or effect plot of a trained model
#'
#' @param mod A trained `DistributionModel` object with \code{fit_best} model within.
#' @param x.var A [character] indicating the variable for which a partial effect is to be calculated.
#' @param constant A [numeric] constant to be inserted for all other variables. Default calculates a mean per variable.
#' @param variable_length [numeric] The interpolation depth (nr. of points) to be used (Default: \code{100}).
#' @param values [numeric] Directly specified values to compute partial effects for. If this parameter is set to anything other
#' than \code{NULL}, the parameter \code{"variable_length"} is ignored (Default: \code{NULL}).
#' @param plot A [`logical`] indication of whether the result is to be plotted?
#' @param type A specified type, either \code{'response'} or \code{'predictor'}. Can be missing.
#' @param ... Other engine specific parameters.
#' @seealso [partial]
#' @details By default the mean is calculated across all parameters that are not \code{x.var}.
#' Instead a *constant* can be set (for instance \code{0}) to be applied to the output.
#' @return A [data.frame] with the created partial response.
#' @aliases partial
#' @examples
#' \dontrun{
#'  # Do a partial calculation of a trained model
#'  partial(fit, x.var = "Forest.cover", plot = TRUE)
#' }
#' @keywords partial
#' @export
#' @name partial
methods::setGeneric(
  "partial",
  signature = methods::signature("mod","x.var"),
  function(mod, x.var, constant = NULL, variable_length = 100, values = NULL, plot = FALSE, type = "response", ...) standardGeneric("partial"))

#' @name partial
#' @rdname partial
#' @usage \S4method{partial}{ANY,character,ANY,numeric,ANY,logical,character}(mod,x.var,constant,variable_length,values,plot,type,...)
methods::setMethod(
  "partial",
  methods::signature(mod = "ANY", x.var = "character"),
  function(mod, x.var, constant = NULL, variable_length = 100, values = NULL, plot = FALSE, type = "response",...) {
    assertthat::assert_that(!missing(x.var),msg = 'Specify a variable name in the model!')
    assertthat::assert_that(inherits(mod, "DistributionModel"),
                            is.character(x.var),
                            is.null(constant) || is.numeric(constant),
                            is.numeric(variable_length),
                            is.numeric(values) || is.null(values),
                            is.logical(plot)
    )
    # Work around to call partial response directly
    if(inherits(mod,'DistributionModel')){
      partial.DistributionModel(mod, x.var, constant, variable_length, values, plot, type, ...)
    } else {
      stop('Partial response calculation not supported!')
    }
  }
)

#' @rdname partial
#' @method partial DistributionModel
#' @keywords partial
#' @export
partial.DistributionModel <- function(mod, ...) mod$partial(...)

#' Obtain spatial partial effects of trained model
#'
#' @description
#' Similar as [partial] this function calculates a partial response of a trained model for a given variable.
#' Differently from [partial] in space. However the result is a [`SpatRaster`] showing the spatial magnitude
#' of the partial response.
#' @param mod A [`DistributionModel-class`] object with trained model.
#' @param x.var A [character] indicating the variable for which a partial effect is to be calculated.
#' @param constant A [numeric] constant to be inserted for all other variables. Default calculates the [mean] per variable.
#' @param plot A [logical] indication of whether the result is to be plotted?
#' @param ... Other engine specific parameters.
#' @seealso [partial]
#' @details By default the [mean] is calculated across all parameters that are not \code{x.var}.
#' Instead a *constant* can be set (for instance \code{0}) to be applied to the output.
#' @returns A [SpatRaster] containing the mapped partial response of the variable.
#' @aliases spartial
#' @examples
#' \dontrun{
#'  # Create and visualize the spartial effect
#'  spartial(fit, x.var = "Forest.cover", plot = TRUE)
#' }
#' @keywords partial
#' @export
#' @name spartial
methods::setGeneric(
  "spartial",
  signature = methods::signature("mod","x.var"),
  function(mod, x.var, constant = NULL, plot = FALSE, ...) standardGeneric("spartial"))

#' @name spartial
#' @rdname spartial
#' @usage \S4method{spartial}{ANY,character,ANY,logical}(mod,x.var,constant,plot,...)
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
#' @keywords partial
#' @export
spartial.DistributionModel <- function(mod, ...) mod$spartial(...)
