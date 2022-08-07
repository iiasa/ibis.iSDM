#' @include utils.R bdproto.R bdproto-prior.R
NULL

#' Create a new monotonic prior
#'
#' @description
#' Function to include prior information as monotonic constrain to a
#' extreme gradient descent boosting model [`engine_xgboost`].
#' Monotonic priors enforce directionality in direction of certain variables, however
#' specifying a monotonic constrain does not guarantee that the variable is not regularized out
#' during model fitting.
#'
#' @param variable A [`character`] matched against existing predictors or latent effects.
#' @param hyper A [`character`] object describing the type of constrain. Available options are \code{'increasing'}, \code{'decreasing'}, \code{'convex'}, \code{'concave'}, \code{'none'}.
#' @param ... Variables passed on to prior object.
#' @references
#' * Chen, T., He, T., Benesty, M., Khotilovich, V., Tang, Y., & Cho, H. (2015). Xgboost: extreme gradient boosting. R package version 0.4-2, 1(4), 1-4.
#' @seealso [`Prior-class`] and [`GDBPrior-class`].
#' @family prior
#' @keywords priors
#' @aliases XGBPrior
#' @name XGBPrior
NULL

#' @name XGBPrior
#' @rdname XGBPrior
#' @exportMethod XGBPrior
#' @export
methods::setGeneric(
  "XGBPrior",
  signature = methods::signature("variable", "hyper"),
  function(variable, hyper = 'increasing', ...) standardGeneric("XGBPrior"))

#' @name XGBPrior
#' @rdname XGBPrior
#' @usage \S4method{XGBPrior}{character, character}(variable, hyper)
methods::setMethod(
  "XGBPrior",
  methods::signature(variable = "character", hyper = "character"),
  function(variable, hyper = 'increasing', ... ) {
    assertthat::assert_that(!missing(variable),!missing(hyper),
                            msg = 'Variable or constrain unset.')
    assertthat::assert_that(
      is.character(variable),
      is.character(hyper)
    )
    assertthat::assert_that(length(variable)==1,msg = 'More than one prior variable supplied. Use XGBPriors')
    # Match supplied constrain in case someone has been lazy
    hyper <- match.arg(hyper, c('increasing', 'decreasing','positive', 'negative', 'none'), several.ok = FALSE)

    # Create new prior object
    bdproto(
      'XGBPrior',
      Prior,
      id = new_id(),
      variable = variable,
      value = hyper
    )
  }
)

#' Helper function when multiple variables are supplied
#' @name XGBPriors
#' @description
#' This is a helper function to specify several [XGBPrior] with the same
#' hyper-parameters, but different variables.
#' @rdname XGBPriors
#' @exportMethod XGBPriors
#' @inheritParams XGBPrior
#' @family prior
#' @keywords priors
#' @export
methods::setGeneric(
  "XGBPriors",
  signature = methods::signature("variable"),
  function(variable, hyper = 'increasing', ...) standardGeneric("XGBPriors"))

#' @name XGBPriors
#' @rdname XGBPriors
#' @usage \S4method{XGBPriors}{character, character}(variable, hyper)
methods::setMethod(
  "XGBPriors",
  methods::signature(variable = "character"),
  function(variable, hyper = 'increasing', ... ) {
    assertthat::assert_that(!missing(variable),!missing(hyper),
                            msg = 'Variable or constrain unset.')
    assertthat::assert_that(
      is.character(variable),
      is.character(hyper), length(hyper)==1
    )
    assertthat::assert_that(length(variable)>1, msg = 'Only one prior variable supplied. Use XGBPrior')
    # Match supplied constrain in case someone has been lazy
    hyper <- match.arg(hyper, c('increasing', 'decreasing','positive', 'negative', 'none'), several.ok = FALSE)

    multiple_priors <- list()
    for(k in variable){
      np <- XGBPrior(variable = k,hyper = hyper)
      multiple_priors[[as.character(np$id)]] <- np
    }
    return(multiple_priors)
  }
)
