#' @include utils.R bdproto.R bdproto-prior.R
NULL

#' Create a new monotonic prior
#'
#' Function to include prior information as monotonic constrain to a gradient descent boosting model
#'
#' @param variable A [`character`] matched against existing predictors or latent effects.
#' @param hyper A [`character`] object describing the type of constrain. Available options are \code{'increasing'},
#'  \code{'decreasing'}, \code{'convex'}, \code{'concave'} or \code{'none'}.
#' @param ... Variables passed on to prior object.
#' @details TBD
#' @references
#' * Hofner, B., Müller, J., & Hothorn, T. (2011). Monotonicity‐constrained species distribution models. Ecology, 92(10), 1895-1901.
#' @seealso [`Prior-class`].
#' @keywords priors
#' @family prior
#' @aliases GDBPrior
#' @name GDBPrior
NULL

#' @name GDBPrior
#' @rdname GDBPrior
#' @exportMethod GDBPrior
#' @export
methods::setGeneric(
  "GDBPrior",
  signature = methods::signature("variable"),
  function(variable, hyper = 'increasing', ...) standardGeneric("GDBPrior"))

#' @name GDBPrior
#' @rdname GDBPrior
#' @usage \S4method{GDBPrior}{character}(variable)
methods::setMethod(
  "GDBPrior",
  methods::signature(variable = "character"),
  function(variable, hyper = 'increasing', ... ) {
    assertthat::assert_that(!missing(variable),!missing(hyper),
                            msg = 'Variable or constrain unset.')
    assertthat::assert_that(
      is.character(variable),
      is.character(hyper)
    )
    assertthat::assert_that(length(variable)==1,msg = 'More than one prior variable supplied. Use GDBPriors')
    # Match supplied constrain in case someone has been lazy
    hyper <- match.arg(hyper, c('increasing', 'decreasing','convex', 'concave','positive', 'negative', 'none'), several.ok = FALSE)

    # Create new prior object
    bdproto(
      'GDBPrior',
      Prior,
      id = new_id(),
      variable = variable,
      value = hyper
    )
  }
)

#' Helper function when multiple variables are supplied
#'
#' @name GDBPriors
#' @rdname GDBPriors
#' @exportMethod GDBPriors
#' @inheritParams GDBPrior
#' @keywords priors
#' @family prior
#' @export
methods::setGeneric(
  "GDBPriors",
  signature = methods::signature("variable"),
  function(variable, hyper = 'increasing', ...) standardGeneric("GDBPriors"))

#' @name GDBPriors
#' @rdname GDBPriors
#' @usage \S4method{GDBPriors}{character}(variable)
methods::setMethod(
  "GDBPriors",
  methods::signature(variable = "character"),
  function(variable, hyper = 'increasing', ... ) {
    assertthat::assert_that(!missing(variable),!missing(hyper),
                            msg = 'Variable or constrain unset.')
    assertthat::assert_that(
      is.character(variable),
      is.character(hyper), length(hyper)==1
    )
    assertthat::assert_that(length(variable)>1, msg = 'Only one prior variable supplied. Use GDBPrior')
    # Match supplied constrain in case someone has been lazy
    hyper <- match.arg(hyper, c('increasing', 'decreasing','convex', 'concave','positive', 'negative', 'none'),several.ok = FALSE)

    multiple_priors <- list()
    for(k in variable){
      np <- GDBPrior(variable = k,hyper = hyper)
      multiple_priors[[as.character(np$id)]] <- np
    }
    return(multiple_priors)
  }
)
