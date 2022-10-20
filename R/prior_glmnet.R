#' @include utils.R bdproto.R bdproto-prior.R
NULL

#' Placeholder for GLMNET
#'
#' @description
#' The [`engine_glmnet`] engine does not support priors in a typical sense, however it is
#' possible to specify so called penalty factors on all variables in the model.
#'
#' The default penalty is \code{1} for each coefficient, i.e. coefficients are penalized equally.
#' In contrast a variable with penalty.factor equal to \code{0} is not penalized at all.
#'
#' For a regularized regression that supports a few more options on the priors, check out the
#' Bayesian [`engine_breg`].
#'
#' @param variable A [`character`] variable passed on to the prior object.
#' @param hyper A [`numeric`] value between \code{0} and \code{1} that state the penalization factor.
#' By default this is set to \code{0}, implying the \code{"variable"} provided is not regularized at all.
#' @param ... Variables passed on to prior object.
#' @seealso [`Prior-class`]
#' @examples
#' \dontrun{
#' # Retain variable
#' p1 <- GLMNETPrior(variable = "forest", hyper = 0)
#' p1
#' # Smaller chance to be regularized
#' p2 <- GLMNETPrior(variable = "forest", hyper = 0.2)
#' p2
#' }
#' @keywords priors
#' @family prior
#' @aliases GLMNETPrior
#' @name GLMNETPrior
NULL

#' @name GLMNETPrior
#' @rdname GLMNETPrior
#' @exportMethod GLMNETPrior
#' @export
methods::setGeneric(
  "GLMNETPrior",
  signature = methods::signature("variable"),
  function(variable, hyper = 0, ...) standardGeneric("GLMNETPrior"))

#' @name GLMNETPrior
#' @rdname GLMNETPrior
#' @usage \S4method{GLMNETPrior}{character, numeric}(variable, hyper)
methods::setMethod(
  "GLMNETPrior",
  methods::signature(variable = "character"),
  function(variable, hyper = 0, ... ) {
    assertthat::assert_that(!missing(variable),
                            is.numeric(hyper),
                            hyper >=0, hyper <=1,
                            msg = 'Variable or hyper unset unset.')

    assertthat::assert_that(length(variable)==1,
                            msg = 'More than one prior variable supplied. Use BREGPriors')

    # Create new prior object
    bdproto(
      'GLMNETPrior',
      Prior,
      id = new_id(),
      variable = variable,
      value = hyper
    )
  }
)

#' Helper function when multiple variables are supplied
#' @name GLMNETPriors
#' @description
#' This is a helper function to specify several [GLMNETPrior] with the same
#' hyper-parameters, but different variables.
#' @rdname GLMNETPriors
#' @exportMethod GLMNETPriors
#' @inheritParams GLMNETPrior
#' @family prior
#' @keywords priors
#' @export
methods::setGeneric(
  "GLMNETPriors",
  signature = methods::signature("variable"),
  function(variable, hyper = 0) standardGeneric("GLMNETPriors"))

#' @name GLMNETPriors
#' @rdname GLMNETPriors
#' @usage \S4method{GLMNETPriors}{character, numeric}(variable, hyper)
methods::setMethod(
  "GLMNETPriors",
  methods::signature(variable = "character"),
  function(variable, hyper = 0) {
    assertthat::assert_that(!missing(variable),
                            msg = 'Variable not set.')
    assertthat::assert_that(
      is.character(variable),
      is.null(hyper) || is.numeric(hyper)
    )
    assertthat::assert_that(length(variable)>1, msg = 'Only one prior variable supplied. Use GLMNETPrior')

    # Check that hyper is between 0 and 1
    assertthat::assert_that(hyper >= 0, hyper <= 1)

    multiple_priors <- list()
    for(k in variable){
      np <- GLMNETPrior(variable = k, hyper = hyper)
      multiple_priors[[as.character(np$id)]] <- np
    }
    return(multiple_priors)
  }
)
