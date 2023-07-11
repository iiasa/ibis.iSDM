#' @include utils.R bdproto.R bdproto-prior.R
NULL

#' Regression penalty priors for GLMNET
#'
#' @description
#' The [`engine_glmnet`] engine does not support priors in a typical sense, however it is
#' possible to specify so called penalty factors as well as lower and upper limits
#' on all variables in the model.
#'
#' The default penalty multiplier is \code{1} for each coefficient X covariate,
#' i.e. coefficients are penalized equally and then informed by an intersection
#' of any absence information with the covariates.
#' In contrast a variable with penalty.factor equal to \code{0} is not penalized at all.
#'
#' In addition, it is possible to specifiy a lower and upper limit for specific coefficients,
#' which constrain them to a certain range.
#' By default those ranges are set to \code{-Inf} and \code{Inf} respectively, but can be reset to a
#' specific value range by altering \code{"lims"} (see examples).
#'
#' For a regularized regression that supports a few more options on the priors, check out the
#' Bayesian [`engine_breg`].
#'
#' @param variable A [`character`] variable passed on to the prior object.
#' @param hyper A [`numeric`] value between \code{0} and \code{1} that state the penalization factor.
#' By default this is set to \code{0}, implying the \code{"variable"} provided is not regularized at all.
#' @param lims A [`numeric`] [`vector`] of the lower and upper limits for each coefficient (Default: \code{c(-Inf, Inf)}).
#' @param ... Variables passed on to prior object.
#' @seealso [`Prior-class`]
#' @examples
#' \dontrun{
#' # Retain variable
#' p1 <- GLMNETPrior(variable = "forest", hyper = 0)
#' p1
#' # Smaller chance to be regularized
#' p2 <- GLMNETPrior(variable = "forest", hyper = 0.2, lims = c(0, Inf))
#' p2
#' }
#' @keywords priors
#' @family prior
#' @aliases GLMNETPrior, Prior
#' @name GLMNETPrior
NULL

#' @name GLMNETPrior
#' @rdname GLMNETPrior
#' @exportMethod GLMNETPrior
#' @export
methods::setGeneric(
  "GLMNETPrior",
  signature = methods::signature("variable"),
  function(variable, hyper = 0, lims = c(-Inf, Inf), ...) standardGeneric("GLMNETPrior"))

#' @name GLMNETPrior
#' @rdname GLMNETPrior
#' @usage \S4method{GLMNETPrior}{character,numeric,ANY}(variable,hyper,lims,...)
methods::setMethod(
  "GLMNETPrior",
  methods::signature(variable = "character"),
  function(variable, hyper = 0, lims = c(-Inf, Inf), ... ) {
    assertthat::assert_that(!missing(variable),
                            is.numeric(hyper),
                            hyper >=0, hyper <=1,
                            is.numeric(lims),
                            length(lims)==2,
                            msg = 'Variable or hyper unset unset.')

    assertthat::assert_that(length(variable)==1,
                            msg = 'More than one prior variable supplied. Use BREGPriors')

    assertthat::assert_that(
      # Check that lower limit is negative or zero
      sign(lims[1])==0 || sign(lims[1])==-1,
      # Check that upper limit is positive or zero
      sign(lims[2])==0 || sign(lims[2])==1,
      msg = "Lower or upper limits specified incorrectly! They have to be negative/positive or zero each."
                            )

    # Sanitize names if specified
    if(getOption('ibis.cleannames')) variable <- sanitize_names(variable)

    # Create new prior object
    bdproto(
      'GLMNETPrior',
      Prior,
      id = new_id(),
      variable = variable,
      value = hyper,
      lims = lims
    )
  }
)

#' Helper function when multiple variables are supplied for a GLMNET prior
#' @name GLMNETPriors
#' @description
#' This is a helper function to specify several [GLMNETPrior] with the same
#' hyper-parameters, but different variables.
#' @rdname GLMNETPriors
#' @exportMethod GLMNETPriors
#' @inheritParams GLMNETPrior
#' @aliases GLMNETPriors
#' @family prior
#' @keywords priors
#' @export
methods::setGeneric(
  "GLMNETPriors",
  signature = methods::signature("variable"),
  function(variable, hyper = 0, lims = c(-Inf, Inf)) standardGeneric("GLMNETPriors"))

#' @name GLMNETPriors
#' @rdname GLMNETPriors
#' @usage \S4method{GLMNETPriors}{character,numeric,ANY}(variable,hyper,lims,...)
methods::setMethod(
  "GLMNETPriors",
  methods::signature(variable = "character"),
  function(variable, hyper = 0, lims = c(-Inf, Inf)) {
    assertthat::assert_that(!missing(variable),
                            msg = 'Variable not set.')
    assertthat::assert_that(
      is.character(variable),
      is.null(hyper) || is.numeric(hyper),
      is.numeric(lims)
    )
    assertthat::assert_that(length(variable)>1, msg = 'Only one prior variable supplied. Use GLMNETPrior')

    # Check that hyper is between 0 and 1
    assertthat::assert_that(hyper >= 0, hyper <= 1)

    multiple_priors <- list()
    for(k in variable){
      np <- GLMNETPrior(variable = k, hyper = hyper, lims = lims)
      multiple_priors[[as.character(np$id)]] <- np
    }
    return(multiple_priors)
  }
)
