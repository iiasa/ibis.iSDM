#' @include class-prior.R
NULL

#' Create a new spike and slab prior for Bayesian generalized linear models
#'
#' @description Function to include prior information via Zellner-style spike
#' and slab prior for generalized linear models used in [engine_breg]. These
#' priors are similar to the horseshoe priors used in regularized [engine_stan]
#' models and penalize regressions by assuming most predictors having an effect
#' of \code{0}.
#'
#' @param variable A [`character`] matched against existing predictors.
#' @param hyper A [`numeric`] estimate of the mean regression coefficients.
#' @param ip A [`numeric`] estimate between 0 and 1 of the inclusion probability
#' of the target variable (Default: \code{NULL}).
#'
#' @details The Zellner-style spike and slab prior for generalized linear models
#' are specified as described in the \pkg{Boom} R-package. Currently supported
#' are two options which work for models with \code{Poisson} and \code{binomial}
#' (\code{Bernoulli}) distributed errors. Two types of priors can be provided on
#' a variable:
#' * \code{"coefficient"} Allows to specify Gaussian priors on the mean coefficients of the model.
#' Priors on the coefficients can be provided via the \code{"hyper"} parameter.
#' Note that variables with such a prior can still be regularized out from the
#' model.
#' * \code{"inclusion.probability"} A [`vector`] giving the prior probability of inclusion for the
#' specified variable. This can be useful when prior information on preference
#' is known but not the strength of it.
#'
#' If coefficients are set, then the inclusion probability is also modified by
#' default. However even when not knowing a particular estimate of a beta
#' coefficients and their direction, one can still provide an estimate of the
#' inclusion probability. In other words:
#' **The hyperparameters 'hyper' and 'ip' can't be both \code{NULL}.**
#'
#' @references
#' * Hugh Chipman, Edward I. George, Robert E. McCulloch, M. Clyde, Dean P. Foster,
#' Robert A. Stine (2001), "The Practical Implementation of Bayesian Model Selection"
#' Lecture Notes-Monograph Series, Vol. 38, pp. 65-134. Institute of Mathematical Statistics.
#'
#' @seealso [`Prior-class`]
#' @family prior
#' @keywords priors
#'
#' @examples
#' \dontrun{
#' # Positive coefficient
#' p1 <- BREGPrior(variable = "forest", hyper = 2, ip = NULL)
#' p1
#' # Coefficient and direction unknown but variable def. important
#' p2 <- BREGPrior(variable = "forest", hyper = NULL, ip = 1)
#' p2
#' }
#'
#' @name BREGPrior
NULL

#' @rdname BREGPrior
#' @export
methods::setGeneric(
  "BREGPrior",
  signature = methods::signature("variable"),
  function(variable, hyper = NULL, ip = NULL) standardGeneric("BREGPrior"))

#' @rdname BREGPrior
methods::setMethod(
  "BREGPrior",
  methods::signature(variable = "character"),
  function(variable, hyper = NULL, ip = NULL) {
    assertthat::assert_that(!missing(variable),
                            msg = 'Variable or constrain unset.')
    assertthat::assert_that(
      is.character(variable),
      is.null(hyper) || is.numeric(hyper),
      is.null(ip) || is.numeric(ip)
    )
    assertthat::assert_that(length(variable)==1,msg = 'More than one prior variable supplied. Use BREGPriors')
    # Both hyper and ip can not both be NULL
    if(is.null(hyper) && is.null(ip)) stop("Set either hyper and/or ip to a non-null value!")

    # Sanitize names if specified
    if(getOption('ibis.cleannames')) variable <- sanitize_names(variable)

    # Check that ip is between 0 and 1
    if(!is.null(ip)) assertthat::assert_that(ip >= 0, ip <= 1)
    # Create new prior object
    bdproto(
      'BREGPrior',
      Prior,
      id = new_id(),
      variable = variable,
      value = hyper,
      prob = ip
    )
  }
)

#' Helper function when multiple variables are supplied for a BREG prior

#' @description This is a helper function to specify several [BREGPrior] with
#' the same hyper-parameters, but different variables.
#'
#' @inheritParams BREGPrior
#'
#' @family prior
#' @keywords priors
#'
#' @name BREGPriors
NULL

#' @rdname BREGPriors
#' @export
methods::setGeneric(
  "BREGPriors",
  signature = methods::signature("variable"),
  function(variable, hyper = NULL, ip = NULL) standardGeneric("BREGPriors"))

#' @rdname BREGPriors
methods::setMethod(
  "BREGPriors",
  methods::signature(variable = "character"),
  function(variable, hyper = NULL, ip = NULL) {
    assertthat::assert_that(!missing(variable),
                            msg = 'Variable not set.')
    assertthat::assert_that(
      is.character(variable),
      is.null(hyper) || is.numeric(hyper),
      is.null(ip) || is.numeric(ip)
    )
    assertthat::assert_that(length(variable)>1, msg = 'Only one prior variable supplied. Use BREGPrior')
    # Both hyper and ip can not both be NULL
    if(is.null(hyper) && is.null(ip)) stop("Set either hyper and/or ip to a non-null value!")

    # Check that ip is between 0 and 1
    if(!is.null(ip)) assertthat::assert_that(ip >= 0, ip <= 1)

    multiple_priors <- list()
    for(k in variable){
      np <- BREGPrior(variable = k, hyper = hyper, ip = ip)
      multiple_priors[[as.character(np$id)]] <- np
    }
    return(multiple_priors)
  }
)
