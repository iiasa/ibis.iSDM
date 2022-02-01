#' @include utils.R bdproto.R bdproto-prior.R
NULL

#' Create a new spike and slab prior for Bayesian generalized linear models
#'
#' @description
#' Function to include prior information via Zellner-style spike and slab prior for
#' generalized linear models used in [engine_breg]. These priors are similar to
#' the horseshoe priors used in regularized [engine_stan] models.
#'
#' @details
#' The Zellner-style spike and slab prior for generalized linear models are specified as described
#' in the \pkg{Boom} R-package. Currently supported are two options which work for models with
#' Poisson and binomial (bernoulli) distributed errors.
#' Two types of priors can be provided on a variable:
#' * \code{"coefficient"} Allows to specify Gaussian priors on the coefficient mean and
#' precision (!) for the model. Either the mean or the mean and precision has to be supplied.
#' * \code{"inclusion.probability"} A [`vector`] giving the prior probability of inclusion for the
#' specified variable. This can be useful when prior information on preference is known but not the
#' strength of it.
#'
#' @param variable A [`character`] matched against existing predictors or latent effects.
#' @param type A [`character`] specifying the type of prior to be set. Can be set either to
#' \code{"coefficient"} or \code{"inclusion.probability"}. See details.
#' @param hyper A [`numeric`] estimate of the mean regression coefficients
#' @param ... Variables passed on to prior object.
#' @references
#' * Hugh Chipman, Edward I. George, Robert E. McCulloch, M. Clyde, Dean P. Foster, Robert A. Stine (2001), "The Practical Implementation of Bayesian Model Selection" Lecture Notes-Monograph Series, Vol. 38, pp. 65-134. Institute of Mathematical Statistics.
#' @seealso [`Prior-class`]
#' @family prior
#' @keywords priors
#' @aliases BREGPrior
#' @name BREGPrior
NULL

#' @name BREGPrior
#' @rdname BREGPrior
#' @exportMethod BREGPrior
#' @export
methods::setGeneric(
  "BREGPrior",
  signature = methods::signature("variable", "type", "hyper"),
  function(variable, type, hyper) standardGeneric("BREGPrior"))

#' @name BREGPrior
#' @rdname BREGPrior
#' @usage \S4method{BREGPrior}{character}(variable)
methods::setMethod(
  "BREGPrior",
  methods::signature(variable = "character", type = "character", hyper = "numeric"),
  function(variable, type, hyper) {
    assertthat::assert_that(!missing(variable),!missing(hyper),
                            msg = 'Variable or constrain unset.')
    assertthat::assert_that(
      is.character(variable),
      is.character(type),
      is.numeric(hyper)
    )
    assertthat::assert_that(length(variable)==1,msg = 'More than one prior variable supplied. Use BREGPriors')
    # Match supplied type in case someone has been lazy
    type <- match.arg(type,  c("coefficient", "inclusion.probability"), several.ok = FALSE)

    # Create new prior object
    bdproto(
      'BREGPrior',
      Prior,
      id = new_id(),
      variable = variable,
      type = type,
      value = hyper
    )
  }
)

#' Helper function when multiple variables are supplied
#' @name BREGPriors
#' @description
#' This is a helper function to specify several [BREGPrior] with the same
#' hyper-parameters, but different variables.
#' @rdname BREGPriors
#' @exportMethod BREGPriors
#' @inheritParams BREGPrior
#' @family prior
#' @keywords priors
#' @export
methods::setGeneric(
  "BREGPriors",
  signature = methods::signature("variable", "type", "hyper"),
  function(variable, type, hyper) standardGeneric("BREGPriors"))

#' @name BREGPriors
#' @rdname BREGPriors
#' @usage \S4method{BREGPriors}{character}(variable, type, hyper)
methods::setMethod(
  "BREGPriors",
  methods::signature(variable = "character", type = "character", hyper = "numeric"),
  function(variable, type, hyper ) {
    assertthat::assert_that(!missing(variable),!missing(hyper),
                            msg = 'Variable or constrain unset.')
    assertthat::assert_that(
      is.character(variable),
      is.character(type),
      is.numeric(hyper)
    )
    assertthat::assert_that(length(variable)>1, msg = 'Only one prior variable supplied. Use BREGPrior')
    # Match supplied type in case someone has been lazy
    type <- match.arg(type,  c("coefficient", "inclusion.probability"), several.ok = FALSE)

    multiple_priors <- list()
    for(k in variable){
      np <- BREGPrior(variable = k,type = type, hyper = hyper)
      multiple_priors[[as.character(np$id)]] <- np
    }
    return(multiple_priors)
  }
)
