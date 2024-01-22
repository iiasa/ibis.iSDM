#' @include class-prior.R
NULL

#' Create a new INLA prior
#'
#' @description For any fixed and random effect INLA supports a range of
#' different priors of exponential distributions.
#'
#' Currently supported for INLA in ibis.iSDM are the following priors that can
#' be specified via \code{"type"}:
#'
#' * \code{"normal"} or \code{"gaussian"}: Priors on normal distributed and set to specified variable. Required parameters
#' are a mean and a precision estimate provided to \code{"hyper"}. Note that
#' precision is not equivalent (rather the inverse) to typical standard
#' deviation specified in Gaussian priors. Defaults are set to a mean of
#' \code{0} and a precision of \code{0.001}.
#'
#' * \code{"clinear"}: Prior that places a constraint on the linear coefficients of a model
#' so as that the coefficient is in a specified interval
#' \code{"c(lower,upper)"}. Specified through hyper these values can be
#' negative, positive or infinite.
#'
#' * \code{"spde"}, specifically \code{'prior.range'} and \code{'prior.sigma'}: Specification of
#' penalized complexity priors which can be added to a SPDE spatial random
#' effect added via [`add_latent_spatial()`]. Here the range of the penalized
#' complexity prior can be specified through \code{'prior.range'} and the
#' uncertainty via \code{'prior.sigma'} both supplied to the options 'type' and
#' 'hyper'.
#'
#' Other priors available in INLA \code{ names(INLA::inla.models()$prior) ) }
#' might also work, but have not been tested!
#'
#' @param variable A [`character`] matched against existing predictors or latent
#' effects.
#' @param type A [`character`] specifying the type of prior to be set.
#' @param hyper A [`vector`] with [`numeric`] values to be used as hyper-parameters.
#' See description. The default values are set to a mean of \code{0} and a precision of \code{0.001}.
#' @param ... Variables passed on to prior object.
#'
#' @note Compared to other engines, INLA does unfortunately does not support
#' priors related to more stringent parameter regularization such as Laplace or
#' Horseshoe priors, which limits the capability of [`engine_inla`] for
#' regularization. That being said many of the default uninformative priors act
#' already regularize the coefficients to some degree.
#'
#' @references
#' * Rue, H., Riebler, A., Sørbye, S. H., Illian, J. B., Simpson, D. P., & Lindgren, F. K. (2017).
#' Bayesian computing with INLA: a review. Annual Review of Statistics and Its Application, 4, 395-421.
#' * Simpson, D., Rue, H., Riebler, A., Martins, T. G., & Sørbye, S. H. (2017).
#' Penalising model component complexity: A principled, practical approach to constructing
#' priors. Statistical science, 32(1), 1-28.
#'
#' @seealso [`Prior-class`].
#' @keywords priors
#' @family prior
#'
#' @name INLAPrior
NULL

#' @rdname INLAPrior
#' @export
methods::setGeneric(
  "INLAPrior",
  signature = methods::signature("variable", "type"),
  function(variable, type = "normal", hyper = c(0, 0.001), ...) standardGeneric("INLAPrior"))

#' @rdname INLAPrior
methods::setMethod(
  "INLAPrior",
  methods::signature(variable = "character", type = "character"),
  function(variable, type = "normal", hyper = c(0, 0.001), ... ) {
    assertthat::assert_that(!missing(variable), !missing(type),
                            msg = 'Variable or type unset.')
    assertthat::assert_that(
      is.character(variable),
      is.character(type),
      is.vector(hyper), length(hyper) == 2,
      all(is.numeric(hyper))
    )

    # Match supplied type in case someone has been lazy
    type <- match.arg(type, c('clinear','prior.range','prior.sigma', names(INLA::inla.models()$prior) ), several.ok = FALSE)

    # Sanitize names if specified
    if(getOption('ibis.cleannames')) variable <- sanitize_names(variable)

    # Create new prior object
    pp <- Prior$new(
      name = 'INLAPrior',
      id = new_id(),
      type = type,
      variable = variable,
      value = hyper,
    )
    # FIXME:
    # https://stats.stackexchange.com/questions/350235/how-to-convert-estimated-precision-to-variance
    # distribution = function(mean, prec) stats::rnorm(n = 1000, mean = mean,
    # sd = prec),
      # Custom function to format INLA priors to be used in hyper
    pp$format = function(self,type){
        assertthat::assert_that(type %in% c('mean','prec','theta'))
        ol <- list()
        ol[[type]] <-  list(prior = self$type, param = self$value)
        return(ol)
    }
    return(pp)
  }
)

#' Helper function when multiple variables and types are supplied for INLA
#'
#' @description This is a helper function to specify several [INLAPrior] objects
#' with the same hyper-parameters, but different variables.
#'
#' @param variables A [`vector`] of [`character`] matched against existing
#' predictors or latent effects.
#' @param type A [`character`] specifying the type of prior to be set.
#' @param hyper A [`vector`] with [`numeric`] values to be used as hyper-parameters.
#' @param ... Variables passed on to prior object.
#'
#' @keywords priors
#' @family prior
#'
#' @name INLAPriors
NULL

#' @rdname INLAPriors
#' @export
methods::setGeneric(
  "INLAPriors",
  signature = methods::signature("variables", "type"),
  function(variables, type, hyper = c(0, 0.001), ...) standardGeneric("INLAPriors"))

#' @rdname INLAPriors
methods::setMethod(
  "INLAPriors",
  methods::signature(variables = "vector", type = "character"),
  function(variables, type, hyper = c(0, 0.001), ... ) {
    assertthat::assert_that(!missing(variables),!missing(type),
                            msg = 'Variable or Type unset.')
    assertthat::assert_that(
      is.vector(variables), is.character(type),
      is.vector(hyper), length(hyper) == 2,
      all(is.numeric(hyper))
    )

    assertthat::assert_that(length(variables)>1, msg = 'Only one prior variable supplied. Use INLAPrior')
    # Match supplied type in case someone has been lazy
    gg <- try({find.package("INLA")}, silent = TRUE)
    if(!inherits(gg, "try-error")){
      pn <- names(INLA::inla.models()$prior)
    } else {
      pn <- c("normal", "gaussian","gamma","flat","pc", "pc.range","pc.prec")
    }
    type <- match.arg(type, pn, several.ok = FALSE)

    multiple_priors <- list()
    for(k in variables){
      np <- INLAPrior(variable = k,type = type, hyper = hyper)
      multiple_priors[[as.character(np$id)]] <- np
    }
    return(multiple_priors)
  }
)


