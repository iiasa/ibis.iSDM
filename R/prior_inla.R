#' @include utils.R bdproto.R bdproto-prior.R
NULL

#' Create a new INLA prior
#'
#' @description
#' For any fixed and random effect Inla supports a range of different priors of exponential distributions.
#' Unfortunately it does not support popular priors related to parameter regularization such as
#' Laplace or Horseshoe priors, which limits the capability of [`engine_inla`] for regularization.
#'
#' Currently this package only implements priors related to \code{"linear"} coefficients
#' that can be added to environmental predictors.
#'
#' In addition it is also supports to specification of penalized complexity priors which can be added
#' to a SPDE spatial random effect added via [`add_latent_spatial()`].
#' Here the range of the penalized complexity prior can be specified through \code{'prior.range'} and
#' the uncertainty via \code{'prior.sigma'} both supplied to the options 'type' and 'hyper'.
#' Other priors available in INLA \code{ names(INLA::inla.models()$prior) ) } might also work, but have not
#' been tested!
#'
#' @param variable A [`character`] matched against existing predictors or latent effects.
#' @param type A [`character`] specifying the type of prior to be set.
#' @param hyper A [`vector`] with [`numeric`] values to be used as hyper-parameters. First entry is treated as mean, the second as precision.
#' Note that precision is not equivalent (rather the inverse) to typical standard deviation specified in Gaussian priors.
#' @param ... Variables passed on to prior object.
#' @note
#' Currently only linear (\code{'linear'}) and constrained \code{'clinear'} linear priors are supported.
#' @references
#' * Rue, H., Riebler, A., Sørbye, S. H., Illian, J. B., Simpson, D. P., & Lindgren, F. K. (2017). Bayesian computing with INLA: a review. Annual Review of Statistics and Its Application, 4, 395-421.
#' * Simpson, D., Rue, H., Riebler, A., Martins, T. G., & Sørbye, S. H. (2017). Penalising model component complexity: A principled, practical approach to constructing priors. Statistical science, 32(1), 1-28.
#' @seealso [`Prior-class`].
#' @keywords priors
#' @family prior
#' @aliases INLAPrior
#' @name INLAPrior
NULL

#' @name INLAPrior
#' @rdname INLAPrior
#' @exportMethod INLAPrior
#' @export
methods::setGeneric(
  "INLAPrior",
  signature = methods::signature("variable", "type"),
  function(variable, type = "linear", hyper = c(0, 0.001), ...) standardGeneric("INLAPrior"))

# TODO: Consider other INLA priors as well as they become relevant for the modelling, `names(INLA::inla.models()$prior)`
#' @name INLAPrior
#' @rdname INLAPrior
#' @usage \S4method{INLAPrior}{character, character}(variable, type)
methods::setMethod(
  "INLAPrior",
  methods::signature(variable = "character", type = "character"),
  function(variable, type = "linear", hyper = c(0, 0.001), ... ) {
    assertthat::assert_that(!missing(variable), !missing(type),
                            msg = 'Variable or type unset.')
    assertthat::assert_that(
      is.character(variable),
      is.character(type),
      is.vector(hyper), length(hyper) == 2,
      all(is.numeric(hyper))
    )
    # FIXME: Further checks for mean and precision and whether they are in realistic bounds
    # assertthat::assert_that(
    #   all(hyper > -1e2) && all(hyper < 1e2), # Check bounds
    #   msg = 'Potentially unrealistic priors specified.'
    # )
    # Match supplied type in case someone has been lazy
    type <- match.arg(type, c('clinear','prior.range','prior.sigma',names(INLA::inla.models()$prior) ), several.ok = FALSE)

    # Check bounds for clinear being correct
    if(type=='clinear') assertthat::assert_that( hyper[1] < hyper[2])

    # Rename prior?
    # FIXME: This assumes that values have been supplied in this order. Implement checks

    # Other supplied arguments
    #args <- as.list(match.call())

    # Create new prior object
    bdproto(
      'INLAPrior',
      Prior,
      id = new_id(),
      type = type,
      variable = variable,
      value = hyper,
      # FIXME:
      # https://stats.stackexchange.com/questions/350235/how-to-convert-estimated-precision-to-variance
      # distribution = function(mean, prec) stats::rnorm(n = 1000, mean = mean, sd = prec),

      # Custom function to format INLA priors to be used in hyper
      format = function(self,type){
        assertthat::assert_that(type %in% c('mean','prec','theta'))
        ol <- list()
        ol[[type]] <-  list(prior = self$type, param = self$value)
        return(ol)
      }
    )
  }
)

#' Helper function when multiple variables and types are supplied
#' @description
#' This is a helper function to specify several [INLAPrior] objects with the same
#' hyper-parameters, but different variables.
#' @name INLAPriors
#' @param variables A [`vector`] of [`characters`] matched against existing predictors or latent effects.
#' @param type A [`character`] specifying the type of prior to be set.
#' @param hyper A [`vector`] with [`numeric`] values to be used as hyper-parameters.
#' @param ... Variables passed on to prior object.
#' @rdname INLAPriors
#' @exportMethod INLAPriors
#' @keywords priors
#' @family prior
#' @export
methods::setGeneric(
  "INLAPriors",
  signature = methods::signature("variables", "type"),
  function(variables, type, hyper = c(0, 0.001), ...) standardGeneric("INLAPriors"))

#' @name INLAPriors
#' @rdname INLAPriors
#' @usage \S4method{INLAPriors}{vector, character}(variables, type)
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
    type <- match.arg(type, names(INLA::inla.models()$prior), several.ok = FALSE)

    multiple_priors <- list()
    for(k in variables){
      np <- INLAPrior(variable = k,type = type, hyper = hyper)
      multiple_priors[[as.character(np$id)]] <- np
    }
    return(multiple_priors)
  }
)


