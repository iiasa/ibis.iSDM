#' @include utils.R bdproto.R bdproto-prior.R
NULL

#' Create a tree-based split probability prior for BART
#'
#' @description
#' Function to include prior information as split probability for the
#' Bayesian additive regression tree model [engine_bart].
#'
#' Priors for bart will be specified as transition probabilities of variables used to generate splits
#' in the regression tree. These can be numeric and between \code{0} and \code{1}.
#'
#' @param variable A [`character`] matched against existing predictors or latent effects.
#' @param hyper A [`numeric`] object with a number being \code{>0} and equal to \code{1}. Defaults to \code{0.75}.
#' @param ... Variables passed on to prior object.
#' @references
#' * Chipman, H., George, E., and McCulloch, R. (2009) BART: Bayesian Additive Regression Trees.
#' * Chipman, H., George, E., and McCulloch R. (2006) Bayesian Ensemble Learning. Advances in Neural Information Processing Systems 19, Scholkopf, Platt and Hoffman, Eds., MIT Press, Cambridge, MA, 265-272.
#' @seealso [`Prior-class`].
#' @family prior
#' @keywords priors
#' @aliases BARTPrior
#' @name BARTPrior
NULL

#' @name BARTPrior
#' @rdname BARTPrior
#' @exportMethod BARTPrior
#' @export
methods::setGeneric(
  "BARTPrior",
  signature = methods::signature("variable"),
  function(variable, hyper = 0.75, ...) standardGeneric("BARTPrior"))

#' @name BARTPrior
#' @rdname BARTPrior
#' @usage \S4method{BARTPrior}{character}(variable)
methods::setMethod(
  "BARTPrior",
  methods::signature(variable = "character"),
  function(variable, hyper = 0.75, ... ) {
    assertthat::assert_that(!missing(variable),!missing(hyper),
                            msg = 'Variable or prior probability unset.')
    assertthat::assert_that(
      is.character(variable),
      is.numeric(hyper),
      (hyper>=0 && hyper <= 1)
    )
    assertthat::assert_that(length(variable)==1,msg = 'More than one prior variable supplied. Use XGBPriors')

    # Create new prior object
    bdproto(
      'BARTPrior',
      Prior,
      id = new_id(),
      variable = variable,
      value = hyper
    )
  }
)

#' Helper function when multiple variables are supplied
#' @name BARTPriors
#' @description
#' This is a helper function to specify several [BARTPrior] objects with the same
#' hyper-parameters, but different variables.
#' @rdname BARTPriors
#' @exportMethod BARTPriors
#' @inheritParams BARTPrior
#' @family prior
#' @keywords priors
#' @export
methods::setGeneric(
  "BARTPriors",
  signature = methods::signature("variable"),
  function(variable, hyper = 0.75, ...) standardGeneric("BARTPriors"))

#' @name BARTPriors
#' @rdname BARTPriors
#' @usage \S4method{BARTPriors}{character}(variable)
methods::setMethod(
  "BARTPriors",
  methods::signature(variable = "character"),
  function(variable, hyper = 0.75, ... ) {
    assertthat::assert_that(!missing(variable),!missing(hyper),
                            msg = 'Variable or prior probability unset.')
    assertthat::assert_that(
      is.character(variable),
      is.numeric(hyper), length(hyper)==1,
      (hyper>=0 && hyper <= 1)
    )
    assertthat::assert_that(length(variable)>1, msg = 'Only one prior variable supplied. Use BARTPrior')

    multiple_priors <- list()
    for(k in variable){
      np <- BARTPrior(variable = k,hyper = hyper)
      multiple_priors[[as.character(np$id)]] <- np
    }
    return(multiple_priors)
  }
)