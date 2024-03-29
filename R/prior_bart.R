#' @include class-prior.R
NULL

#' Create a tree-based split probability prior for BART
#'
#' @description Function to include prior information as split probability for
#' the Bayesian additive regression tree model added via [engine_bart].
#'
#' Priors for [engine_bart] have to be specified as transition probabilities of
#' variables which are internally used to generate splits in the regression
#' tree. Specifying a prior can thus help to 'enforce' a split with a given
#' variable. These can be numeric and coded as values between \code{0} and
#' \code{1}.
#'
#' @param variable A [`character`] matched against existing predictors or latent
#' effects.
#' @param hyper A [`numeric`] object with a number being \code{>0} and equal to
#' \code{1}. Defaults to \code{0.75}.
#' @param ... Variables passed on to prior object.
#'
#' @note Even if a given variable is included as split in the regression or
#' classification tree, this does not necessarily mean that the prediction
#' changes if the value is non-informative (as the split can occur early on). It
#' does however affect any variable importance estimates calculated from the
#' model.
#'
#' @seealso [`Prior-class`].
#' @family prior
#' @keywords priors
#'
#' @references
#' * Chipman, H., George, E., and McCulloch, R. (2009) BART: Bayesian Additive
#' Regression Trees.
#' * Chipman, H., George, E., and McCulloch R. (2006) Bayesian Ensemble Learning.
#' Advances in Neural Information Processing Systems 19, Scholkopf, Platt and Hoffman,
#' Eds., MIT Press, Cambridge, MA, 265-272.
#'
#' @name BARTPrior
NULL

#' @rdname BARTPrior
#' @export
methods::setGeneric(
  "BARTPrior",
  signature = methods::signature("variable"),
  function(variable, hyper = 0.75, ...) standardGeneric("BARTPrior"))

#' @rdname BARTPrior
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

    # Sanitize names if specified
    if(getOption('ibis.cleannames', default = TRUE)) variable <- sanitize_names(variable)

    # Create new prior object
    pp <- Prior$new(
      name = 'BARTPrior',
      id = new_id(),
      variable = variable,
      value = hyper
    )
    return(pp)
  }
)

#' Helper function when multiple variables are supplied for a BART prior
#'
#' @description This is a helper function to specify several [BARTPrior] objects
#' with the same hyper-parameters, but different variables.
#'
#' @inheritParams BARTPrior
#'
#' @family prior
#' @keywords priors
#'
#' @name BARTPriors
NULL

#' @rdname BARTPriors
#' @export
methods::setGeneric(
  "BARTPriors",
  signature = methods::signature("variable"),
  function(variable, hyper = 0.75, ...) standardGeneric("BARTPriors"))

#' @rdname BARTPriors
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
