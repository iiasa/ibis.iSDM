#' @include class-prior.R
NULL

#' Monotonic constrained priors for boosted regressions
#'
#' @description Monotonic constrains for gradient descent boosting models do not
#' work in the same way as other priors where a specific coefficient or
#' magnitude of importance is specified. Rather monotonic constraints
#' **enforce** a specific directionality of regression coefficients so that for
#' instance a coefficient has to be positive or negative.
#'
#' __Important:__ Specifying a monotonic constrain for the [engine_gdb] does not
#' guarantee that the variable is retained in the model as it can still be
#' regularized out.
#'
#' @param variable A [`character`] matched against existing predictors variables.
#' @param hyper A [`character`] object describing the type of constrain. Available
#' options are \code{'increasing'}, \code{'decreasing'}, \code{'convex'}, \code{'concave'},
#' \code{'positive'}, \code{'negative'} or \code{'none'}.
#' @param ... Variables passed on to prior object.
#'
#' @note Similar priors can also be defined for the [`engine_xgboost`] via
#' [`XGBPrior()`].
#'
#' @references
#' * Hofner, B., Müller, J., & Hothorn, T. (2011). Monotonicity‐constrained species
#'  distribution models. Ecology, 92(10), 1895-1901.
#'
#' @seealso [`Prior-class`], [`XGBPrior`]
#' @keywords priors
#' @family prior
#'
#' @name GDBPrior
NULL

#' @rdname GDBPrior
#' @export
methods::setGeneric(
  "GDBPrior",
  signature = methods::signature("variable"),
  function(variable, hyper = 'increasing', ...) standardGeneric("GDBPrior"))

#' @rdname GDBPrior
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

    # Sanitize names if specified
    if(getOption('ibis.cleannames', default = TRUE)) variable <- sanitize_names(variable)

    # Create new prior object
    pp <- Prior$new(
      name = 'GDBPrior',
      id = new_id(),
      variable = variable,
      value = hyper
    )
    return(pp)
  }
)

#' Helper function when multiple variables are supplied for a GDB prior
#'
#' @description This is a helper function to specify several [GLMNETPrior] with
#' the same hyper-parameters, but different variables.
#'
#' @inheritParams GDBPrior
#'
#' @keywords priors
#' @family prior
#'
#' @name GDBPriors
NULL

#' @rdname GDBPriors
#' @export
methods::setGeneric(
  "GDBPriors",
  signature = methods::signature("variable"),
  function(variable, hyper = 'increasing', ...) standardGeneric("GDBPriors"))

#' @rdname GDBPriors
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
    assertthat::assert_that(length(variable)>1,
                            msg = 'Only one prior variable supplied. Use GDBPrior')
    # Match supplied constrain in case someone has been lazy
    hyper <- match.arg(hyper,
                       c('increasing', 'decreasing','convex',
                         'concave','positive', 'negative', 'none'),
                       several.ok = FALSE)

    multiple_priors <- list()
    for(k in variable){
      np <- GDBPrior(variable = k,hyper = hyper)
      multiple_priors[[as.character(np$id)]] <- np
    }
    return(multiple_priors)
  }
)
