#' @include class-prior.R
NULL

#' Create a new monotonic prior for boosted regressions
#'
#' @description Function to include prior information as monotonic constrain to
#' an extreme gradient descent boosting model [`engine_xgboost`]. Monotonic
#' priors enforce directionality in direction of certain variables, however
#' specifying a monotonic constrain does not guarantee that the variable is not
#' regularized out during model fitting.
#'
#' @param variable A [`character`] matched against existing predictors or latent
#' effects.
#' @param hyper A [`character`] object describing the type of constrain. Available
#' options are \code{'increasing'}, \code{'decreasing'}, \code{'positive'},
#' \code{'negative'}, \code{'none'}.
#' @param ... Variables passed on to prior object.
#'
#' @references
#' * Chen, T., He, T., Benesty, M., Khotilovich, V., Tang, Y., & Cho, H. (2015).
#' Xgboost: extreme gradient boosting. R package version 0.4-2, 1(4), 1-4.
#'
#' @seealso [`Prior-class`] and [`GDBPrior`].
#' @family prior
#' @keywords priors
#'
#' @examples
#' \dontrun{
#'  pp <- XGBPrior("forest", "increasing")
#' }
#'
#' @name XGBPrior
NULL

#' @rdname XGBPrior
#' @export
methods::setGeneric(
  "XGBPrior",
  signature = methods::signature("variable", "hyper"),
  function(variable, hyper = 'increasing', ...) standardGeneric("XGBPrior"))

#' @rdname XGBPrior
methods::setMethod(
  "XGBPrior",
  methods::signature(variable = "character", hyper = "character"),
  function(variable, hyper = 'increasing', ... ) {
    assertthat::assert_that(!missing(variable),!missing(hyper),
                            msg = 'Variable or constrain unset.')
    assertthat::assert_that(
      is.character(variable),
      is.character(hyper)
    )
    assertthat::assert_that(length(variable)==1,msg = 'More than one prior variable supplied. Use XGBPriors')
    # Match supplied constrain in case someone has been lazy
    hyper <- match.arg(hyper, c('increasing', 'decreasing','positive', 'negative', 'none'), several.ok = FALSE)

    # Sanitize names if specified
    if(getOption('ibis.cleannames', default = TRUE)) variable <- sanitize_names(variable)

    # Create new prior object
    pp <- Prior$new(
      name = 'XGBPrior',
      id = new_id(),
      variable = variable,
      value = hyper
    )
    return(pp)
  }
)

#' Helper function when multiple variables are supplied for XGBoost priors

#' @description This is a helper function to specify several [XGBPrior] with the
#' same hyper-parameters, but different variables.
#'
#' @inheritParams XGBPrior
#'
#' @family prior
#' @keywords priors
#'
#' @name XGBPriors
NULL

#' @rdname XGBPriors
#' @export
methods::setGeneric(
  "XGBPriors",
  signature = methods::signature("variable"),
  function(variable, hyper = 'increasing', ...) standardGeneric("XGBPriors"))

#' @rdname XGBPriors
methods::setMethod(
  "XGBPriors",
  methods::signature(variable = "character"),
  function(variable, hyper = 'increasing', ... ) {
    assertthat::assert_that(!missing(variable),!missing(hyper),
                            msg = 'Variable or constrain unset.')
    assertthat::assert_that(
      is.character(variable),
      is.character(hyper), length(hyper)==1
    )
    assertthat::assert_that(length(variable)>1, msg = 'Only one prior variable supplied. Use XGBPrior')
    # Match supplied constrain in case someone has been lazy
    hyper <- match.arg(hyper, c('increasing', 'decreasing','positive', 'negative', 'none'), several.ok = FALSE)

    multiple_priors <- list()
    for(k in variable){
      np <- XGBPrior(variable = k,hyper = hyper)
      multiple_priors[[as.character(np$id)]] <- np
    }
    return(multiple_priors)
  }
)

#' Create a new interaction prior for XGBoost
#'
#' @description Function to include prior information as interaction constraints
#' in an extreme gradient descent boosting model [`engine_xgboost`]. Interaction
#' priors define groups of variables that are allowed to interact in the same
#' tree path. Variables outside the same group are not allowed to interact.
#'
#' @param variables A [`character`] vector matched against existing predictors or
#' latent effects after XGBoost preprocessing.
#' @param ... Variables passed on to prior object.
#'
#' @details XGBoost interaction constraints are only supported by tree boosters.
#' They can be combined with monotonic constraints supplied through [`XGBPrior`].
#'
#' @seealso [`Prior-class`], [`XGBPrior`] and [`engine_xgboost`].
#' @family prior
#' @keywords priors
#'
#' @examples
#' \dontrun{
#'  pp <- XGBInteractionPrior(c("forest", "temperature"))
#' }
#'
#' @name XGBInteractionPrior
NULL

#' @rdname XGBInteractionPrior
#' @export
methods::setGeneric(
  "XGBInteractionPrior",
  signature = methods::signature("variables"),
  function(variables, ...) standardGeneric("XGBInteractionPrior"))

#' @rdname XGBInteractionPrior
methods::setMethod(
  "XGBInteractionPrior",
  methods::signature(variables = "character"),
  function(variables, ... ) {
    assertthat::assert_that(!missing(variables),
                            msg = 'Interaction variables unset.')
    assertthat::assert_that(
      is.character(variables),
      length(variables) > 0,
      all(nchar(variables) > 0),
      msg = 'Supply at least one non-empty interaction variable.'
    )
    # Sanitize names if specified
    if(getOption('ibis.cleannames', default = TRUE)) variables <- sanitize_names(variables)
    assertthat::assert_that(!anyDuplicated(variables),
                            msg = 'Duplicated variables in interaction prior.')

    # Create a unique synthetic variable name so PriorList deduplication does
    # not replace monotone priors for the same real predictor.
    id <- new_id()
    variable <- paste0("xgb_interaction_", as.character(id))

    # Create new prior object
    pp <- Prior$new(
      name = 'XGBInteractionPrior',
      id = id,
      variable = variable,
      value = variables
    )
    return(pp)
  }
)

#' Helper function when multiple interaction groups are supplied for XGBoost
#'
#' @description This is a helper function to specify several
#' [XGBInteractionPrior] objects from a list of variable groups.
#'
#' @param groups A [`list`] of [`character`] vectors. Each vector is one allowed
#' interaction group.
#' @param ... Variables passed on to prior object.
#'
#' @family prior
#' @keywords priors
#'
#' @name XGBInteractionPriors
NULL

#' @rdname XGBInteractionPriors
#' @export
methods::setGeneric(
  "XGBInteractionPriors",
  signature = methods::signature("groups"),
  function(groups, ...) standardGeneric("XGBInteractionPriors"))

#' @rdname XGBInteractionPriors
methods::setMethod(
  "XGBInteractionPriors",
  methods::signature(groups = "list"),
  function(groups, ... ) {
    assertthat::assert_that(!missing(groups),
                            msg = 'Interaction groups unset.')
    assertthat::assert_that(
      length(groups) > 0,
      all(vapply(groups, is.character, logical(1))),
      msg = 'Supply a non-empty list of character vectors.'
    )

    multiple_priors <- list()
    for(group in groups){
      np <- XGBInteractionPrior(variables = group)
      multiple_priors[[as.character(np$id)]] <- np
    }
    return(multiple_priors)
  }
)
