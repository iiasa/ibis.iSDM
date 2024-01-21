#' @include class-prior.R
NULL

#' Create a new STAN prior
#'
#' @description Function to create a new prior for [engine_stan] models. Priors
#' currently can be set on specific environmental predictors.
#'
#' @param variable A [`character`] matched against existing predictors or latent
#' effects.
#' @param type A [`character`] specifying the type of prior to be set.
#' @param hyper A [`vector`] with [`numeric`] values to be used as hyper parameters.
#' First entry is treated as mean (Default: \code{0}), the second as the standard
#' variation (Default: \code{2}) of a Gaussian distribution on the respective coefficient.
#' @param ... Variables passed on to prior object.
#'
#' @references
#' * Lemoine, N. P. (2019). Moving beyond noninformative priors: why and how to
#' choose weakly informative priors in Bayesian analyses. Oikos, 128(7), 912-928.
#' * Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B., Betancourt, M.,
#'  ... & Riddell, A. (2017). Stan: A probabilistic programming language. Journal
#'  of statistical software, 76(1), 1-32.
#'
#' @seealso [`Prior-class`].
#' @family prior
#' @keywords priors
#'
#' @examples
#' \dontrun{
#'  pp <- STANPrior("forest", "normal", c(0,1))
#' }
#'
#' @name STANPrior
NULL

#' @rdname STANPrior
#' @export
methods::setGeneric(
  "STANPrior",
  signature = methods::signature("variable", "type"),
  function(variable, type, hyper = c(0, 2), ...) standardGeneric("STANPrior"))

#' @rdname STANPrior
methods::setMethod(
  "STANPrior",
  methods::signature(variable = "character", type = "character"),
  function(variable, type, hyper = c(0, 2), ... ) {
    assertthat::assert_that(!missing(variable), !missing(type),
                            msg = 'Variable or type unset.')
    assertthat::assert_that(
      is.character(variable),
      is.character(type),
      is.vector(hyper), length(hyper) == 2,
      all(is.numeric(hyper))
    )
    # Match supplied type in case someone has been lazy
    type <- match.arg(type, c('gaussian', 'normal'), several.ok = FALSE)

    # Sanitize names if specified
    if(getOption('ibis.cleannames')) variable <- sanitize_names(variable)

    # Create new prior object
    bdproto(
      'STANPrior',
      Prior,
      id = new_id(),
      type = type,
      variable = variable,
      value = hyper,
      # plot
      plot = function(self){
        return(NULL) # TODO: Could implement some QoL plotting here
      }
    )
  }
)

#' Helper function when multiple variables and types are supplied for STAN
#'

#' @description This is a helper function to specify several [STANPrior] with
#' the same hyper-parameters, but different variables.
#'
#' @param variables A [`vector`] of [`character`] matched against existing
#' predictors or latent effects.
#' @param type A [`character`] specifying the type of prior to be set.
#' @param hyper A [`vector`] with [`numeric`] values to be used as hyper-parameters.
#' @param ... Variables passed on to prior object
#'
#' @family prior
#' @keywords priors
#'
#' @name STANPriors
NULL

#' @rdname STANPriors
#' @export
methods::setGeneric(
  "STANPriors",
  signature = methods::signature("variables", "type"),
  function(variables, type, hyper = c(0, 2), ...) standardGeneric("STANPriors"))

#' @rdname STANPriors
methods::setMethod(
  "STANPriors",
  methods::signature(variables = "vector", type = "character"),
  function(variables, type, hyper = c(0, 2), ... ) {
    assertthat::assert_that(!missing(variables),!missing(type),
                            msg = 'Variable or Type unset.')
    assertthat::assert_that(
      is.vector(variables), is.character(type),
      is.vector(hyper), length(hyper) == 2,
      all(is.numeric(hyper))
    )

    assertthat::assert_that(length(variables)>1, msg = 'Only one prior variable supplied. Use STANPrior')
    # Match supplied type in case someone has been lazy
    type <- match.arg(type, c("gaussian", "normal"), several.ok = FALSE)

    multiple_priors <- list()
    for(k in variables){
      np <- STANPrior(variable = k, type = type, hyper = hyper)
      multiple_priors[[as.character(np$id)]] <- np
    }
    return(multiple_priors)
  }
)
