#' @include utils.R bdproto.R bdproto-prior.R bdproto-priorlist.R
NULL

#' Set priors to an existing distribution object
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param prior A [`Prior-class`] object containing a single prior
#' @param priors A [`PriorList-class`] object containing multiple priors
#' @param ... Other parameters passed down

#' @details TBD
#' @section Notes:
#' @aliases add_predictors
#' @references
#'
#' @examples
#' \dontrun{
#'  TBD
#' }
#' @name set_priors
NULL

#' @name set_priors
#' @rdname set_priors
#' @exportMethod set_priors
#' @export
methods::setGeneric(
  "set_priors",
  signature = methods::signature("x"),
  function(x, prior = NULL, priors = NULL, ...) standardGeneric("set_priors"))

#' @name set_priors
#' @rdname set_priors
#' @usage \S4method{set_priors}{BiodiversityDistribution}(x)
methods::setMethod(
  "set_priors",
  methods::signature(x = "BiodiversityDistribution"),
  function(x, prior = NULL, priors = NULL, ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.null(prior) || inherits(prior, c('INLAPrior', 'GDBPrior')),
                            is.null(priors) || inherits(priors, "PriorList")
    )

    #FIXME: Ideally check whether provided prior is already in object. Either per id or per variable match

    if(!is.null(prior)){
      x <- x$set_priors(
        priors(prior)
      )
    }
    if(!is.null(priors)){
      x <- x$set_priors( priors )
    }
    # Return x with newly added priors
    x
  }
)
