#' @include utils.R bdproto.R bdproto-prior.R bdproto-priorlist.R
NULL

#' Add priors to an existing distribution object
#'
#' @description
#' This function simply allows to add priors to an existing [distribution] object.
#' The supplied priors must be a [`PriorList-class`] object created through
#' calling [priors].
#' @note
#' Alternatively priors to environmental predictors can also directly added as parameter
#' via [add_predictors]
#' @param x [distribution] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param priors A [`PriorList-class`] object containing multiple priors.
#' @param ... Other parameters passed down.
#' @family prior
#' @aliases add_priors
#' @examples
#' \dontrun{
#'  x <- distribution(background)
#' }
#' @name add_priors
NULL

#' @name add_priors
#' @rdname add_priors
#' @exportMethod add_priors
#' @export
methods::setGeneric(
  "add_priors",
  signature = methods::signature("x"),
  function(x, priors = NULL, ...) standardGeneric("add_priors"))

#' @name add_priors
#' @rdname add_priors
#' @usage \S4method{add_priors}{BiodiversityDistribution}(x)
methods::setMethod(
  "add_priors",
  methods::signature(x = "BiodiversityDistribution"),
  function(x, priors = NULL ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.null(priors) || inherits(priors, "PriorList") || inherits(priors, 'INLAPrior') || inherits(priors, 'GDBPrior')
    )

    # Convert to prior list object
    if(inherits(priors, 'INLAPrior') || inherits(priors, 'GDBPrior')){
      priors <- priors(priors)
    }
    if(!is.null(priors)){
      x <- x$set_priors( priors )
    }
    # Return x with newly added priors
    x
  }
)

#' @name set_priors
#' @inherit add_priors
#' @inheritParams add_priors
#' @keywords deprecated
methods::setGeneric(
  "set_priors",
  signature = methods::signature("x"),
  function(x, priors = NULL, ...) standardGeneric("set_priors"))

#' @name set_priors
#' @inherit add_priors
#' @inheritParams add_priors
#' @keywords deprecated
#' @usage \S4method{set_priors}{BiodiversityDistribution}(x)
methods::setMethod(
  "set_priors",
  methods::signature(x = "BiodiversityDistribution"),
  function(x, priors = NULL ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.null(priors) || inherits(priors, "PriorList") || inherits(priors, 'INLAPrior') || inherits(priors, 'GDBPrior')
    )
    message('Deprecated. Use add_priors ')
    add_priors(x, priors)
  }
)

#' Remove existing priors from an existing distribution object
#'
#' @description
#' This function allows to remove priors from an existing [distribution] object.
#' In order to remove a set prior, the name of the prior has to be specified.
#' @param x [distribution] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param names A [`vector`] or [`character`] object for priors to be removed.
#' @param ... Other parameters passed down
#' @aliases rm_priors
#' @family prior
#' @examples
#' \dontrun{
#'  TBD
#' }
#' @name rm_priors
NULL

#' @name rm_priors
#' @rdname rm_priors
#' @exportMethod rm_priors
#' @export
methods::setGeneric(
  "rm_priors",
  signature = methods::signature("x"),
  function(x, names = NULL, ...) standardGeneric("rm_priors"))

#' @name rm_priors
#' @rdname rm_priors
#' @usage \S4method{rm_priors}{BiodiversityDistribution}(x)
methods::setMethod(
  "rm_priors",
  methods::signature(x = "BiodiversityDistribution"),
  function(x, names = NULL ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.null(names) || is.vector(names) || is.character(names)
    )
    x <- x$rm_priors(names)
    # Return x with newly added priors
    x
  }
)
