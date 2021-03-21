#' @include utils.R bdproto.R bdproto-priorlist.R bdproto-prior.R
NULL

#' Creates a new PriorList object that contains Priors
#'
#' @param x A [`Prior-class`] object added to the list
#' @param ... One or multiple additional [`Prior-class`] object added to the list
#' @details TBD
#' @section Notes:
#' @references
#'
#' @seealso [`Prior-class`], [`PriorList-class`]
#'
#' @aliases priorlist, priors
#' @name priorlist
NULL

#' @name priorlist
#' @rdname priorlist
#' @exportMethod priorlist
#' @export
methods::setGeneric(
  "priorlist",
  signature = methods::signature('x'),
  function(x,...) standardGeneric("priorlist"))

#' @name priorlist
#' @rdname priorlist
#' @usage \S4method{priorlist}{Prior}(x)
methods::setMethod(
  "priorlist",
  methods::signature(x = "ANY"),
  function(x,...) {
    # Get supplied matched arguments
    assertthat::assert_that(
      inherits(x,'Prior'),
      msg = 'Supply a prior object.'
    )
    # Capture dot objects and check that they are all of class Prior
    dots <- list(...)
    assertthat::assert_that(
      all( sapply(dots, function(x) inherits(x, 'Prior')) ),
      msg = 'Only prior objects are supported.'
    )

    # Prior object list
    ll <- list()
    ll[[as.character(x$id)]] <- x
    for(var in dots) ll[[as.character(var$id)]] <- var

    # Create new priorlist object and supply all objects here
    bdproto(
      NULL,
      PriorList,
      priors = ll
      )

  }
)

#' Creates a new PriorList object that contains Priors
#'
#' @aliases priorlist, priors
#' @name priors
#' @rdname priors
#' @exportMethod priors
#' @export
methods::setGeneric(
  "priors",
  signature = methods::signature('x'),
  function(x,...) standardGeneric("priors")
)
# FIXME: Ensure that everything is properly inherited from priorList
#' @name priors
#' @rdname priors
#' @usage \S4method{priors}{Prior}(x)
methods::setMethod(
  "priors",
  methods::signature(x = "ANY"),
  function(x,...) priorlist(x, ...)
)
