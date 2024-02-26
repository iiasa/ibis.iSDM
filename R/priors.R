#' @include class-priorlist.R class-prior.R
NULL

#' Creates a new PriorList object
#'
#' @description A [`PriorList`] object is essentially a list that contains
#' individual [`Prior-class`] objects. In order to use priors for any of the
#' engines, the respective [`Prior-class`] has to be identified (e.g.
#' [`INLAPrior`]) and embedded in a [`PriorList`] object. Afterwards these
#' objects can then be added to a [distribution] object with the [add_priors]
#' function.
#'
#' @param x A [`Prior-class`] object added to the list.
#' @param ... One or multiple additional [`Prior-class`] object added to the list.
#'
#' @returns A [`PriorList`] object.
#'
#' @seealso [`Prior-class`], [`PriorList-class`]
#' @family prior
#'
#' @examples
#' p1 <- GDBPrior(variable = "Forest", hyper = "positive")
#' p2 <- GDBPrior(variable = "Urban", hyper = "decreasing")
#' priors(p1, p2)
#'
#' @name priors
NULL

#' @rdname priors
#' @export
methods::setGeneric(
  "priors",
  signature = methods::signature('x'),
  function(x,...) standardGeneric("priors"))

#' @rdname priors
methods::setMethod(
  "priors",
  methods::signature(x = "ANY"),
  function(x,...) {
    # Get supplied matched arguments
    assertthat::assert_that(
      inherits(x,'Prior') || is.list(x),
      msg = 'Supply a prior object.'
    )
    # Capture dot objects and check that they are all of class Prior
    dots <- list(...)
    # if(is.list(dots[[1]])) dots <- unlist(dots)

    assertthat::assert_that(
      all( sapply(dots, function(x) inherits(x, 'Prior')) ),
      msg = 'Only prior objects are supported.'
    )

    # Single entry
    if(is.list(x)){
      ll <- x
      # Set names of list entries to their id
      names(ll) <- sapply(ll, function(z) z$id )
    } else {
      # Prior object list
      ll <- list()
      ll[[as.character(x$id)]] <- x
      for(var in dots) ll[[as.character(var$id)]] <- var
    }
    # Ensure that names of the list are identical to the prior ids
    assertthat::assert_that(
      all( sapply(ll, function(z) z$id ) %in% names(ll) )
    )
    assertthat::assert_that(all( sapply(ll, function(z) nchar(z$variable)) > 0 ),
                            msg = 'Empty variable name. Possible error?')

    # Assess duplicates of priors in order they appear in the list
    vars <- vapply(ll, function(z) z$variable, character(1))
    types <- sapply(ll, function(z) z$type)
    # Replace any Gaussian with normal
    if(any(types == "gaussian")) types[which(types=="gaussian")] <- "normal"

    # Check that there are no NA values among the variables
    assertthat::assert_that(!anyNA(vars))

    # Remove duplicated SPDE priors for INLA
    if(any(vars == "spde")){
      if(any(duplicated(cbind(vars,types)))){
        ll <- ll[-which(duplicated(cbind(vars,types),fromLast = TRUE))]
        warning(paste0('Ignoring duplicated prior(s) for: ', paste0(vars[which(duplicated(vars))],collapse = ', ')))
      }
    } else {
      # Remove duplicated variables, taking only the one added
      if(anyDuplicated(vars)>0){
        if(getOption('ibis.setupmessages',default = TRUE)) myLog('[Setup]','yellow','Found duplicated prior variables. Taking last one added.')
        ll <- ll[-which(duplicated(vars,fromLast = TRUE))]
      }
    }

    # Create new priorlist object and supply all objects here
    pl <- PriorList$new(priors = ll)
    return(pl)
  }
)

#' Creates a new PriorList object
#'
#' @description A [`PriorList`] object is essentially a list that contains
#' individual [`Prior-class`] objects. In order to use priors for any of the
#' engines, the respective [`Prior-class`] has to be identified (e.g.
#' [`INLAPrior`]) and embedded in a [`PriorList`] object. Afterwards these
#' objects can then be added to a [distribution] object with the [add_priors]
#' function.
#'
#' @param ... One or multiple additional [`Prior-class`] object added to the list.
#'
#' @returns A [`PriorList`] object.
#'
#' @seealso [`Prior-class`], [`PriorList-class`]
#' @family prior
#'
#' @examples
#' \dontrun{
#' p1 <- INLAPrior(variable = "Forest",type = "normal", hyper = c(1,1e4))
#' p2 <- INLAPrior(variable = "Urban",type = "normal", hyper = c(0,1e-2))
#' priors(p1, p2)
#' }
#'
#' @name priors
NULL

#' @rdname priors
#' @export
methods::setGeneric(
  "priors",
  signature = methods::signature('x'),
  function(x,...) standardGeneric("priors"))

#' @rdname priors
methods::setMethod(
  "priors",
  methods::signature(x = "ANY"),
  function(x,...) {
    # Get supplied matched arguments
    assertthat::assert_that(
      inherits(x,'Prior') || is.list(x),
      msg = 'Supply a prior object.'
    )
    # Capture dot objects and check that they are all of class Prior
    dots <- list(...)
    if(length(dots)>0) dots <- unlist(dots)

    assertthat::assert_that(
      all( sapply(dots, function(x) inherits(x, 'Prior')) ),
      msg = 'Only prior objects are supported.'
    )

    if(is.list(x)){
      ll <- x
      # Set names of list entries to their id
      names(ll) <- sapply(ll, function(z) z$id )
      # Combine all priors objects together if multiple are supplied
      if(length(dots)>0) for(var in dots) ll[[as.character(var$id)]] <- var
    } else {
      # Prior object list
      ll <- list()
      ll[[as.character(x$id)]] <- x
      for(var in dots) ll[[as.character(var$id)]] <- var
    }
    # Ensure that names of the list are identical to the prior ids
    assertthat::assert_that(
      all( sapply(ll, function(z) z$id ) %in% names(ll) )
    )
    assertthat::assert_that(all( sapply(ll, function(z) nchar(z$variable)) > 0 ),
                            msg = 'Empty variable name. Possible error?')

    # Assess duplicates of priors in order they appear in the list
    vars <- vapply(ll, function(z) z$variable, character(1))
    types <- sapply(ll, function(z) z$type)
    # Replace any Gaussian with normal
    if(any(types == "gaussian")) types[which(types=="gaussian")] <- "normal"

    # Check that there are no NA values among the variables
    assertthat::assert_that(!anyNA(vars))

    # Remove duplicated SPDE priors for INLA
    if(any(vars == "spde")){
      if(any(duplicated(cbind(vars,types)))){
        ll <- ll[-which(duplicated(cbind(vars,types),fromLast = TRUE))]
        warning(paste0('Ignoring duplicated prior(s) for: ', paste0(vars[which(duplicated(vars))],collapse = ', ')))
      }
    } else {
      # Remove duplicated variables, taking only the one added
      if(anyDuplicated(vars)>0){
        if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Setup]','yellow','Found duplicated prior variables. Taking last one added.')
        ll <- ll[-which(duplicated(vars,fromLast = TRUE))]
      }
    }

    # Create new priorlist object and supply all objects here
    pl <- PriorList$new(priors = ll)
    return(pl)
  }
)
