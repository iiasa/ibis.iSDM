#' @include utils.R bdproto-biodiversitydistribution.R utils-spatial.R
NULL

#' Train the model from a given engine
#'
#' Train a [distribution()] model with the specified engine.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object).
#'
#' @param ... further arguments passed on.
#'
#' @details
#'
#' @return A distribution prediction object
#'
#' @name train
#'
#' @exportMethod train
#'
#' @aliases train, train-method
#'
#' @export
NULL

#' @name train
#' @rdname train
#' @exportMethod train
#' @export
methods::setGeneric(
  "train",
  signature = methods::signature("x"),
  function(x, runname = NULL,...) standardGeneric("train"))

#' @name train
#' @rdname train
#' @usage \S4method{train}{BiodiversityDistribution}(x)
methods::setMethod(
  "train",
  methods::signature(x = "BiodiversityDistribution"),
  function(x, runname = NULL, ...) {
    # Make load checks
    assertthat::assert_that(
      inherits(x, "BiodiversityDistribution"),
      is.null(runname) || is.character(runname)
    )
    # Now make checks on completeness of the object
    assertthat::assert_that(!is.Waiver(x$engine),
                            msg = 'No engine set for training the distribution model.')
    assertthat::assert_that( x$show_biodiversity_length() > 0,
                             msg = 'No biodiversity data specified.')

    # assign default priors
    if(is.Waiver( x$priors )){
      # TODO: Define prior objects
      #x$set_prior( add_default_priors() )
    }
    # Assign default formula
    if(is.Waiver(x$equation)){
      # FIXME: Needs tuning so that different formulas are supported
      x$set_equation('observation ~ .')
    }

    # Generate an id
    out$id <- new_id()

    # Train the object
    out <- x$engine$train()

    # return output object
    return(out)
  }
)
