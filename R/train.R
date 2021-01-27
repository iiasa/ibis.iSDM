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

    # --- #
    # Set model object for fitting
    model <- list()

    # Get biodiversity data
    model[['data']] <- list()
    types <- names( x$biodiversity$get_types() )
    for(ty in types) model[['data']][[ty]] <- x$biodiversity$get_data(ty)

    # Get covariates
    if(is.Waiver(x$predictor_names())) {
      # Dummy covariate of background raster
      dummy <- x$background; names(dummy) <- 'dummy'
      model[['predictors']] <- dummy
      } else { model[['predictors']] <- x$predictors$get_data() }

    # assign default priors
    if(is.Waiver( x$priors )){
      # TODO: Define prior objects
      message('TODO: Define prior objects')
      model['priors'] <- NULL
      #x$set_prior( add_default_priors() )
    }
    # Assign default formula
    if(is.Waiver(x$equation)){
      message('TODO: Support custom formula')
      model$
      x$get_biodiversity_equations()
    }

    # Generate an id
    out$id <- new_id()

    # Train the object
    out <- x$engine$train()

    # return output object
    return(out)
  }
)
