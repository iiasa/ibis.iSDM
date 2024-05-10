#' @include class-biodiversitydistribution.R class-predictors.R
NULL

#' Add predictions from a fitted model to a Biodiversity distribution object
#'
#' @description This function is a convenience wrapper to add the output from a
#' previous fitted [`DistributionModel`] to another [`BiodiversityDistribution-class`]
#' object. Obviously only works if a prediction was fitted in the model. Options
#' to instead add thresholds, or to transform / derivate the model outputs are
#' also supported.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param model A [`DistributionModel`] object.
#' @param transform A [`vector`] stating whether predictors should be preprocessed
#' in any way (Options: \code{'none'},\code{'pca'}, \code{'scale'}, \code{'norm'})
#' @param derivates A Boolean check whether derivate features should be considered
#' (Options: \code{'none'}, \code{'thresh'}, \code{'hinge'}, \code{'quad'}) )
#' @param threshold_only A [`logical`] flag indicating whether to add thresholded
#' layers from the fitted model (if existing) instead (Default: \code{FALSE}).
#' @param priors A [`PriorList-class`] object. Default is set to \code{NULL} which
#' uses default prior assumptions.
#' @param ... Other parameters passed down
#'
#' @details
#' A transformation takes the provided rasters and for instance rescales them or
#' transforms them through a principal component analysis ([prcomp]). In
#' contrast, derivates leave the original provided predictors alone, but instead
#' create new ones, for instance by transforming their values through a
#' quadratic or hinge transformation. Note that this effectively increases the
#' number of predictors in the object, generally requiring stronger
#' regularization by the used [`Engine`]. Both transformations and derivates can
#' also be combined. Available options for transformation are:
#' * \code{'none'} - Leaves the provided predictors in the original scale.
#' * \code{'pca'} - Converts the predictors to principal components. Note that this
#' results in a renaming of the variables to principal component axes!
#' * \code{'scale'} - Transforms all predictors by applying [scale] on them.
#' * \code{'norm'} - Normalizes all predictors by transforming them to a scale from 0 to 1.
#' * \code{'windsor'} - Applies a windsorization to the target predictors. By default
#' this effectively cuts the predictors to the 0.05 and 0.95, thus helping to
#' remove extreme outliers.
#'
#' Available options for creating derivates are:
#' * \code{'none'} - No additional predictor derivates are created.
#' * \code{'quad'} - Adds quadratic transformed predictors.
#' * \code{'interaction'} - Add interacting predictors. Interactions need to be specified (\code{"int_variables"})!
#' * \code{'thresh'} - Add threshold transformed predictors.
#' * \code{'hinge'} - Add hinge transformed predictors.
#' * \code{'bin'} - Add predictors binned by their percentiles.
#'
#' @examples
#' \dontrun{
#'  # Fit first model
#'  fit <- distribution(background) |>
#'         add_predictors(covariates) |>
#'         add_biodiversity_poipa(species) |>
#'         engine_glmnet() |>
#'         train()
#'
#'  # New model object
#'  obj <- distribution(background) |>
#'         add_predictors_model(fit)
#'  obj
#' }
#'
#' @name add_predictors_model
NULL

#' @rdname add_predictors_model
#' @export
methods::setGeneric(
  "add_predictors_model",
  signature = methods::signature("x", "model"),
  function(x, model, transform = 'scale', derivates = 'none',
           threshold_only = FALSE, priors = NULL, ...) standardGeneric("add_predictors_model"))

#' @rdname add_predictors_model
methods::setMethod(
  "add_predictors_model",
  methods::signature(x = "BiodiversityDistribution", model = "ANY"),
  function(x, model, transform = 'scale', derivates = 'none',
           threshold_only = FALSE, priors = NULL, ...) {
    # Try and match transform and derivatives arguments
    transform <- match.arg(transform, c('none','pca', 'scale', 'norm', 'windsor') , several.ok = TRUE)
    derivates <- match.arg(derivates, c('none','thresh', 'hinge', 'quadratic', 'bin', 'interaction') , several.ok = TRUE)

    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            inherits(model, "DistributionModel"),
                            is.logical(threshold_only),
                            all(transform == 'none') || all( transform %in% c('pca', 'scale', 'norm', 'windsor') ),
                            all(derivates == 'none') || all( derivates %in% c('thresh', 'hinge', 'quadratic', 'bin', 'interaction') ),
                            is.null(priors) || inherits(priors,'PriorList')
    )
    # Messenger
    if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Setup]','green','Adding predictors from fitted model...')

    # Make a clone copy of the object
    y <- x$clone(deep = TRUE)

    # If priors have been set, save them in the distribution object
    if(!is.null(priors)) {
      y <- y$set_priors(priors)
    }

    # Get prediction from model object
    assertthat::assert_that(
      "prediction" %in% model$show_rasters(),
      msg = "No prediction found in model object!"
    )
    if(threshold_only){
      tr <- grep('threshold', model$show_rasters(),ignore.case = TRUE,value = TRUE)
      if(length(tr)==1){
        prediction <- model$get_data(tr)
      } else {
        if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Setup]','yellow','No threshold found in fitted model. Using prediction...')
        prediction <- model$get_data()
      }
    } else {
      prediction <- model$get_data()
    }
    assertthat::assert_that(is.Raster(prediction))
    # Set names
    names(prediction) <- paste0(make.names(model$model$runname),"__",names(prediction))

    # Standardization and scaling
    if('none' %notin% transform){
      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Setup]','green','Transforming prediction...')
      for(tt in transform) prediction <- predictor_transform(prediction, option = tt)
    }
    assertthat::assert_that(is.Raster(prediction))

    # Calculate derivates if set
    if('none' %notin% derivates){
      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Setup]','green','Creating prediction derivates...')
      # Specific condition for interaction
      new_prediction <- terra::rast()
      for(dd in derivates){
        suppressWarnings(
          new_prediction <- c(new_prediction, predictor_derivate(prediction,
                                                                 option = dd) )
        )
      }
      # Add
      prediction <- c(prediction, new_prediction)
      rm(new_prediction)
    }

    # Assign an attribute to this object to keep track of it
    attr(prediction,'transform') <- transform

    # Sanitize names if specified
    if(getOption('ibis.cleannames', default = TRUE)) names(prediction) <- sanitize_names(names(prediction))

    # Get existing predictors
    if(!is.Waiver(x$predictors)){
      env <- y$predictors$get_data()
      env <- c(env, prediction)
    }

    # Finally set the data to the BiodiversityDistribution object
    pd <- PredictorDataset$new(id = new_id(),
                               data = env,
                               ...)
    y$set_predictors(pd)
  }
)
