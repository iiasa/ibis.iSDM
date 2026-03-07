#' Add a control to a BiodiversityModel object to train ensembles of small models
#'
#' @description
#' Estimating species distributions with few occurrences can be challenging. This is
#' however often the case for many rare species, where seldom more than 5-10 occurrences
#' are observed, or species that are estimated over a relatively coarse grid (e.g.,
#' aggregating observations over a coarse grid). In these cases, it is often that
#' more complex algorithms tend to overpredict/overfit any given distributions (Breiner et al. 2015,
#' 2018).
#' A solution to this is to use ensembles of small models (ESM) that are essentially
#' multiple separate models fitted on a subset of the covariates supplied. After inference,
#' each model is then combined via [`ensemble()`] to produce a final distribution and coefficients.
#'
#' In a range of studies it has been shown that this is an effective way to model
#' species with few occurrences, and that it can outperform more complex models (Breiner
#' et al. 2018, Erickson and Smith 2023).
#'
#' See Details for the implementation in this package.
#'
#' @details
#' To make ensembles of small models (ESM) work, the [`distribution`] object gets
#' assigned a flag that indicates that the model should be trained as an ensemble of small models
#'  within [`train()`]. The number of variables in each small model can be controlled
#' via the `n_covs` parameter.
#'
#' After inference and prediction, the coefficients from all models are extracted
#' and also combined in an ensemble. A dummy [`engine_glm()`] model is then trained
#' with the exact supplied coefficients per covariate. Note that this necessarily
#' implies that the coefficients are not directly comparable to those of a single model.
#' Further only linear modelling is supported to (re)capture the coefficients.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param n_covs A [`numeric`] on the number of covariates in each small model (Default: \code{2}).
#'
#' @note
#' Apart from ensembles of small models, consider also filtering the predictors
#' prior to supplying them a [`distribution`] object. This can be done using via
#' the function [`predictor_filter()`].
#'
#' @returns Adds a control option to a [`distribution`] object that introduces
#' the modelling of ensembles of small models in [`train()`].
#'
#' @references
#' * Breiner, F. T., Guisan, A., Bergamini, A., & Nobis, M. P. (2015).
#' Overcoming limitations of modelling rare species by using ensembles of small models.
#' Methods in Ecology and Evolution, 6(10), 1210-1218.
#' * Breiner, F. T., Nobis, M. P., Bergamini, A., & Guisan, A. (2018).
#' Optimizing ensembles of small models for predicting the distribution of species with few occurrences.
#' Methods in Ecology and Evolution, 9(4), 802-808.
#' * Erickson, K. D., & Smith, A. B. (2023).
#' Modeling the rarest of the rare: a comparison between multi‐species distribution models,
#' ensembles of small models, and single‐species models at extremely low sample sizes.
#' Ecography, 2023(6), e06500.
#'
#' @seealso [train()]
#' @keywords engine control
#'
#' @examples
#' \dontrun{
#'  x <- distribution(background) |>
#'    add_biodiversity_poipa(species) |>
#'    add_predictors(covariates) |>
#'    add_control_esm(n_covs = 2)
#' }
#'
#' @name add_control_esm
NULL

#' @rdname add_control_esm
#' @export
methods::setGeneric(
  "add_control_esm",
  signature = methods::signature("x"),
  function(x, n_covs = 2) standardGeneric("add_control_esm"))

#' @rdname add_control_esm
methods::setMethod(
  "add_control_esm",
  methods::signature(x = "BiodiversityDistribution"),
  function(x, n_covs = 2) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.numeric(n_covs),
                            n_covs > 0
    )

    # Make a clone
    y <- x$clone(deep = TRUE)

    # Set the control object
    y <- y$set_control(type = "train", method = "esm",
                       value = c("n_covs" = n_covs))

    # Return the cloned object
    return(y)
  }
)
