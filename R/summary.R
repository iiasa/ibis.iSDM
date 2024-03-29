#' Summarises a trained model or predictor object
#'
#' @description This helper function summarizes a given object, including
#' [DistributionModel], [PredictorDataset] or [PriorList] objects and others.
#' This can be a helpful way to summarize what is contained within and the
#' values of specified models or objects.
#'
#' When unsure, it is usually a good strategy to run [summary] on any object.
#'
#' @param object Any prepared object.
#' @param ... not used.
#'
#' @seealso [base::summary()].
#' @keywords summary
#'
#' @examples
#' \dontrun{
#' # Example with a trained model
#' x <- distribution(background) |>
#'         # Presence-absence data
#'         add_biodiversity_poipa(surveydata) |>
#'         # Add predictors and scale them
#'         add_predictors(env = predictors) |>
#'         # Use glmnet and lasso regression for estimation
#'         engine_glmnet(alpha = 1)
#'  # Train the model
#'  mod <- train(x)
#'  summary(mod)
#'
#'  # Example with a prior object
#'  p1 <- BREGPrior(variable = "forest", hyper = 2, ip = NULL)
#'  p2 <- BREGPrior(variable = "cropland", hyper = NULL, ip = 1)
#'  pp <- priors(p1,p2)
#'  summary(pp)
#' }
#'
#'
#' @name summary
NULL

#' @rdname summary
#' @method summary distribution
#' @export
summary.distribution <- function(object, ...) object$summary()

#' @rdname summary
#' @method summary DistributionModel
#' @export
summary.DistributionModel <- function(object, ...) object$summary()

#' @rdname summary
#' @method summary PredictorDataset
#' @export
summary.PredictorDataset <- function(object, ...) object$summary()

#' @rdname summary
#' @method summary BiodiversityScenario
#' @export
summary.BiodiversityScenario <- function(object, ...) object$summary()

#' @rdname summary
#' @method summary PriorList
#' @export
summary.PriorList <- function(object, ...) object$summary()

#' @rdname summary
#' @method summary Settings
#' @export
summary.Settings <- function(object, ...) object$summary()

#' Obtains the coefficients of a trained model
#'
#' @description Similar as [`summary`], this helper function obtains the
#' coefficients from a given [DistributionModel] object.
#'
#' @param object Any prepared object.
#' @param ... not used.
#'
#' @note For models trained with machine-learning approaches (e.g. [`engine_bart`]
#' etc) this function will return variable importance estimates rather than
#' linear coefficients. Similar can be said for trained non-linear models.
#'
#' @seealso [stats::coef()].
#' @keywords coef
#'
#' @name coef
NULL

#' @rdname coef
#' @export
coef.DistributionModel <- function(object, ...) object$get_coefficients()
