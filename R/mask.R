#' @include utils.R
NULL

#' Mask data with an external layer
#'
#' @description This is a helper function that takes an existing object created
#' by the ibis.iSDM package and an external layer, then intersects both. It
#' currently takes either a [DistributionModel], [BiodiversityDatasetCollection],
#' [PredictorDataset] or [BiodiversityScenario] as input.
#'
#' As mask either a [`sf`] or [`SpatRaster`] object can be chosen. The mask will
#' be converted internally depending on the object.
#'
#' @param x Any object belonging to [DistributionModel], [BiodiversityDatasetCollection],
#' [PredictorDataset] or [BiodiversityScenario].
#' @param mask A [`sf`] or [`SpatRaster`] object.
#' @param inverse A [`logical`] flag whether to take inverse of the mask instead
#' (Default: \code{FALSE}).
#' @param ... Passed on arguments
#'
#' @return A respective object of the input type.
#'
#' @seealso [terra::mask()]
#' @keywords utils
#'
#' @examples
#' \dontrun{
#' # Build and train a model
#' mod <- distribution(background) |>
#'   add_biodiversity_poipo(species) |>
#'   add_predictors(predictors) |>
#'   engine_glmnet() |>
#'   train()
#'
#' # Constrain the prediction by another object
#' mod <- mask(mod, speciesrange)
#'
#' }
#'
#' @name mask
NULL

#' @rdname mask
#' @method mask DistributionModel
#' @export
mask.DistributionModel <- function(x, mask, inverse = FALSE, ...) x$mask(mask,inverse,...)

#' @rdname mask
#' @method mask BiodiversityDatasetCollection
#' @export
mask.BiodiversityDatasetCollection <- function(x, mask, inverse = FALSE, ...) x$mask(mask,inverse,...)

#' @rdname mask
#' @method mask PredictorDataset
#' @export
mask.PredictorDataset <- function(x, mask, inverse = FALSE, ...) x$mask(mask,inverse,...)

#' @rdname mask
#' @method mask BiodiversityScenario
#' @export
mask.BiodiversityScenario <- function(x, mask, inverse = FALSE, ...) x$mask(mask,inverse,...)
