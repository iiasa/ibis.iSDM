#' Prediction with `Boom` package for breg models
#'
#' @param obj A [list] containing the fitted model
#' @param newdata A [`data.frame`] with all the predictor used for model fitting
#' @param fam A [`character`] denoting the family used for estimation. This
#' is necessary as prediction methods differ among the various functions.
#' @param params A [`list`] with paramters for estimation. Normally created during
#' model fitting.
#' @param w A [`numeric`] [`vector`] containing the exposure variables for PPMs. Can
#' be \code{NULL} if the model is not a PPM.
#' @returns A [`data.frame`] with the respective prediction
#' @keywords utils, internal
#' @noRd
predict_boom <- function(obj, newdata, fam, params, w = NULL) {
  assertthat::assert_that(
    is.list(obj),
    is.data.frame(newdata) || inherits(newdata, "SpatialPixelsDataFrame"),
    is.character(fam),
    is.list(params),
    is.null(w) || is.numeric(w)
  )
  check_package("BoomSpikeSlab")

  # Make a prediction
  if(fam == "poisson"){
    suppressWarnings(
      pred_breg <- BoomSpikeSlab::predict.poisson.spike(
        object = obj,
        newdata = newdata,
        exposure = w,
        burn = ceiling(params$iter*0.2),
        type = params$type,
        mean.only = FALSE # Return full posterior
      )
    )
  } else if(fam == "binomial"){
    suppressWarnings(
      pred_breg <- BoomSpikeSlab::predict.logit.spike(
        object = obj,
        newdata = newdata,
        burn = ceiling(params$iter*0.2),
        type = params$type,
        mean.only = FALSE # Return full posterior
      )
    )
  } else {
    suppressWarnings(
      pred_breg <- BoomSpikeSlab::predict.lm.spike(
        object = obj,
        newdata = newdata,
        burn = ceiling(params$iter*0.2),
        type = params$type,
        mean.only = FALSE # Return full posterior
      )
    )
  }
  return(pred_breg)
}
