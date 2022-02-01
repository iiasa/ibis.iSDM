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
        burn = BoomSpikeSlab::SuggestBurn(obj), # ceiling(params$iter*0.1)
        type = params$type,
        mean.only = FALSE # Return full posterior
      )
    )
  } else if(fam == "binomial"){
    suppressWarnings(
      pred_breg <- BoomSpikeSlab::predict.logit.spike(
        object = obj,
        newdata = newdata,
        burn = BoomSpikeSlab::SuggestBurn(obj),
        type = params$type,
        mean.only = FALSE # Return full posterior
      )
    )
  } else {
    suppressWarnings(
      pred_breg <- BoomSpikeSlab::predict.lm.spike(
        object = obj,
        newdata = newdata,
        burn = BoomSpikeSlab::SuggestBurn(obj),
        type = params$type,
        mean.only = FALSE # Return full posterior
      )
    )
  }
  return(pred_breg)
}

#' Spline predictors for non-linearity
#'
#' @description This function takes an existing [data.frame]
#' of predictors and transforms them through B- and I-splines (if specified).
#'
#' @param df A [`data.frame`] containing the original predictors.
#' @param transform A [`vector`] specifying the type of transformation to be applied.
#' Either \code{B} and \code{I} splines are available (Default: \code{c("B", "I")}).
#' @returns
#' A [list] with 2 entries:
#' * A [`data.frame`] of the original df and transformed data
#' * A [`vector`] with the new predictor names
#' @keywords utils, internal
#' @noRd
predictor_spline <- function( df, transform = c("B", "I")){
  assertthat::assert_that(
    is.data.frame(df), nrow(df) > 0,
    is.vector(transform)
  )
  # Match transformation parameter
  transform <- match.arg(transform, c("B", "I"),several.ok = TRUE)

  # Define output list
  out <- list()
  out[["predictors_names"]] <- names(df)

  # B-splies
  if("B" %in% transform){
    BoomSpikeSlab::BsplineBasis(df$Girth)
  }
  # I-splines (monotonic)
  if("B" %in% transform){
    BoomSpikeSlab::IsplineBasis(df$Girth)
  }

}
