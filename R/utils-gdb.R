#' Use a fitted model for creating a new class prediction in raster form
#'
#' @param fit A fitted [`mboost`] model with [`binomial`] distribution
#' @param nd A new data.frame with all predictiors used in fit
#' @param template A [`Raster`] object that can be used as spatial template.
#' @returns A [`RasterLayer`] containing a presence-absence prediction.
#' @keywords utils
#' @noRd
predict_gdbclass <- function(fit, nd, template){
  assertthat::assert_that(
    inherits(fit, 'mboost'),
    !is.null( grep('Binomial', fit$family@name,ignore.case = TRUE) ),
    is.data.frame(nd),
    inherits(template,'RasterLayer')
  )
  # Redo a template to be sure
  template <- emptyraster(template)

  # Remove missing data in newdata data.frame
  nd$cellid <- rownames(nd)
  nd <- subset(nd, complete.cases(nd))

  suppressWarnings(
    pred_gdb <- mboost::predict.mboost(object = fit, newdata = nd,
                                       type = 'class', aggregate = 'sum')
  )
  # Fill output
  prediction <- emptyraster(template)
  prediction[as.numeric(nd$cellid)] <- pred_gdb
  prediction[prediction < max(as.numeric(pred_gdb),na.rm = T)] <- 0; prediction[prediction >0] <- 1
  names(prediction) <- 'presence'
  return(prediction)
}

#' Check training predictor complexity
#'
#' @description
#' This internal function tests an existing set of covariates (usually the training data)
#' on the number of unique values within.
#' If fewer values than a given threshold (\code{'tr'}) is detected, then the predictor is removed, thus
#' reducing complexity.
#' @note Maybe in the future a more cleverer solution could be thought of, for instance using a singular value decompoistion?
#' @param model A [`list`] of a model object containing the various predictors and biodiversity occurence information.
#' @param tr A [`numeric`] value describing a threshold of minimum unique values to be retained.
#' @returns A [`vector`] with the variables that fullfill the threshold.
#' @keywords internal
#' @noRd
rm_insufficient_covs <- function(model, tr = 5){
  assertthat::assert_that(
    is.list(model),
    is.numeric(tr),
    tr > 0
  )

  # Check that biodiversity information is present
  assertthat::assert_that(
    assertthat::has_name(model, "observations"),
    assertthat::has_name(model, "predictors_names")
  )

  # Now get all continuous ones
  vars_num <- model$predictors_types$predictors[model$predictors_types$type=="numeric"]
  vars_fac <- model$predictors_types$predictors[model$predictors_types$type=="factor"]

  vars_uniques <- apply(model$predictors[,vars_num], 2, function(x) length(unique(x,na.rm = TRUE)) )

  # Get all variables smaller than the threshold and return the original data.frame without them
  sufficient <- which(vars_uniques >= tr)

  # Get the factor variables in it as well
  if(length(vars_fac)>0){
    vars_uniques <- apply(model$predictors[,vars_fac], 2, function(x) length(unique(x,na.rm = TRUE)) )
    # Get all factor variables with at least 2 levels
    sufficient_fac <- which(vars_uniques >= 2)
    if(length(sufficient_fac)>0){
      sufficient <- c(names(sufficient), names(sufficient_fac))
    }
  } else {
    sufficient <- names(sufficient)
  }
  assertthat::assert_that(all(sufficient %in% model$predictors_names)) # This should return a character of covariate names
  if(length(sufficient)==0){
    return(NULL)
  } else {
    return(sufficient)
  }
}

#' Calculate weights for Point Process models
#'
#' @param df The [`data.frame`] for which weights are to be calculated.
#' @param presence A [`vector`] with the observed species. Has to be in range \code{0} to \code{Inf}
#' @param bg A background [`raster`] layer
#' @param weight A [`numeric`] weight to be used in down-weighted regressions.
#' @param type Accepting either “Infinitely weighted logistic regression” \code{'IWLR'} for use with binomial
#' logistic regressions or
#' “Down-weighted Poisson regression” \code{"DWPR"} (Default).
#' @references
#' * Renner, I.W., Elith, J., Baddeley, A., Fithian, W., Hastie, T., Phillips, S.J., Popovic, G. and Warton, D.I., 2015. Point process models for presence‐only analysis. Methods in Ecology and Evolution, 6(4), pp.366-379.
#' * Fithian, W. & Hastie, T. (2013) Finite-sample equivalence in statistical models for presence-only data. The Annals of Applied Statistics 7, 1917–1939
#' @return A vector with the weights
#' @keywords utils
#' @noRd
ppm_weights <- function(df, pa, bg, weight = 1e-6, type = "DWPR"){
  assertthat::assert_that(
    is.data.frame(df),
    length(unique(pa)) > 1,
    nrow(df) == length(pa),
    is.numeric(weight),
    is.character(type)
  )
  type <- match.arg(type, c("DWPR", "IWLR"),several.ok = FALSE)

  # number of non-NA cells
  nc = cellStats(!is.na(bg), sum)

  # Set output weight as default
  if(type == "DWPR"){
    w <- rep( weight, nrow(df) )
    w[which(pa == 0)] = nc / sum(pa == 0)
  } else {
    w = (10^6)^(1 - pa)
  }

  assertthat::assert_that(
    length(unique(w)) > 1,
    length(w) == nrow(df)
  )
  return(w)
}
