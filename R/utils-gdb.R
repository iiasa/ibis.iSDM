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

#' Calculate weights for Point Process models
#'
#' @param df The [`data.frame`] for which weights are to be calculated
#' @param presence A [`vector`] with the observed species. Has to be in range 0 to Inf
#' @param bg A background [`raster`] layer
#' @param weight A [`numeric`] with the base weight. Default to 1e-6
#' @references Renner, I.W., Elith, J., Baddeley, A., Fithian, W., Hastie, T., Phillips, S.J., Popovic, G. and Warton, D.I., 2015. Point process models for presence‐only analysis. Methods in Ecology and Evolution, 6(4), pp.366-379.
#' @return A vector with the weights
#' @keywords utils
#' @noRd
ppm_weights <- function(df, pa, bg, weight = 1e-4){
  assertthat::assert_that(
    is.data.frame(df),
    nrow(df) == length(pa),
    is.numeric(weight)
  )
  # number of non-NA cells
  nc = cellStats(!is.na(bg), sum)

  # Set output weight as default
  w <- rep( weight, nrow(df) )
  w[which(pa == 0)] = nc / sum(pa == 0)

  assertthat::assert_that(
    length(unique(w)) > 1,
    length(w) == nrow(df)
  )
  return(w)
}
