#' Built formula for glmnet model
#'
#' @description
#' This function builds a formula for a `engine_glmnet()` model.
#' @param obj A [`list()`] object containing the prepared model data for a given biodiversity dataset.
#' @note Function is not meant to be run outside the train() call.
#' @author Martin Jung
#' @keywords internal
#' @noRd
built_formula_glmnet <- function(obj){
  assertthat::assert_that(
    is.list(obj),
    length(obj) > 0,
    assertthat::has_name(obj, "observations"),
    assertthat::has_name(obj, "equation"),
    assertthat::has_name(obj, "predictors_names"),
    msg = "Error in model object. This function is not meant to be called outside ouf train()."
  )

  # Default equation found
  if(obj$equation =='<Default>' || is.Waiver(obj$equation)){
    # Construct formula with all variables
    form <- paste( 'observed', ifelse(obj$family=='poisson', '/w', ''), ' ~ ')
    # Add linear predictors
    form <- paste(form, paste0(obj$predictors_names, collapse = ' + '))
    # Convert to formula
    form <- to_formula(form)
  } else{
    # FIXME: Also make checks for correctness in supplied formula, e.g. if variable is contained within object
    if(getOption('ibis.setupmessages')) myLog('[Estimation]','yellow','Use custom model equation')
    form <- to_formula(obj$equation)
    assertthat::assert_that(
      all( all.vars(form) %in% c('observed', obj[['predictors_names']]) )
    )
  }

  return(form)
}

#' Determine best lambda
#'
#' @description
#' For glmnet fits the variables are often overregularized. This helper function picks
#' the best lambda estimate from the model.
#' By default use the one within 1 SE of minimum lambda, unless it falls on the very first value,
#' likely indicating an overregularized model. In this case take the minimum value
#' of all lambda's.
#' @param obj A \code{"glmnet"} object.
#' @keywords internal, utils
#' @noRd
determine_lambda <- function(obj){
  assertthat::assert_that(
    inherits(obj, "cv.glmnet"),
    is.numeric(obj$lambda.1se),
    is.numeric(obj$lambda.min),
    is.numeric(obj$lambda)
  )
  if(obj$lambda.1se != obj$lambda[1])	{
    la <- obj$lambda.1se
  } else if (obj$lambda.min != obj$lambda[1]){
    la <- obj$lambda.min
  } else {
    la <- tail(obj$lambda ,1)
  }
  return(la)
}

#' Summarize cross-validated glmnet model
#'
#' @description
#' This helper function summarizes the coefficients from a glmnet model.
#' The optimal lambda is determined through the [`determine_lambda`] function.
#' @param obj An object created with \code{'cv.glmnet'}.
#' @keywords internal, utils
#' @noRd
tidy_glmnet_summary <- function(obj){
  assertthat::assert_that(
    inherits(obj, "cv.glmnet")
  )
  # Determine best lambda
  lambda <- determine_lambda(obj)

  # Summarise coefficients within 1 standard deviation
  ms <- coef(obj, s = lambda) |>
    as.matrix() |> as.data.frame()
  names(ms) <- "mean"
  ms$variable <- rownames(ms)
  ms <- ms[,c("variable", "mean")]
  ms <- subset(ms, mean != 0) # Remove regularized coefficients for some clean up.
  if(nrow(ms)>0){
    # Reorder
    ms <- ms[order(ms$mean,decreasing = TRUE),] # Sort
    rownames(ms) <- NULL
  } else {
    ms <- data.frame()
  }
  return(ms)
}
