#' Built formula for XGBOOST model
#'
#' @description This function built a formula for a `engine_xgboost()` model.
#' @param obj A [`list()`] object containing the prepared model data for a given
#'   biodiversity dataset.
#' @author Martin Jung
#' @note Function is not meant to be run outside the train() call.
#' @keywords internal
#' @noRd
built_formula_xgboost <- function(obj){
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
    # XGboost does not explicitly work with formulas, thus all supplied objects
    # are assumed to be part a covariate
    form <- new_waiver()
    # Note: Priors are added in the fitted distribution object through the model
    # object
  } else{
    # If custom supplied formula, check that variable names match the supplied
    # predictors
    if(getOption('ibis.setupmessages')) myLog('[Estimation]','yellow','Use custom model equation')
    form <- to_formula(obj$equation)

    # Get all variables and check
    varn <- obj$predictors_names[which( obj$predictors_names %in% formula_terms(form))]
    form <- to_formula( paste0("observed ~", paste0(varn, collapse = " + ")) )
    assertthat::assert_that(
      is.formula(form), length(all.vars(form))>1
    )
  }

  return(form)
}

#' Split up a factor variable into binary components
#'
#' @param df A [`vector`] object containing the factor variables
#' @param name Name for the new object
#' @keywords utils, internal
#' @noRd
explode_factor <- function(df, name = "facvar"){
  assertthat::assert_that(
    is.data.frame(df) || is.factor(df),
    all(is.factor(df)),
    is.character(name)
  )
  z <- as.data.frame(
    outer(df, levels(df), function(w, f) ifelse(w == f, 1, 0))
    )
  names(z) <- paste(name, levels(df), sep = ".")
  return(z)
}
