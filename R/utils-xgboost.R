#' Built formula for XGBOOST model
#'
#' @description This function built a formula for a `engine_xgboost()` model.
#'
#' @param obj A [`list()`] object containing the prepared model data for a given
#' biodiversity dataset.
#'
#' @note Function is not meant to be run outside the train() call.
#'
#' @author Martin Jung
#'
#' @noRd
#'
#' @keywords internal
built_formula_xgboost <- function(obj){
  assertthat::assert_that(
    is.list(obj),
    length(obj) > 0,
    assertthat::has_name(obj, "observations"),
    assertthat::has_name(obj, "equation"),
    assertthat::has_name(obj, "predictors_names"),
    msg = "Error in model object. This function is not meant to be called outside of train()."
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
    if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','yellow','Use custom model equation')
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
#'
#' @keywords utils
#'
#' @noRd
#'
#' @keywords internal
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

#' Format XGBoost priors for xgb.train
#'
#' @param priors A [`PriorList`] object or Waiver.
#' @param feature_names Feature names after XGBoost preprocessing.
#' @param booster XGBoost booster type.
#'
#' @keywords internal
#' @noRd
format_xgboost_priors <- function(priors, feature_names, booster = "gbtree"){
  assertthat::assert_that(
    is.Waiver(priors) || inherits(priors, "PriorList"),
    is.character(feature_names),
    length(feature_names) > 0,
    is.character(booster),
    length(booster) == 1
  )

  out <- list()
  if(is.Waiver(priors) || priors$length() == 0) return(out)

  prior_objects <- priors$priors
  prior_classes <- vapply(prior_objects, function(x) x$get_name(), character(1))
  monotone_priors <- prior_objects[prior_classes == "XGBPrior"]
  interaction_priors <- prior_objects[prior_classes == "XGBInteractionPrior"]

  if(booster == "gblinear"){
    if(length(monotone_priors) > 0 || length(interaction_priors) > 0){
      warning(
        "XGBoost monotone and interaction constraints are only supported for tree boosters. Ignoring xgboost priors for gblinear.",
        call. = FALSE
      )
    }
    return(out)
  }

  if(length(monotone_priors) > 0){
    monotone_variables <- vapply(monotone_priors, function(x) x$variable, character(1))
    missing_monotone <- monotone_variables[!monotone_variables %in% feature_names]
    assertthat::assert_that(
      length(missing_monotone) == 0,
      msg = paste0("XGBoost monotone prior variables missing from model features: ",
                   paste(missing_monotone, collapse = ", "))
    )

    mc <- rep(0, length(feature_names))
    names(mc) <- feature_names
    for(prior in monotone_priors){
      mc[prior$variable] <- switch(prior$get(),
                                   'increasing' = 1, 'positive' = 1,
                                   'decreasing' = -1, 'negative' = -1,
                                   0
      )
    }
    out$monotone_constraints <- mc
  }

  if(length(interaction_priors) > 0){
    interaction_groups <- unname(lapply(interaction_priors, function(x) x$get()))
    interaction_variables <- unlist(interaction_groups, use.names = FALSE)
    missing_interactions <- unique(interaction_variables[
      !interaction_variables %in% feature_names
    ])
    assertthat::assert_that(
      length(missing_interactions) == 0,
      msg = paste0("XGBoost interaction prior variables missing from model features: ",
                   paste(missing_interactions, collapse = ", "))
    )

    out$interaction_constraints <- unname(lapply(interaction_groups, function(group){
      as.integer(match(group, feature_names) - 1L)
    }))
  }

  return(out)
}
