#' Setup a prior for `Boom` engine model
#'
#' @description
#' For more information see helpfile of [`BREGPrior`] and the respective package
#' helpfiles
#' @param form A [`formula`] object.
#' @param data A [`data.frame`] with all variables (response and predictors) in the formula.
#' Needs to contain the observed y_hat variable as well.
#' @param priors A [`PriorList`] object with [`BREGPrior`] priors.
#' @param family A [`character`] object giving either `poisson` or `binomial`.
#' @param exposure A [`numeric`] vector giving the exposure for `poisson` family priors.
#' @returns A [`SpikeSlabPriorBase`] object for use with a [`Boom`] engine trained model
#' @keywords utils, internal
#' @noRd
setup_prior_boom <- function(form, data, priors, family, exposure = NULL){
  assertthat::assert_that(
    is.formula(form),
    is.data.frame(data),
    inherits(priors, "PriorList"),
    is.null(exposure) || is.numeric(exposure),
    is.character(family)
  )
  family <- match.arg(family, c("poisson", "binomial"), several.ok = FALSE)

  # Check that all terms are in the data.frame
  assertthat::assert_that(
    all( attr(terms.formula(form),"term.labels") %in% names(data) ),
    "observed" %in% names(data)
  )
  # Create model matrix
  mm <- model.matrix.default(object = form, data = data)

  # Expected model size:
  # Conservatively just assume the priors are relevant (Default is 1)
  esize <- priors$length()

  # Optional coefficient estimates
  if(any(is.numeric( priors$get(priors$varnames()) ))){
    # If any are set, define optional coefficient estimates
    co <- vector(length = ncol(mm));co[] <- 0;names(co) <- colnames(mm)
    co["(Intercept)"] <- mean( mean(data[["observed"]]) )
    for(val in priors$varnames()){
      co[val] <- priors$get(val, what = "value")
    }
  } else { co <- NULL}

  # Option inclusion probabilities
  if(any(is.numeric( priors$get(priors$varnames(), what = "prob") ))){
    # Probabily of inclusion for each variable
    if(esize < (ncol(mm)) ) {
      #Specify default as each having equal probability
      co.ip <- rep(esize/ncol(mm), ncol(mm))
    } else {
      co.ip <- rep(1, ncol(mm))
    }
    names(co.ip) <- colnames(mm)
    # Now set priors for those where set
    for(val in priors$varnames()){
      co.ip[val] <- priors$get(val, what = "prob")
    }
  } else { co.ip <- NULL}

  # Create priors depending on input family
  if(family == "poisson"){
    assertthat::assert_that(!is.null(exposure))
    pp <- PoissonZellnerPrior(predictors = mm,
                              counts = data$observed,
                              exposure = exposure,
                              expected.model.size = esize,
                              optional.coefficient.estimate = co,
                              prior.inclusion.probabilities = co.ip
    )
  } else {
    # For binomial
    pp <- LogitZellnerPrior(predictors = mm,
                            successes = data$observed,
                            trials = exposure,
                            expected.model.size = esize,
                            optional.coefficient.estimate = co,
                            prior.inclusion.probabilities = co.ip
    )
  }
  # Return the created prior object
  assertthat::assert_that(inherits(pp, "SpikeSlabPriorBase"))
  return(pp)
}

#' Prediction with `Boom` package for breg models
#'
#' @note By Default 20% of the iterations are considered as burnin.
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
