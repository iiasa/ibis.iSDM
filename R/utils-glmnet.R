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
    if(getOption('ibis.setupmessages')) myLog('[Estimation]','yellow','Use custom model equation')
    form <- to_formula(obj$equation)
    # If response is missing, add manually
    if(attr(stats::terms(form), "response")==0){
      form <- stats::update.formula(form, "observed ~ .")
    }
    # Security checks
    assertthat::assert_that(
      is.formula(form),
      attr(stats::terms(form), "response")==1, # Has Response
      all( all.vars(form) %in% c('observed', obj[['predictors_names']]) )
    )
  }

  return(form)
}

#' Default regularization constant
#'
#' @description
#' This function was taken from the [`maxnet`] R-package to get some more informed
#' default lambda values for the regularization.
#'
#' @param p A [`vector`] of \code{1} (for presence) or \code{0} (for background).
#' @param m A [`model.matrix`] object.
#' @references
#' * Phillips S (2021). _maxnet: Fitting 'Maxent' Species Distribution Models with 'glmnet'_. R package version 0.1.4, <https://CRAN.R-project.org/package=maxnet>.
#' @source maxnet
#' @keywords internal, utils
#' @noRd
default.regularization <- function(p, m){
  isproduct <- function(x) grepl(":", x) & !grepl("\\(", x)
  isquadratic <- function(x) grepl("^I\\(.*\\^2\\)", x)
  ishinge <- function(x) grepl("^hinge\\(", x)
  isthreshold <- function(x) grepl("^thresholds\\(", x)
  iscategorical <- function(x) grepl("^categorical\\(", x)
  regtable <- function(name, default) {
    if (ishinge(name))
      return(list(c(0, 1), c(0.5, 0.5)))
    if (iscategorical(name))
      return(list(c(0, 10, 17), c(0.65, 0.5, 0.25)))
    if (isthreshold(name))
      return(list(c(0, 100), c(2, 1)))
    default
  }
  lregtable <- list(c(0, 10, 30, 100), c(1, 1, 0.2, 0.05))
  qregtable <- list(c(0, 10, 17, 30, 100), c(1.3, 0.8, 0.5,
                                             0.25, 0.05))
  pregtable <- list(c(0, 10, 17, 30, 100), c(2.6, 1.6, 0.9,
                                             0.55, 0.05))
  mm <- m[p == 1, ]
  np <- nrow(mm)
  lqpreg <- lregtable
  if (sum(isquadratic(colnames(mm))))
    lqpreg <- qregtable
  if (sum(isproduct(colnames(mm))))
    lqpreg <- pregtable
  classregularization <- sapply(colnames(mm), function(n) {
    t <- regtable(n, lqpreg)
    stats::approx(t[[1]], t[[2]], np, rule = 2)$y
  })/sqrt(np)
  ishinge <- grepl("^hinge\\(", colnames(mm))
  hmindev <- sapply(1:ncol(mm), function(i) {
    if (!ishinge[i])
      return(0)
    avg <- mean(mm[, i])
    std <- max(stats::sd(mm[, i]), 1/sqrt(np))
    std * 0.5/sqrt(np)
  })
  tmindev <- sapply(1:ncol(mm), function(i) {
    ifelse(isthreshold(colnames(mm)[i]) && (sum(mm[, i]) ==
                                              0 || sum(mm[, i]) == nrow(mm)), 1, 0)
  })
  pmax(0.001 * (apply(m, 2, max) - apply(m, 2, min)), hmindev,
       tmindev, apply(as.matrix(mm), 2, sd) * classregularization)
}

#' Determine best lambda
#'
#' @description
#' For glmnet fits the variables are often overregularized. This helper function picks
#' the best lambda estimate from the model.
#' By default use the one within 1 SE of minimum lambda, unless it falls on the very first value,
#' likely indicating an overregularized model. In this case take the minimum value of all lambda's.
#' @param obj A \code{"glmnet"} object.
#' @aliases determine_lambda
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
  ms <- stats::coef(obj, s = lambda) |>
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
