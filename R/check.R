#' @include utils.R
NULL

#' Check objects in the package for common errors or issues
#'
#' @description Not always is there enough data or sufficient information to
#' robustly infer the suitable habitat or niche of a species. As many SDM
#' algorithms are essentially regression models, similar assumptions about model
#' convergence, homogeneity of residuals and inferrence usually apply (although
#' often ignored). This function simply checks the respective input object for
#' common issues or mistakes.
#'
#' @param obj A [`BiodiversityDistribution`], [`DistributionModel`] or
#' [`BiodiversityScenario`] object.
#' @param stoponwarning [`logical`] Should check return a stop if warning is raised?
#' (Default: \code{FALSE}).
#'
#' @details Different checks are implemented depending on the supplied object
#'
#' * [`BiodiversityDistribution`]
#' - Checks if there are less than 200 observations
#' - TODO: Add rm_insufficient_covs link
#'
#' * [`DistributionModel`]
#' - Check model convergence
#' - Check if model is found
#' - Check if coefficients exist
#' - Check if there are unusal outliers in prediction (using 10median absolute deviation)
#' - Check if threshold is larger than layer
#'
#' * [`BiodiversityScenario`]
#' -
#'
#' @note This function will likely be expanded with additional checks in the
#' future. If you have ideas, please let them know per issue.
#'
#' @returns Message outputs
#'
#' @keywords misc
#'
#' @examples
#' \dontrun{
#'  # Where mod is an estimated DistributionModel
#'  check(mod)
#' }
#'
#' @name check
NULL

#' @rdname check
#' @export
methods::setGeneric(
  "check",
  signature = methods::signature("obj", "stoponwarning"),
  function(obj, stoponwarning = FALSE) standardGeneric("check"))

#' @rdname check
methods::setMethod(
  "check",
  methods::signature(obj = "ANY"),
  function(obj, stoponwarning = FALSE) {
    assertthat::assert_that(
      is.logical(stoponwarning)
    )

    # Messages
    ms <- list(Warnings = character())

    # Type specific checks
    if(inherits(obj,"BiodiversityDistribution")){

      # Check if necessary data is there
      if(is.Waiver(obj$get_biodiversity_types())){
        ms$Warnings <- append(ms$Warnings, "No biodiversity data added.")
      } else {
        # Are there enough observations?
        if( sum( obj$biodiversity$get_observations() ) < 100 ){
          ms$Warnings <- append(ms$Warnings, "Less that 100 observations found.")
        }

        # Count number of observations and predictors
        t1 <- obj$biodiversity$get_observations() |> as.numeric()
        t2 <- length( obj$get_predictor_names() )
        if(!is.Waiver(t2)){
          if(t2 >= t1){
            ms$Warnings <- append(ms$Warnings,
                                  "More predictors than observations.")
          }
        }
      }

      # Check predictors
      if(is.Waiver(obj$get_predictor_names())){
        ms$Warnings <- append(ms$Warnings, "No predictors added.")
      }

    }
    if(inherits(obj, "DistributionModel")){

      fit <- obj$get_data("fit_best")
      # Did the model exist?
      if(is.Waiver(fit)){
        ms$Warnings <- append(ms$Warnings, "No fitted model found!")
      }

      # Sufficient iterations for converged?
      if(!obj$has_converged()){
        ms$Warnings <- append(ms$Warnings, "Model likely has not converged!")
      }

      # No coefficients?
      if(nrow(obj$get_coefficients())==0){
        ms$Warnings <- append(ms$Warnings, "No coefficients in the model!")
      }

      # Check if some coefficients are clear outliers
      if(nrow(obj$get_coefficients())>0){
        co <- obj$get_coefficients()
        # Remove large outlier from coefficient and check if that results in NA (outlier)
        co <- rm_outlier_revjack(co[[2]], procedure = "missing")
        if(anyNA(co)){
          ms$Warnings <- append(ms$Warnings, "Likely unstable coefficient in model (outlier)!")
        }
      }

      # Check if threshold exists, if so if it is equal to maximum
      if(!is.Waiver(obj$get_thresholdvalue())){
        # Get prediction
        pred <- obj$get_data()
        if(obj$get_thresholdvalue() >= terra::global(pred[[1]],"max",na.rm=TRUE)[,1]){
          ms$Warnings <- append(ms$Warnings, "Threshold larger than prediction!")
        }
      }

      # Check for positive outliers in Prediction
      if(!is.Waiver( obj$get_data() )){
        pred <- obj$get_data()[['mean']]
        # Calculate outliers using the mad
        pmed <- terra::global(pred, median, na.rm = TRUE)[,1]
        abs_dev <- abs(pred[]-pmed)[,1]
        # Calculate mad
        pmad <- 1.4826 * median(abs_dev,na.rm = TRUE)
        test <- terra::values(pred)[,1] > (pmed + (10 * pmad)) # 10x
        # For a single outlier value, raise warning
        if(length(test)==1){
          ms$Warnings <- append(ms$Warnings, "Single outlier in prediction. Likely overfit!")
        }
      }
    }


    if(inherits(obj, "BiodiversityScenario")){
      invisible()
    }

    # Check for warnings
    if(stoponwarning){
      # Check if there are any warnings
      if(any(length(ms$Warnings)>0)){
        stop("Warning raised during checks!")
      }
    }
    # Compile
    for(entry in ms$Warnings){
      myLog("[Checks]",col = "yellow", entry)
    }
  }
)
