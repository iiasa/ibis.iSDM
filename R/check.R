#' @include utils.R
NULL

#' Check objects in the package for common errors or issues
#'
#' @description
#' Not always is there enough data or sufficient information to robustly
#' infer the suitable habitat or niche of a species. As many SDM algorithms are
#' essentially regression models, similar assumptions about model convergence,
#' homogeneity of residuals and inferrence usually apply (although often
#' ignored).
#' This function simply checks the respective input object for common issues or
#' mistakes.
#'
#' @details
#' Different checks are implemented depending on the supplied object
#'
#' * [`BiodiversityDistribution`]
#' -
#'
#' * [`DistributionModel`]
#' -
#'
#' * [`BiodiversityScenario`]
#' -
#'
#' @note
#' This function will likely be expanded with additional checks in the future.
#' If you have ideas, please let them know per issue.
#'
#' @param obj A [`BiodiversityDistribution`], [`DistributionModel`] or [`BiodiversityScenario`] object.
#' @param stoponwarning [`logical`] Should check return a stop if warning is raised? (Default: \code{FALSE}).
#' @name check
#' @returns Message outputs
#' @keywords misc
#' @examples
#' \dontrun{
#'  # Where mod is an estimated DistributionModel
#'  check(mod)
#' }
#' @export
NULL

#' @name check
#' @rdname check
#' @exportMethod check
#' @export
methods::setGeneric(
  "check",
  signature = methods::signature("obj", "stoponwarning"),
  function(obj, stoponwarning = FALSE) standardGeneric("check"))

#' @name check
#' @rdname check
#' @usage \S4method{check}{ANY, logical}(obj, stoponwarning)
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

      # Are there enough observations?
      if( sum( obj$biodiversity$get_observations() ) ){
        ms$Warnings <- append(ms$Warnings, "Not enough observations found.")
      }

    }
    if(inherits(obj, "DistributionModel")){

    }
    if(inherits(obj, "BiodiversityScenario")){

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
      myLog("[Checks]",col = "green", entry)
    }
  }
)
