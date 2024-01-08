#' Small helper function to obtain predictions from an object
#'
#' @description This function is a short helper function to return the fitted
#' data from a `DistributionModel` or `BiodiversityScenario` object. It can be
#' used to easily obtain for example the estimated prediction from a model or
#' the projected scenario from a `scenario()` object.
#'
#' @note This function is essentially identical to querying the internal
#' function \code{x$get_data()} from the object. However it does attempt some
#' lazy character matching if what is supplied.
#'
#' @param obj Provided [`DistributionModel`] or [`BiodiversityScenario`] object.
#' @param what A [`character`] of specific layer to be returned if existing
#'   (Default: \code{NULL}).
#' @keywords utils
#' @examples
#' \dontrun{
#'  # Assumes previously computed model
#'  get_data(fit)
#' }
#' @returns A `SpatRaster` or "stars" object depending on the input.
#' @export
methods::setGeneric("get_data",
                    signature = methods::signature("obj"),
                    function(obj, what = NULL) standardGeneric("get_data"))

#' @name get_data
#' @rdname get_data
methods::setMethod(
  "get_data",
  methods::signature("ANY"),
  function(obj, what = NULL){
    assertthat::assert_that(
      inherits(obj, "DistributionModel") || inherits(obj, "BiodiversityScenario"),
      (is.null(what) || is.character(what))
    )

    # Check for waiver and return respectively
    if(is.Waiver(obj$get_data())){
      return( new_waiver() )
    }

    # Depending on the object return the data
    if(inherits(obj, "DistributionModel")){
      if(is.null(what)){
        what <- "prediction"
      } else {
        vals <- obj$show_rasters()
        what <- grep(what, vals, value = TRUE,ignore.case = TRUE)
        if(length(what)==0){
          myLog("Modification", "red", "Specific layer not found. Returning whole object.")
          what <- NULL
        }
      }
      results <- obj$get_data(x = what)
      if(is.Waiver(results)) stop("No fitted prediction in object.")
    } else if(inherits(obj, "BiodiversityScenario")){
      results <-obj$get_data()
      if(!is.null(what)){
        vals <- names(results)
        what <- grep(what, vals, value = TRUE,ignore.case = TRUE)
        if(length(what)==0){
          myLog("Modification", "red", "Specific layer not found. Returning whole object.")
          what <- NULL
        } else {
          results <- results |> dplyr::select(what)
        }
      }
    }
    # Quick assertion
    assertthat::assert_that(is.Raster(results) || inherits(results, "stars"),
                            msg = "Something went wrong with the output?")
    return(results)
  }
)
