#' @include class-biodiversitydistribution.R class-predictors.R class-biodiversityscenario.R
NULL

#' Function to add GLOBIOM-DownScalr derived predictors to a Biodiversity distribution
#' object
#'
#' @description
#' **This function is defunct! Use the BNRTools package for formatting input data!**
#'
#' @inheritParams add_predictors
#' @param ... Other parameters passed down
#'
#' @details See [`add_predictors()`]
#'
#' @seealso [add_predictors]
#' @references [https://github.com/iiasa/BNRTools](https://github.com/iiasa/BNRTools)
#'
#' @examples
#' \dontrun{
#'  obj <- distribution(background) |>
#'         add_predictors_globiom(fname = "", transform = 'none')
#'  obj
#' }
#'
#' @name add_predictors_globiom
NULL

#' @rdname add_predictors_globiom
#' @export
methods::setGeneric(
  "add_predictors_globiom",
  signature = methods::signature("x"),
  function(x, ...) standardGeneric("add_predictors_globiom"))

#' @rdname add_predictors_globiom
methods::setMethod(
  "add_predictors_globiom",
  methods::signature(x = "BiodiversityDistribution"),
  function(x, ... ) {

    cli::cli_abort(c("add_predictors_globiom() is defunct.",
                     "v" = "Use the BNRTools package directly!"), call = NULL)

  }
)

#' @rdname add_predictors_globiom
methods::setMethod(
  "add_predictors_globiom",
  methods::signature(x = "BiodiversityScenario"),
  function(x, ... ) {

    cli::cli_abort(c("add_predictors_globiom() is defunct.",
                     "v" = "Use the BNRTools package directly!"), call = NULL)
  }
)
