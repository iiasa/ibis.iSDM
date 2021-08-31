#' @include utils.R bdproto-biodiversityscenario.R
NULL

#' Create a new scenario based on trained model parameters
#'
#' @param fit A trained [`BiodiversityDistribution`] model
#'
#' @aliases scenario
#' @exportMethod scenario
#' @name scenario
#'
#' @examples
#' \dontrun{
#' print('test')
#' }
#' @export
methods::setGeneric("scenario",
                    signature = methods::signature("fit"),
                    function(fit) standardGeneric("scenario"))

#' @name scenario
#' @usage \S4method{scenario}{ANY}(fit)
#' @rdname scenario
methods::setMethod(
  "scenario",
  methods::signature(fit = "ANY"),
  function(fit) {
    # Check that arguments are valid
    assertthat::assert_that(!missing(fit) || inherits(fit,'DistributionModel'),
                            msg = 'No trained model supplied!')

    # Get model object name and id
    modelobject <- deparse(substitute(fit))
    modelid <- fit$id

    # Create BiodiversityScenario object
    bdproto(NULL, BiodiversityScenario,
            modelobject = modelobject,
            modelid = modelid

    )
  })
