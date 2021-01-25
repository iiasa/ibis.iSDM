#' @include utils.R bdproto-biodiversitydistribution.R
NULL

#' Setup biodiversity distribution modelling procedure
#'
#' @param background Specification of the modelling background. Must be a
#' [`raster`], [`sf`] or [`extent`] object
#' @param formula Model formula. Has to be [`character`] or [`formula`] object (optional)
#' @param ... not used.
#'
#' @details TBD. Say something about PPMs, INLA and co
#'
#' @return [`BiodiversityDistribution-class`] object containing
#'   data for building a biodiversity distribution modelling problem.
#'
#' @seealso [bdproto]
#'
#' @aliases distribution
#'
#' @exportMethod distribution
#'
#' @name distribution
#'
#' @examples
#' \dontrun{
#' print('test')
#' }
#' @export
methods::setGeneric("distribution",
                    signature = methods::signature("background"),
                    function(background, formula = NULL, ...) standardGeneric("distribution"))

#' @name distribution
#' @usage \S4method{distribution}{Raster}(background, formula, ...)
#' @rdname distribution
methods::setMethod(
  "distribution",
  methods::signature(background = "Raster"),
  function(background, formula = NULL, ...) {
    # Check that arguments are valid
    assertthat::assert_that(
      inherits(background,'Raster'),
      inherits(formula,'formula') || is.null(formula) || is.character(formula),
      raster::nlayers(background) == 1
    )

    # Convert to formula object
    if(!is.null(formula)) {
        formula = as.formula(formula)
    } else {
      # Asign a new waiver object
        formula = new_waiver()
    }

    # Create BiodiversityDistribution object
    bdproto(NULL, BiodiversityDistribution,
            background = background,
            equation = formula
            )
  })
