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
                    function(background, ...) standardGeneric("distribution"))

#' @name distribution
#' @usage \S4method{distribution}{Raster}(background, formula, ...)
#' @rdname distribution
methods::setMethod(
  "distribution",
  methods::signature(background = "Raster"),
  function(background, ...) {
    # Check that arguments are valid
    assertthat::assert_that(
      inherits(background,'Raster'),
      raster::nlayers(background) == 1
    )

    # Create BiodiversityDistribution object
    bdproto(NULL, BiodiversityDistribution,
            background = background,
            ...
            )
  })
