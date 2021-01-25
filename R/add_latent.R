#' @include utils.R bdproto.R waiver.R
NULL

#' Add latent spatial factor to the model
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param ... Other parameters passed down
#'
#' @details Different for INLA. Otherwise a CAR object
#' @section Notes:
#' @references
#'
#' @examples
#' \dontrun{
#'  TBD
#' }
#' @name add_latent_spatial
NULL

#' @name add_latent_spatial
#' @rdname add_latent_spatial
#' @exportMethod add_latent_spatial
#' @export
methods::setGeneric(
  "add_latent_spatial",
  signature = methods::signature("x"),
  function(x,...) standardGeneric("add_latent_spatial"))

#' @name add_latent_spatial
#' @rdname add_latent_spatial
#' @usage \S4method{add_latent_spatial}{BiodiversityDistribution}(x)
methods::setMethod(
  "add_latent_spatial",
  methods::signature(x = "BiodiversityDistribution"),
  function(x, ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"))
    # Finally add to data to the BiodiversityDistribution object
    x$set_latent(type = '<Spatial>')
    # TODO:
    # Create a prototype with prior ranges and sigmas for the matern matrix
    return(x)
  }
)

#TODO: Add other dummy methods for latent factors
