#' @include utils.R bdproto.R waiver.R
NULL

#' Add latent spatial factor to the model
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param spatial_model a [`character`] describing what kind of spatial to be fitted (Option: 'spde' | 'iCAR')
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
  function(x, spatial_model = 'spde', ...) standardGeneric("add_latent_spatial"))

#' @name add_latent_spatial
#' @rdname add_latent_spatial
#' @usage \S4method{add_latent_spatial}{BiodiversityDistribution}(x)
methods::setMethod(
  "add_latent_spatial",
  methods::signature(x = "BiodiversityDistribution"),
  function(x, spatial_model = 'spde', ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            spatial_model %in% c('spde','iCAR') )
    if(spatial_model=='iCAR') stop('Needs to be debugged. ID of mesh not linked.')
    # Finally add to data to the BiodiversityDistribution object
    x$set_latent(type = '<Spatial>', spatial_model)
    # TODO:
    # Create a prototype with prior ranges and sigmas for the matern matrix
    return(x)
  }
)

#TODO: Add other dummy methods for latent factors
