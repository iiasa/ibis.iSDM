#' @include utils.R bdproto.R waiver.R
NULL

#' Add latent spatial factor to the model
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param spatial_model A [`character`] describing what kind of spatial to be fitted (Option: 'spde' | 'iCAR')
#' @param priors A [`Prior-List`] object supplied to the latent effect. NULL equating default priors
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
  function(x, spatial_model = 'spde', priors = NULL, ...) standardGeneric("add_latent_spatial"))

#' @name add_latent_spatial
#' @rdname add_latent_spatial
#' @usage \S4method{add_latent_spatial}{BiodiversityDistribution}(x)
methods::setMethod(
  "add_latent_spatial",
  methods::signature(x = "BiodiversityDistribution"),
  function(x, spatial_model = 'spde', priors = NULL, ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            spatial_model %in% c('spde','iCAR'),
                            is.null(priors) || inherits(priors, 'PriorList')
                            )
    if(spatial_model=='iCAR') stop('Needs to be debugged. ID of mesh not linked.')

        # If priors have been set, save them in distribution object
    if(!is.null(priors)) {
      assertthat::assert_that(priors$varnames() == 'spde' && (!is.null(priors$exists('spde','prior.range')) || !is.null(priors$exists('spde','prior.sigma')) ),
                              msg = 'Priors for spatial latent effect misspeficied (required spde | prior.range/prior.sigma)'  )
      # Add priors
      x <- x$set_priors(priors)
    }
    # Add to data to the BiodiversityDistribution object
    x$set_latent(type = '<Spatial>', spatial_model)
  }
)

#TODO: Add other dummy methods for latent factors we might want to specify
