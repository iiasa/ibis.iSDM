#' @include utils.R bdproto.R waiver.R
NULL

#' Add latent spatial effect to the model equation
#'
#' @description
#' In general we understand under latent spatial effects the occurrence of spatial dependency in the observations, which
#' might either be caused by spatial biases, similarities in the underlying sampling processes or unmeasured latent
#' covariates, e.g. those that have not been quantified.
#'
#' This package supports a range of different spatial effects, however they differ from another by their impact on
#' the estimated prediction. Some effects simply add the spatial dependence as covariate, others make use
#' of spatial random effects to account for spatial dependence in the predictions. By default these effects are
#' added to each dataset as covariate or shared spatial field (e.g. SPDE). See details for an
#' explanation of the available options.
#'
#' @details
#' There are several different options some of which depend on the engine used. In case a
#' unsupported method for an engine is chosen this is modified to the next similar method
#'
#' Available are:
#'
#' [*] \code{"spde"} - stochastic partial differential equation (SPDE) for [`INLA-engine`] and [`INLABRU-engine`].
#' SPDE effects aim at capturing the variation of the response variable in space, once all of the covariates are accounted for.
#' Examining the spatial distribution of the spatial error can reveal which covariates might be missing. For example,
#' if elevation is positively correlated with the response variable, but is not included in the model,
#' we could see a higher posterior mean in areas with higher elevation. Note that calculations of
#' SPDE's can be computationally costly.
#' * \code{"car"} - conditional autocorrelative errors (CAR) for [`INLA-engine`]. Not yet implemented in full.
#' * \code{"kde"} - additional covariate of the kernel density of input point observations.
#' * \code{"poly"} - spatial trend correction by adding coordinates as polynominal transformation. Available for all Engines.
#' * \code{"nnd"} - nearest neighbour distance. This function calculates the euclidean distance from each grid cell
#' to the nearest other grid cell with values. Applied across all datasets in the [`BiodiversityDistribution-class`]) object.
#' Available for all Engines.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param method A [`character`] describing what kind of spatial effect is to be added to the model. See details.
#' @param priors A [`Prior-List`] object supplied to the latent effect. Relevant only for [`engine_inla`] and \code{NULL} equates the use of default priors.
#' @param separate_spde A [`logical`] parameter indicating whether, in the case of SPDE effects, separate effects
#' for each likelihood are being fitted. Default (\code{FALSE}) uses a copy of the first added likelihood.
#' @param ... Other parameters passed down
#'
#' @references
#' * Fletcher, R., & Fortin, M. (2018). Spatial ecology and conservation modeling. Springer International Publishing.
#' * Mendes, P., Velazco, S. J. E., de Andrade, A. F. A., & Júnior, P. D. M. (2020). Dealing with overprediction in species distribution models: How adding distance constraints can improve model accuracy. Ecological Modelling, 431, 109180.
#'
#' @keywords latent
#' @examples
#' \dontrun{
#'  distribution(background) %>% add_latent_spatial(method = "poly")
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
  function(x, method = 'spde', priors = NULL, separate_spde = FALSE, ...) standardGeneric("add_latent_spatial"))

#' @name add_latent_spatial
#' @rdname add_latent_spatial
#' @usage \S4method{add_latent_spatial}{BiodiversityDistribution}(x)
methods::setMethod(
  "add_latent_spatial",
  methods::signature(x = "BiodiversityDistribution"),
  function(x, method = 'spde', priors = NULL, separate_spde = FALSE, ...) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.character(method) && !missing(method),
                            is.null(priors) || inherits(priors, 'PriorList'),
                            is.logical(separate_spde)
                            )

    # Match the spatial method
    method <- match.arg(method, c('spde', 'car', 'poly', 'nnd', 'kde'), several.ok = FALSE)

    # Messager
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding latent spatial terms...')

    # If priors have been set, save them in distribution object
    if(!is.null(priors) & method == 'spde') {
      assertthat::assert_that(any(priors$varnames() == 'spde') && (!is.null(priors$exists('spde','prior.range')) || !is.null(priors$exists('spde','prior.sigma')) ),
                              msg = 'Priors for spatial latent effect misspeficied (required for spde | prior.range/prior.sigma)'  )
      # Add priors
      x <- x$set_priors(priors)
    }
    # Add to data to the BiodiversityDistribution object
    x$set_latent(type = '<Spatial>', method, separate_spde)
  }
)

#' Function to remove a latent effect
#'
#' @description
#' This is just a wrapper function for removing specified offsets from a [`BiodiversityDistribution-class`]) object.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @seealso add_latent_spatial
#' @keywords latent, internal
#' @name rm_latent
NULL

#' @name rm_latent
#' @rdname rm_latent
#' @exportMethod rm_latent
#' @export
methods::setGeneric(
  "rm_latent",
  signature = methods::signature("x"),
  function(x) standardGeneric("rm_latent"))

#' @name rm_latent
#' @rdname rm_latent
#' @usage \S4method{rm_latent}{BiodiversityDistribution}(x)
methods::setMethod(
  "rm_latent",
  methods::signature(x = "BiodiversityDistribution"),
  function(x) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution") )
    # If no offset can be found, just return proto object
    if(is.Waiver(x$latentfactors)){ return(x) }

    # Messenger
    if(getOption('ibis.setupmessages')) myLog('[Setup]','yellow','Removing latent effects.')

    # Now remove the offset
    x$rm_latent()
  }
)
