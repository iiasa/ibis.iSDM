#' @include utils.R bdproto.R waiver.R
NULL

#' Add latent spatial effect to the model equation
#'
#' @description Adding latent spatial error terms, spatial trends or other autocovariates
#' to a model can account for spatial dependence in response or covariates
#' @details There are several different options which depend on the engine used. In case a
#' unsupported method for an engine is chosen this is modified to the next similar method
#' Available are:
#' "spde" - stochastic partial differential equation (SPDE) for [`INLA-engine`]. SPDE effects aim at
#' capturing the variation of the response variable in space, once all of the covariates are accounted for.
#' Examining the spatial distribution of the spatial error can reveal which covariates might be missing.
#' For example, if elevation is positively correlated with the response variable, but is not included in the model,
#' we could see a higher posterior mean in areas with higher elevation.
#'
#' "car" - conditional autocorrelative errors (CAR) for [`INLA-engine`]
#' "poly"- spatial trend correction by adding spatial coordinates as polynominals in the model. Available for all Engines.
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param method A [`character`] describing what kind of spatial effect is to be added to the model.
#' @param priors A [`Prior-List`] object supplied to the latent effect. NULL equating default priors
#' @param ... Other parameters passed down
#'
#' @references Fletcher, R., & Fortin, M. (2018). Spatial ecology and conservation modeling. Springer International Publishing.
#' @references Mendes, P., Velazco, S. J. E., de Andrade, A. F. A., & Júnior, P. D. M. (2020). Dealing with overprediction in species distribution models: How adding distance constraints can improve model accuracy. Ecological Modelling, 431, 109180.
#'
#' @examples
#' \dontrun{
#'  distribution(background) %>% add_latent_spatial()
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
  function(x, method = 'spde', priors = NULL, ...) standardGeneric("add_latent_spatial"))

#' @name add_latent_spatial
#' @rdname add_latent_spatial
#' @usage \S4method{add_latent_spatial}{BiodiversityDistribution}(x)
methods::setMethod(
  "add_latent_spatial",
  methods::signature(x = "BiodiversityDistribution"),
  function(x, method = 'spde', priors = NULL) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.character(method) && !missing(method),
                            is.null(priors) || inherits(priors, 'PriorList')
                            )

    # Match the spatial method
    method <- match.arg(method, c('spde', 'car', 'poly'), several.ok = FALSE)

    # Messager
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding latent spatial terms...')

    # If priors have been set, save them in distribution object
    if(!is.null(priors) & method == 'spde') {
      assertthat::assert_that(priors$varnames() == 'spde' && (!is.null(priors$exists('spde','prior.range')) || !is.null(priors$exists('spde','prior.sigma')) ),
                              msg = 'Priors for spatial latent effect misspeficied (required for spde | prior.range/prior.sigma)'  )
      # Add priors
      x <- x$set_priors(priors)
    }
    # Add to data to the BiodiversityDistribution object
    x$set_latent(type = '<Spatial>', method)
  }
)
