#' @include utils.R bdproto.R bdproto-prior.R bdproto-priorlist.R bdproto-distributionmodel.R
NULL

#' Add priors to an existing distribution object
#'
#' @description
#' This function simply allows to add priors to an existing [distribution] object.
#' The supplied priors must be a [`PriorList-class`] object created through
#' calling [priors].
#' @note
#' Alternatively priors to environmental predictors can also directly added as parameter
#' via [add_predictors]
#' @param x [distribution] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param priors A [`PriorList-class`] object containing multiple priors.
#' @param ... Other parameters passed down.
#' @family prior
#' @aliases add_priors
#' @examples
#' \dontrun{
#'  pp <-  GLMNETPrior("forest")
#'  x <- distribution(background) |>
#'   add_priors(pp)
#'
#' }
#' @name add_priors
NULL

#' @name add_priors
#' @rdname add_priors
#' @exportMethod add_priors
#' @export
methods::setGeneric(
  "add_priors",
  signature = methods::signature("x"),
  function(x, priors = NULL, ...) standardGeneric("add_priors"))

#' @name add_priors
#' @rdname add_priors
#' @usage \S4method{add_priors}{BiodiversityDistribution, ANY}(x, priors, ...)
methods::setMethod(
  "add_priors",
  methods::signature(x = "BiodiversityDistribution"),
  function(x, priors = NULL, ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.null(priors) || inherits(priors, "PriorList") || class(x)[1] %in% getOption("ibis.priors")
    )

    # Convert to prior list object
    if( class(x)[1] %in% getOption("ibis.priors") ){
      priors <- priors(priors)
    }
    if(!is.null(priors)){
      x <- x$set_priors( priors )
    }
    # Return x with newly added priors
    x
  }
)

#' @name set_priors
#' @inherit add_priors
#' @inheritParams add_priors
#' @keywords deprecated
methods::setGeneric(
  "set_priors",
  signature = methods::signature("x"),
  function(x, priors = NULL, ...) standardGeneric("set_priors"))

#' @name set_priors
#' @inherit add_priors
#' @inheritParams add_priors
#' @keywords deprecated
#' @usage \S4method{set_priors}{BiodiversityDistribution, ANY}(x, priors, ...)
methods::setMethod(
  "set_priors",
  methods::signature(x = "BiodiversityDistribution"),
  function(x, priors = NULL, ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.null(priors) || inherits(priors, "PriorList") || inherits(priors, 'INLAPrior') || inherits(priors, 'GDBPrior')
    )
    message('Deprecated. Use add_priors ')
    add_priors(x, priors, ...)
  }
)

#' Remove existing priors from an existing distribution object
#'
#' @description
#' This function allows to remove priors from an existing [distribution] object.
#' In order to remove a set prior, the name of the prior has to be specified.
#' @param x [distribution] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param names A [`vector`] or [`character`] object for priors to be removed.
#' @param ... Other parameters passed down
#' @aliases rm_priors
#' @family prior
#' @examples
#' \dontrun{
#'  # Add prior
#'  pp <-  GLMNETPrior("forest")
#'  x <- distribution(background) |>
#'   add_priors(pp)
#'  # Remove again
#'  x <- x |> rm_priors("forest")
#' }
#' @name rm_priors
NULL

#' @name rm_priors
#' @rdname rm_priors
#' @exportMethod rm_priors
#' @export
methods::setGeneric(
  "rm_priors",
  signature = methods::signature("x"),
  function(x, names = NULL, ...) standardGeneric("rm_priors"))

#' @name rm_priors
#' @rdname rm_priors
#' @usage \S4method{rm_priors}{BiodiversityDistribution, ANY}(x, names, ...)
methods::setMethod(
  "rm_priors",
  methods::signature(x = "BiodiversityDistribution"),
  function(x, names = NULL, ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.null(names) || is.vector(names) || is.character(names)
    )
    x <- x$rm_priors(names)
    # Return x with newly added priors
    return(x)
  }
)

#### Get priors from model - ####

#' Create priors from an existing distribution model
#'
#' @description
#' Often it can make sense to fit an additional model to get a grasp on the range of
#' values that "beta" parameters can take. This function takes an existing [`BiodiversityDistribution-class`] object
#' and creates [`PriorList-class`] object from them. The resulting object can be used to add for instance [priors]
#' to a new model.
#' @note
#' Not all engines support priors in similar ways. See the vignettes and help pages on that topic!
#' @param mod A fitted [`DistributionModel-class`] object. If instead a [`BiodiversityDistribution-class`] object
#' is passed to this function, it simply returns the contained priors used for estimation (if any).
#' @param target_engine A [`character`] for which the priors should be created.
#' @param ... Other parameters passed down.
#' @family prior
#' @aliases get_priors
#' @examples
#' \dontrun{
#'  mod <- distribution(background) |>
#'     add_predictors(covariates) |>
#'     add_biodiversity_poipo(points) |>
#'     engine_inlabru() |>
#'     train()
#'  get_priors(mod, target_engine = "BART")
#' }
#' @name get_priors
NULL

#' @name get_priors
#' @rdname get_priors
#' @exportMethod get_priors
#' @export
methods::setGeneric(
  "get_priors",
  signature = methods::signature("mod", "target_engine"),
  function(mod, target_engine, ...) standardGeneric("get_priors"))

#' @name get_priors
#' @rdname get_priors
#' @usage \S4method{get_priors}{ANY, character}(mod, target_engine, ...)
methods::setMethod(
  "get_priors",
  methods::signature(mod = "ANY", target_engine = "character"),
  function(mod, target_engine, ...) {
    assertthat::assert_that(inherits(mod, "DistributionModel") || inherits(mod, "BiodiversityDistribution"),
                            is.character(target_engine)
    )
    if(inherits(mod, "BiodiversityDistribution")){ return( mod$get_priors() ) }
    check_package("scales")
    # Catch character for engines
    target_engine <- match.arg(toupper(target_engine), choices = c("BART", "<BART>",
                                                                   "GDB", "<GDB>",
                                                                   "INLA", "<INLA>",
                                                                   "BART", "<BART>",
                                                                   "STAN", "<STAN>",
                                                                   "BREG", "<BREG>",
                                                                   "GLMNET", "<GLMNET>",
                                                                   "XGBOOST","<XGBOOST>"), several.ok = FALSE)

    # Get the coefficients
    cofs <- mod$get_coefficients()
    if(nrow(cofs)==0) return(NULL) # Something went wrong here

    # Check type of coefficients
    has_weights <- length( grep("weight", names(cofs),ignore.case = TRUE) ) > 0
    has_beta <- length( grep("beta", names(cofs),ignore.case = TRUE) ) > 0
    has_sigma <- length( grep("sigma", names(cofs),ignore.case = TRUE) ) > 0

    # Now depending on the target engine, and formulate priors for the target engine
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Create prior object from coefficients.')

    pl <- list()
    for(i in 1:nrow(cofs)){
      sub <- cofs[i,]

      # --- BREG ---
      if(target_engine %in% c("BREG", "<BREG>")){
        if(has_weights){
          # Only weights available
          pl[[i]] <- BREGPrior(variable = sub$Feature,
                               # Hyper
                               hyper = NULL,
                               ip = sub$Weight
          )
        } else {
          # Has beta coefficients
          # Also sigma? If so, alter inclusion probability for sellner-Prior  ?
          pl[[i]] <- BREGPrior(variable = sub$Feature,
                               # Hyper
                               hyper = sub$Beta,
                               ip = NULL
          )
        }
      }
      # --- BART ---
      if(target_engine %in% c("BART", "<BART>")){
        if(has_weights){
          pl[[i]] <- BARTPrior(variable = sub$Feature,
                               hyper = scales::rescale(cofs$Weights, to = c(1e-2, 1))[i] )
        } else {
          # Beta coefficients, convert them to absolute, log-transform and then scale
          val <- scales::rescale( log(abs(cofs$Beta)), to = c(1e-02, 1))[i]
          pl[[i]] <- BARTPrior(variable = sub$Feature,
                               hyper = val )
        }
      }
      # --- GDB ---
      if(target_engine %in% c("GDB", "<GDB>")){
        if(has_beta){
          pl[[i]] <- GDBPrior(variable = sub$Feature,
                              hyper = ifelse(sub$Beta >0, "increasing", "decreasing"))
        } else {
          # Unfortunately can't set any priors based on weights alone, thus set to none
          pl[[i]] <- GDBPrior(variable = sub$Feature,
                              hyper = "none")
        }
      }
      # --- GLMNET ---
      if(target_engine %in% c("GLMNET", "<GLMNET>")){
        if(has_beta){
          pl[[i]] <- GLMNETPrior(variable = sub$Feature,
                                 # If absolute beta coefficient is larger than 0.05, use a rescaled value
                                 # for defining the regularization constant
                                 hyper = ifelse(abs(sub$Beta)>0.05, abs(scales::rescale(abs(cofs$Beta), to = c(0, 1))-1)[i], 1),
                                 lims = c(-Inf, Inf))
        }
      }
      # --- XGBOOST ---
      if(target_engine %in% c("XGBOOST", "<XGBOOST>")){
        if(has_beta){
          val <- ifelse(sub$Beta >0,
                        ifelse(sub$Beta > 0.05, "increasing", "positive"),
                        ifelse(sub$Beta < -0.05, "decreasing", "negative")
                        )
          pl[[i]] <- XGBPrior(variable = sub$Feature, hyper = val)
        } else {
          # Unfortunately can't set any priors based on weights alone, thus set to none
          pl[[i]] <- XGBPrior(variable = sub$Feature,
                              hyper = "none")
        }
      }
      # --- STAN ---
      if(target_engine %in% c("STAN", "<STAN>")){
        if(has_beta && has_sigma){
          pl[[i]] <- STANPrior(variable = sub$Feature,
                               type = "gaussian",
                               hyper = c(sub$Beta, sub$Sigma)
          )
        } else if(has_beta && !has_sigma) {
          pl[[i]] <- STANPrior(variable = sub$Feature,
                               type = "gaussian",
                               hyper = c(sub$Beta, 0.05)
          )
        }
      }
      # --- INLA ---
      if(target_engine %in% c("INLA", "<INLA>")){
        if(has_beta && has_sigma){
          # Specify with precision priors
          pl[[i]] <- INLAPrior(variable = sub$Feature,
                               type = "gaussian",
                               hyper = c(sub$Beta, (sub$Sigma)^-2)
          )
        } else if(has_beta && !has_sigma) {
          # Now sigma. Set to a default
          pl[[i]] <- INLAPrior(variable = sub$Feature,
                               type = "gaussian",
                               hyper = c(sub$Beta, 0.001)
          )
        } # INLA not working with weights afaik
      }
    }

    # Final checks and conversion
    assertthat::assert_that(is.list(pl))
    if(length(pl)>0){
      op <- priors(pl)

      assertthat::assert_that(
        inherits(op, "PriorList"),
        op$length() == nrow(cofs)
      )
    } else { op <- NULL }

    return( op )
  }
)

