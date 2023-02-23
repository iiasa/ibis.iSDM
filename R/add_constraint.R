#' @include utils.R bdproto-biodiversityscenario.R
NULL

#' Add a constraint to an existing \code{scenario}
#'
#' @description This function adds a constrain to a [`BiodiversityScenario-class`] object to
#' constrain (future) projections. These constrains can for instance be constraints on a possible
#' dispersal distance, connectivity between identified patches or limitations on species adaptability.
#' **Most constrains require pre-calculated thresholds to present in the [`BiodiversityScenario-class`] object!**
#' @param mod A [`BiodiversityScenario`] object with specified predictors.
#' @param method A [`character`] indicating the type of constraints to be added to the scenario. See details.
#' @param value For many dispersal [`constrain`] this is set as [`numeric`] value specifying a
#' fixed constrain or constant in units \code{"m"} (Default: \code{NULL}). For kissmig the value needs to
#' give the number of iteration steps (or within year migration steps).
#' For adaptability constraints this parameter specifies the extent (in units of standard deviation) to which extrapolations
#' should be performed.
#' @param type A [`character`] indicating the type used in the method. See for instance [kissmig::kissmig].
#' @param layer A [`Raster`] object that can be used for boundary constraints (Default: \code{NULL}).
#' @param pext [`numeric`] indicator for [`kissmig`] of the probability a colonized cell becomes uncolonized,
#' i.e., the species gets locally extinct (Default: \code{0.1}).
#' @param pcor [`numeric`] probability that corner cells are considered in the 3x3 neighbourhood (Default: \code{0.2}).
#' @param ... passed on parameters. See also the specific methods for adding constraints.
#'
#' @seealso [`add_constraint_dispersal`], [`add_constraint_connectivity`], [`add_constraint_adaptability`],[`add_constraint_boundary`]
#' @details
#' Currently this method functions as a wrapper to support the definition of further modelling constraints.
#' Supported are the options for dispersal and connectivity constraints:
#' * \code{sdd_fixed} - Applies a fixed uniform dispersal distance per modelling timestep.
#' * \code{sdd_nexpkernel} - Applies a dispersal distance using a negative exponential kernel from its origin.
#' * \code{kissmig} - Applies the kissmig stochastic dispersal model. Requires [kissmig] package. Applied at each modelling time step.
#' * \code{migclim} - Applies the dispersal algorithm MigClim to the modelled objects. Requires [MigClim] package.
#' * \code{hardbarrier} - Defines a hard barrier to any dispersal events.
#' * \code{nichelimit} - Specifies a limit on the environmental niche to only allow a modest amount of extrapolation beyond the known occurrences. This
#' can be particular useful to limit the influence of increasing marginal responses and avoid biologically unrealistic projections.
#' * \code{boundary} - Applies a hard boundary constraint on the projection.
#'
#' A comprehensive overview of the benefits of including dispersal constrains in species distribution models
#' can be found in Bateman et al. (2013).
#' @references
#' * Bateman, B. L., Murphy, H. T., Reside, A. E., Mokany, K., & VanDerWal, J. (2013). Appropriateness of full‐, partial‐and no‐dispersal scenarios in climate change impact modelling. Diversity and Distributions, 19(10), 1224-1234.
#' * Nobis MP and Normand S (2014) KISSMig - a simple model for R to account for limited migration in analyses of species distributions. Ecography 37: 1282-1287.
#' * Mendes, P., Velazco, S. J. E., de Andrade, A. F. A., & Júnior, P. D. M. (2020). Dealing with overprediction in species distribution models: How adding distance constraints can improve model accuracy. Ecological Modelling, 431, 109180.
#' @name add_constraint
#' @family constraint
#' @aliases add_constraint
#' @keywords scenario
#' @exportMethod add_constraint
#' @export
NULL
methods::setGeneric("add_constraint",
                    signature = methods::signature("mod"),
                    function(mod, method, ...) standardGeneric("add_constraint"))

#' @name add_constraint
#' @rdname add_constraint
#' @usage \S4method{add_constraint}{BiodiversityScenario, character}(mod, method)
methods::setMethod(
  "add_constraint",
  methods::signature(mod = "BiodiversityScenario"),
  function(mod, method, ...) {
    assertthat::assert_that(
      inherits(mod, "BiodiversityScenario"),
      !is.Waiver(mod$get_predictors()),
      is.character(method)
    )
    # Match method
    method <- match.arg(arg = method,
                        choices = c("sdd_fixed", "sdd_nexpkernel", "kissmig", "migclim",
                                    "hardbarrier","boundary",
                                    "nichelimit"), several.ok = FALSE)

    # Now call the respective functions individually
    o <- switch(method,
                  # Fixed dispersal
                  "sdd_fixed" = add_constraint_dispersal(mod, method = "sdd_fixed", ...),
                  # Short-distance dispersal
                  "sdd_nexpkernel" = add_constraint_dispersal(mod, method = "sdd_nexpkernel", ...),
                  # Add kissmig dispersal
                  "kissmig" = add_constraint_dispersal(mod, method = "kissmig", ...),
                  # Using the migclim package
                  "migclim" = add_constraint_dispersal(mod, method = "migclim", ...),
                  # --- #
                 "hardbarrier" = add_constraint_connectivity(mod, method = "hardbarrier", ...),
                  # --- #
                  "nichelimit" = add_constraint_adaptability(mod, method = "nichelimit", ...),
                  # --- #
                  "boundary" = add_constraint_boundary(mod, ...)
                  )
    return(o)
  }
)

# ------------------------ #
#### Dispersal constraints ####

#' @title Adds a dispersal constrain to a scenario object
#' @name add_constraint_dispersal
#' @aliases add_constraint_dispersal
#' @inheritParams add_constraint
#' @family constraint
#' @keywords scenario
#' @exportMethod add_constraint_dispersal
#' @export
NULL
methods::setGeneric("add_constraint_dispersal",
                    signature = methods::signature("mod"),
                    function(mod, method, value = NULL, type = NULL, ...) standardGeneric("add_constraint_dispersal"))

#' @name add_constraint_dispersal
#' @rdname add_constraint_dispersal
#' @usage \S4method{add_constraint_dispersal}{BiodiversityScenario, character, numeric}(mod, method, value)
methods::setMethod(
  "add_constraint_dispersal",
  methods::signature(mod = "BiodiversityScenario"),
  function(mod, method, value = NULL, type = NULL, ...){
    assertthat::assert_that(
      inherits(mod, "BiodiversityScenario"),
      !is.Waiver(mod$get_predictors()),
      is.character(method),
      is.null(value) || is.numeric(value),
      is.null(type) || is.character(type)
    )
    # Match method
    method <- match.arg(arg = method,
                        choices = c("sdd_fixed", "sdd_nexpkernel", "kissmig", "migclim"), several.ok = FALSE)

    # Other arguments supplied
    dots <- list(...)
    argnames <- names(dots)

    # Check if there is already a dispersal constrain, if yes raise warning
    if(!is.Waiver(mod$get_constraints())){
      # If there are any dispersal constrains in there, raise warning
      if(any("dispersal" %in% names(mod$get_constraints()))){
        if(getOption('ibis.setupmessages')) myLog('[Estimation]', 'yellow', 'Overwriting existing dispersal constraint.')
      }
    }

    # Add processing method #
    # --- #
    cr <- list()
    if(method == "sdd_fixed"){
      # Short-distance dispersal (Fixed)
      assertthat::assert_that(
        is.numeric(value), msg = "Fixed short distance dispersal needs an annual mean disperal distance value."
      )
      cr[['dispersal']] <- list(method = method,
                                params = c("mean_dispersal_distance" = value))
    } else if(method == "sdd_nexpkernel") {
      # Negative exponential kernel
      assertthat::assert_that(
        is.numeric(value), msg = "Short distance negative exponential kernal dispersal needs an annual mean disperal distance value."
      )
      cr[['dispersal']] <- list(method = method,
                                params = c("mean_dispersal_distance" = value))
    } else if(method == "kissmig"){
      # Check parameters to be correct
      check_package("kissmig")
      # Gather some default parameters
      if(is.null(type)) type <- "DIS" else match.arg(type, c("DIS", "FOC", "LOC", "NOC"), several.ok = FALSE)
      assertthat::assert_that(
        is.numeric(value),
        value > 0, msg = "For kissmig the value needs to give the number of iteration steps (or within time migration steps)."
      )
      # probability [0,1] a colonized cell becomes uncolonized between iteration steps, i.e., the species gets locally extinct
      if("pext" %in% argnames) pext <- dots[["pext"]] else pext <- 0.1
      # probability [0,1] corner cells are considered in the 3x3 cell neighborhood. Following Nobis & Nomand 2014, 0.2 is recommended for circular spread
      if("pcor" %in% argnames) pcor <- dots[["pcor"]] else pcor <- 0.2

      if(getOption('ibis.setupmessages')) myLog('[Estimation]', 'green', 'KISSMIG options: iterations=',value,'| pext=', pext,'| pcor=', pcor)

      cr[['dispersal']] <- list(method = method,
                                params = c("iteration" = value,
                                           "type" = type,
                                           "signed" = FALSE,
                                           "pext" = pext,
                                           "pcor" = pcor
                                           ))

    }
    if(method == "migclim"){
      # Using the MigClim package for calculating any transitions and
      # This requires prior calculated Thresholds!
      out <- add_constraint_MigClim(mod = mod, ...)
      return(out)
    } else {
      # --- #
      new <- mod$set_constraints(cr)
      return(
        bdproto(NULL, new)
      )
    }

  }
)

#' Short-distance fixed dispersal function
#' @param baseline_threshold The [`RasterLayer`] with presence/absence information from a previous year
#' @param new_suit A new [`RasterLayer`] object
#' @param value A [`numeric`] value of the fixed dispersal threshold. In unit 'meters'.
#' @param resistance A resistance [`RasterLayer`] object with values to be omitted during distance calculation (Default: NULL)
#' @noRd
#' @keywords internal
.sdd_fixed <- function(baseline_threshold, new_suit, value, resistance = NULL){
  assertthat::assert_that(
    is.Raster(baseline_threshold), is.Raster(new_suit),
    raster::compareRaster(baseline_threshold, new_suit),
    is.numeric(value),
    is.null(resistance) || is.Raster(resistance),
    # Check that baseline threshold raster is binomial
    length(unique(baseline_threshold))==2
  )

  # Set resistance layer to 0 if set to zero.
  if(is.Raster(resistance)){
    baseline_threshold[resistance == 1] <- 2
    # Set resistance to the value omitted
    resistance <- 2
  }
  # Grow baseline raster by the amount of value at max
  # Furthermore divide by value to get a normalized distance
  dis <- raster::gridDistance(baseline_threshold, origin = 1, omit = resistance)
  ras_dis <- raster::clamp(dis, lower = 0, upper = value) / value
  # Invert
  ras_dis <- abs(ras_dis - 1)

  # Now multiply the net suitability projection with this mask
  # Thus removing any grid cells outside
  out <- new_suit * ras_dis
  return(out)
}

#' Short-distance negative exponential kernel dispersal function
#' @param baseline_threshold The [`RasterLayer`] with presence/absence information from a previous year
#' @param new_suit A new [`RasterLayer`] object
#' @param value A [`numeric`] value of the fixed dispersal threshold. In unit 'meters'.
#' @param normalize Should a normalising constant be used for the exponential dispersal parameter. (Default: FALSE)
#' @param resistance A resistance [`RasterLayer`] object with values to be omitted during distance calculation (Default: NULL)
#' @noRd
#' @keywords internal
.sdd_nexpkernel <- function(baseline_threshold, new_suit, value, normalize = FALSE, resistance = NULL){
  assertthat::assert_that(
    is.Raster(baseline_threshold), is.Raster(new_suit),
    raster::compareRaster(baseline_threshold, new_suit),
    is.numeric(value),
    is.null(resistance) || is.Raster(resistance),
    # Check that baseline threshold raster is binomial
    length(unique(baseline_threshold))==2
  )

  # Set resistance layer to 0 if set to zero.
  if(is.Raster(resistance)){
    baseline_threshold[resistance == 1] <- 2
    # Set resistance to the value omitted
    resistance <- 2
  }

  # Inverse of mean dispersal distance
  alpha <- 1/value

  # Grow baseline raster by using an exponentially weighted kernel
  ras_dis <- raster::gridDistance(baseline_threshold, origin = 1, omit = resistance)
  if(normalize){
    # Normalized (with a constant) negative exponential kernel
    ras_dis <- raster::calc(ras_dis, fun = function(x) (1 / (2 * pi * value ^ 2)) * exp(-x / value) )
  } else {
    ras_dis <- raster::calc(ras_dis, fun = function(x) exp(-alpha * x))
  }

  # Now multiply the net suitability projection with this mask
  # Thus removing any non-suitable grid cells (0) and changing the value of those within reach
  out <- new_suit * ras_dis
  return(out)
}

#' Keep it simple migration calculation.
#' @param baseline_threshold The [`RasterLayer`] with presence/absence information from a previous year.
#' @param new_suit A new [`RasterLayer`] object.
#' @param params A [vector] or [list] with passed on parameter values.
#' @param resistance A resistance [`RasterLayer`] object with values to be omitted during distance calculation (Default: NULL).
#' @noRd
#' @keywords internal
.kissmig_dispersal <- function(baseline_threshold, new_suit, params, resistance = NULL){
  assertthat::assert_that(
    is.Raster(baseline_threshold), is.Raster(new_suit),
    raster::compareRaster(baseline_threshold, new_suit),
    is.vector(params) || is.list(params),
    is.null(resistance) || is.logical(resistance) || is.Raster(resistance),
    # Check that baseline threshold raster is binomial
    length(unique(baseline_threshold))==2
  )

  check_package('kissmig')
  if(!isNamespaceLoaded("kissmig")) { attachNamespace("kissmig");requireNamespace("kissmig") }

  # Set suitability layer to 0 if resistance layer is set
  if(is.Raster(resistance)){
    new_suit[resistance>0] <- 0
  }

  # Simulate kissmig for a given threshold and suitability raster
  km <- kissmig::kissmig(O = baseline_threshold,
                         # Rescale newsuit to 0-1
                         S = predictor_transform(new_suit, 'norm'),
                         it = as.numeric( params['iteration'] ),
                         type = params['type'],
                         pext = as.numeric(params['pext']),
                         pcor = as.numeric(params['pcor'])
                        )
  if(is.factor(km)) km <- raster::deratify(km, complete = TRUE)

  # Now multiply the net suitability projection with this mask
  # Thus removing any non-suitable grid cells (0) and changing the value of those within reach
  ns <- new_suit * km

  return(
    raster::stack(km, ns)
  )
}

# ------------------------ #
#### Connectivity constraints ####

#' @title Adds a connectivity constraint to a scenario object
#' @name add_constraint_connectivity
#' @aliases add_constraint_connectivity
#' @inheritParams add_constraint
#' @param resistance A [`RasterLayer`] object describing a resistance surface or barrier for use in connectivity constrains (Default: \code{NULL}).
#' @family constraint
#' @keywords scenario
#' @exportMethod add_constraint_connectivity
#' @export
NULL
methods::setGeneric("add_constraint_connectivity",
                    signature = methods::signature("mod"),
                    function(mod, method, value = NULL, resistance = NULL, ...) standardGeneric("add_constraint_connectivity"))

#' @name add_constraint_connectivity
#' @rdname add_constraint_connectivity
#' @usage \S4method{add_constraint_connectivity}{BiodiversityScenario, character, numeric, RasterLayer}(mod, method, value, resistance)
methods::setMethod(
  "add_constraint_connectivity",
  methods::signature(mod = "BiodiversityScenario"),
  function(mod, method, value = NULL, resistance = NULL, ...){
    assertthat::assert_that(
      inherits(mod, "BiodiversityScenario"),
      !is.Waiver(mod$get_predictors()),
      is.character(method),
      is.null(value) || is.numeric(value),
      is.Raster(resistance) || is.null(resistance)
    )
    # Match method
    method <- match.arg(arg = method,
                        choices = c("hardbarrier"), several.ok = FALSE)

    # Check if there is already a dispersal constrain, if yes raise warning
    if(!is.Waiver(mod$get_constraints())){
      # If there are any dispersal constrains in there, raise warning
      if(any( "connectivity" %in% names(mod$get_constraints()) )){
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','yellow','Overwriting existing connectivity constraint')
      }
    }

    # Add processing method #
    # --- #
    co <- list()
    if(method == "hardbarrier"){
      # Short-distance dispersal (Fixed)
      assertthat::assert_that(
        is.Raster(resistance),
        !is.null(resistance), msg = "Set a hard barrier via the resistance parameter."
      )
      # Check that resistance layer is a binary mask
      assertthat::assert_that(length(unique(resistance))<=2,
                              cellStats(resistance,'max')>0,
                              msg = "Resistance layer should be a binary mark.")
      co[['connectivity']] <- list(method = method,
                                params = c("resistance" = resistance))
    }
    # --- #
    new <- mod$set_constraints(co)
    return(
      bdproto(NULL, new)
    )
  }
)

# ------------------------ #
#### Adaptability constraints ####

#' @title Adds an adaptability constraint to a scenario object
#' @description
#' Adaptability constraints assume that suitable habitat for species in (future) projections might be unsuitable if
#' it is outside the range of conditions currently observed for the species.
#'
#' Currently only `nichelimit` is implemented, which adds a simple constrain on the predictor parameter space, which
#' can be defined through the \code{"value"} parameter. For example by setting it to \code{1} (Default), any projections
#' are constrained to be within the range of at maximum 1 standard deviation from the range of covariates used for model
#' training.
#' @name add_constraint_adaptability
#' @aliases add_constraint_adaptability
#' @inheritParams add_constraint
#' @param names A [`character`] vector with names of the predictors for which an adaptability threshold should be set (Default: \code{NULL} for all).
#' @param value A [`numeric`] value in units of standard deviation (Default: \code{1}).
#' @param increment A [`numeric`] constant that is added to value at every time step (Default: \code{0}).
#' Allows incremental widening of the niche space, thus opening constraints.
#' @family constraint
#' @keywords scenario
#' @exportMethod add_constraint_adaptability
#' @export
NULL
methods::setGeneric("add_constraint_adaptability",
                    signature = methods::signature("mod"),
                    function(mod, method = "nichelimit", names = NULL, value = 1, increment = 0, ...) standardGeneric("add_constraint_adaptability"))

#' @name add_constraint_adaptability
#' @rdname add_constraint_adaptability
#' @usage \S4method{add_constraint_adaptability}{BiodiversityScenario, character, character, numeric, numeric}(mod, method, names, value, increment)
methods::setMethod(
  "add_constraint_adaptability",
  methods::signature(mod = "BiodiversityScenario"),
  function(mod, method = "nichelimit", names = NULL, value = 1, increment = 0, ...){
    assertthat::assert_that(
      inherits(mod, "BiodiversityScenario"),
      !is.Waiver(mod$get_predictors()),
      is.character(method),
      is.null(names) || is.character(names),
      is.null(value) || is.numeric(value),
      is.numeric(increment)
    )
    # Match method
    method <- match.arg(arg = method,
                        choices = c("nichelimit"), several.ok = FALSE)

    # Add processing method #
    # --- #
    co <- list()
    if(method == "nichelimit"){
      # Add a constrain on parameter space, e.g. max 1 SD from training data covariates
      assertthat::assert_that(
        is.numeric(value),
        is.null(names) || is.character(names),
        value > 0, msg = "Specify a value threshold (SD) and names of predictors, for which
        we do not expect the species to persist."
      )
      if(is.null(names)) names <- character()
      co[['adaptability']] <- list(method = method,
                                   params = c("names" = names, "value" = value,
                                              "increment" = increment))
    }
    # --- #
    new <- mod$set_constraints(co)
    return(
      bdproto(NULL, new)
    )
  }
)

#' Adaptability constrain by applying a limit on extrapolation beyond the niche
#'
#' @param newdata A [`data.frame`] with the information about new data layers.
#' @param model A [`list`] created by the modelling object containing the full predictors and biodiversity predictors.
#' @param names A [`character`] or \code{NULL} of the names of predictors.
#' @param value A [`numeric`] value in units of standard deviation (Default: \code{1}).
#' @param increment A [`numeric`] constant that is added to value at every time step (Default: \code{0}).
#' Allows incremental widening of the niche space, thus opening constraints.
#' @param increment_step A [`numeric`] indicating the number of time increment should be applied.
#' @keywords internal
#' @noRd
.nichelimit <- function(newdata, model, names = NULL, value = 1, increment = 0, increment_step = 1){
  assertthat::assert_that(
    is.data.frame(newdata),
    is.list(model),
    is.numeric(as.numeric(value)),
    is.null(names) || is.na(names) || is.character(names),
    is.numeric(as.numeric(increment)),
    is.numeric(as.numeric(increment_step))
  )
  # Check that names are present if set
  if(is.null(names) || is.na(names)) names <- model$predictors_names
  if(is.character(names) ) assertthat::assert_that(all(names %in% model$predictors_names))
  # Convert numeric paramters to numeric to be sure
  value <- as.numeric(value)
  increment <- as.numeric(increment)
  increment_step <- as.numeric(increment_step)
  # --- #
  # Now calculate the range across each target predictor and occurrence dataset
  df <- data.frame()
  for(id in names(model$biodiversity)){
    sub <- model$biodiversity[[id]]
    # Which are presence data
    is_presence <- which(sub$observations[['observed']] > 0)
    df <- rbind(df,
                sub$predictors[is_presence, names])
  }
  rr <- sapply(df, function(x) range(x, na.rm = TRUE))   # Calculate ranges
  rsd <- sapply(df, function(x) sd(x, na.rm = TRUE))   # Calculate standard deviation

  # Apply value and increment if set
  rsd <- rsd * (value + (increment*increment_step))
  rr[1,] <- rr[1,] - rsd; rr[2,] <- rr[2,] + rsd

  # Now 'clamp' all predictor values beyond these names to 0, e.g. partial out
  nd <- newdata
  for(n in names){
    # Calc min
    min_ex <- which(nd[,n] < rr[1,n])
    max_ex <- which(nd[,n] > rr[2,n])
    if(length(min_ex)>0) nd[min_ex,n] <- NA
    if(length(max_ex)>0) nd[max_ex,n] <- NA
    # FIXME Or rather do a smooth logistic decay for less extreme?
  }
  return(nd)
}

# ------------------------ #
#### Boundary constraints ####

#' @title Adds a boundary constraint to a scenario object
#' @description
#' The purpose of boundary constraints is to limit a future projection within a specified area
#' (such as for example a range or ecoregion). This can help to limit unreasonable projections into geographic space.
#'
#' Similar to boundary constraints it is also possible to define a \code{"zone"} for the scenario projections, similar
#' as was done for model training. The difference to a boundary constraint is that the boundary constraint is applied posthoc
#' as a hard cut on any projection, while the zones would allow any projection (and other constraints) to be applied within
#' the zone.
#' **Note: Setting a boundary constraint for future projections effectively potentially suitable areas!**
#' @name add_constraint_boundary
#' @aliases add_constraint_boundary
#' @inheritParams add_constraint
#' @param layer A [`RasterLayer`] object with the same extent as the model background. Has to be binary and
#' is used for a posthoc masking of projected grid cells.
#' @family constraint
#' @keywords scenario
#' @exportMethod add_constraint_boundary
#' @export
NULL
methods::setGeneric("add_constraint_boundary",
                    signature = methods::signature("mod", "layer"),
                    function(mod, layer, ...) standardGeneric("add_constraint_boundary"))

#' @name add_constraint_boundary
#' @rdname add_constraint_boundary
#' @usage \S4method{add_constraint_boundary}{BiodiversityScenario, sf, character}(mod, layer, method)
methods::setMethod(
  "add_constraint_boundary",
  methods::signature(mod = "BiodiversityScenario", layer = "sf"),
  function(mod, layer, method = "boundary", ...){
    assertthat::assert_that(
      inherits(mod, "BiodiversityScenario"),
      inherits(layer, "sf"),
      is.character(method)
    )

    # Rasterize the layer
    # First try and dig out a layer from a predictor dataset if found
    if(inherits( mod$get_predictors(), "PredictorDataSet")){
      ras <- mod$get_predictors()$get_data() |> stars_to_raster()
      ras <- ras[[1]]
    } else {
      # Try and get the underlying model and its predictors
      ras <- mod$get_model()$get_data()
    }
    assertthat::assert_that(is.Raster(ras))
    bb <- try({ raster::rasterize(layer, ras, 1)},silent = TRUE)
    if(inherits(bb, "try-error")) stop("Provide a rasterized layer of the boundary constraint!")

    # Call again
    o <- add_constraint_boundary(mod, layer = bb, method = method, ..)

    return( o )
  }
)

#' @name add_constraint_boundary
#' @rdname add_constraint_boundary
#' @usage \S4method{add_constraint_boundary}{BiodiversityScenario, ANY, character}(mod, layer, method)
methods::setMethod(
  "add_constraint_boundary",
  methods::signature(mod = "BiodiversityScenario", layer = "ANY"),
  function(mod, layer, method = "boundary", ...){
    assertthat::assert_that(
      inherits(mod, "BiodiversityScenario"),
      is.Raster(layer),
      is.character(method)
    )

    # Check that layer is a single RasterLayer
    if(!inherits(layer, "RasterLayer")){
      assertthat::assert_that(raster::nlayers(layer) == 1)
      layer <- layer[[1]]
    }

    # Add processing method #
    # --- #
    co <- list()
    if(method == "boundary"){
      # Add a constrain on parameter space, e.g. max 1 SD from training data covariates
      assertthat::assert_that(
        length( unique( layer )) <=2
      )
      # If length of values is greater than 1, remove everything else by setting it to NA
      if( length( unique( layer )) >1 ){
        layer[layer<1] <- NA
      }
      co[['boundary']] <- list(method = method,
                                   params = c("layer" = layer))
    }
    # --- #
    new <- mod$set_constraints(co)
    return(
      bdproto(NULL, new)
    )
  }
)
