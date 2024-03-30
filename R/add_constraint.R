#' @include class-biodiversityscenario.R
NULL

#' Add a constraint to an existing \code{scenario}
#'
#' @description This function adds a constrain to a
#' [`BiodiversityScenario-class`] object to constrain (future) projections.
#' These constrains can for instance be constraints on a possible dispersal
#' distance, connectivity between identified patches or limitations on species
#' adaptability.
#'
#' **Most constrains require pre-calculated thresholds to present in the [`BiodiversityScenario-class`] object!**
#'
#' @param mod A [`BiodiversityScenario`] object with specified predictors.
#' @param method A [`character`] indicating the type of constraints to be added
#' to the scenario. See details for more information.
#' @param ... passed on parameters. See also the specific methods for adding
#' constraints.
#'
#' @details Constraints can be added to scenario objects to increase or decrease
#' the suitability of a given area for the target feature. This function acts
#' as a wrapper to add these constraints. Currently supported are the
#' following options:
#'
#' **Dispersal**:
#' * \code{sdd_fixed} - Applies a fixed uniform dispersal distance per modelling timestep.
#' * \code{sdd_nexpkernel} - Applies a dispersal distance using a negative exponential kernel from its origin.
#' * \code{kissmig} - Applies the kissmig stochastic dispersal model. Requires \code{`kissmig`} package. Applied at each modelling time step.
#' * \code{migclim} - Applies the dispersal algorithm MigClim to the modelled objects. Requires \code{"MigClim"} package.
#'
#' A comprehensive overview of the benefits of including dispersal constrains in
#' species distribution models can be found in Bateman et al. (2013).
#'
#' **Connectivity**:
#' * \code{hardbarrier} - Defines a hard barrier to any dispersal events. By
#' definition this sets all values larger
#' than \code{0} in the barrier layer to \code{0} in the projection. Barrier has
#' to be provided through the \code{"resistance"} parameter.
#' * \code{resistance} - Allows the provision of a static or dynamic layer that is
#' multiplied with the projection at each time step. Can for example be used to
#' reduce the suitability of any given area (using pressures not included in the model).
#' The respective layer(s) have to be provided through the \code{"resistance"} parameter.
#' Provided layers are incorporated as \code{abs(resistance - 1)} and multiplied with
#' the prediction.
#'
#' **Adaptability**:
#' * \code{nichelimit} - Specifies a limit on the environmental niche to only allow
#' a modest amount of extrapolation beyond the known occurrences. This can be particular
#' useful to limit the influence of increasing marginal responses and avoid biologically
#' unrealistic projections.
#'
#' **Boundary and size**:
#' * \code{boundary} - Applies a hard boundary constraint on the projection, thus
#' disallowing an expansion of a range outside the provide layer. Similar as specifying
#' projection limits (see [`distribution`]), but can be used to specifically constrain a
#' projection within a certain area (e.g. a species range or an island).
#' * \code{minsize} - Allows to specify a certain size that must be satisfied in
#' order for a thresholded patch to be occupied. Can be thought of as a minimum
#' size requirement. See `add_constraint_minsize()` for the required parameters.
#' * \code{threshold} - Applies the set threshold as a constrain directly on the
#' suitability projections. Requires a threshold to be set.
#'
#' @returns Adds constraints data to a [`BiodiversityScenario`] object.
#'
#' @references
#' * Bateman, B. L., Murphy, H. T., Reside, A. E., Mokany, K., & VanDerWal, J. (2013).
#' Appropriateness of full‐, partial‐and no‐dispersal scenarios in climate change impact
#'  modelling. Diversity and Distributions, 19(10), 1224-1234.
#' * Nobis MP and Normand S (2014) KISSMig - a simple model for R to account for
#' limited migration in analyses of species distributions. Ecography 37: 1282-1287.
#' * Mendes, P., Velazco, S. J. E., de Andrade, A. F. A., & Júnior, P. D. M. (2020).
#' Dealing with overprediction in species distribution models: How adding distance
#' constraints can improve model accuracy. Ecological Modelling, 431, 109180.
#'
#' @family constraint
#' @keywords scenario
#'
#' @examples
#' \dontrun{
#' # Assumes that a trained 'model' object exists
#'  mod <- scenario(model) |>
#'   add_predictors(env = predictors, transform = 'scale', derivates = "none") |>
#'   add_constraint_dispersal(method = "kissmig", value = 2, pext = 0.1) |>
#'   project()
#' }
#'
#' @name add_constraint
NULL

#' @rdname add_constraint
#' @export
methods::setGeneric("add_constraint",
                    signature = methods::signature("mod"),
                    function(mod, method, ...) standardGeneric("add_constraint"))

#' @rdname add_constraint
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
                                    "hardbarrier","resistance",
                                    "boundary", "minsize", "threshold",
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
                  "resistance" = add_constraint_connectivity(mod, method = "resistance", ...),
                  # --- #
                  "nichelimit" = add_constraint_adaptability(mod, method = "nichelimit", ...),
                  # --- #
                  "boundary" = add_constraint_boundary(mod, ...),
                  # --- #
                  "threshold" = add_constraint_threshold(mod, ...),
                  # --- #
                  "minsize" = add_constraint_minsize(mod, ...)
    )
    return(o)
  }
)

# ------------------------ #
#### Dispersal constraints ####

#' Add dispersal constraint to an existing \code{scenario}
#'
#' @inheritParams add_constraint
#' @param value For many dispersal \code{"constrain"} this is set as [`numeric`]
#' value specifying a fixed constrain or constant in units \code{"m"}
#' (Default: \code{NULL}). For kissmig the value needs to give the number of
#' iteration steps (or within year migration steps). For adaptability
#' constraints this parameter specifies the extent (in units of standard
#' deviation) to which extrapolations should be performed.
#' @param type A [`character`] indicating the type used in the method. See for
#' instance \code{`kissmig`}.
#'
#' @details
#' **Dispersal**:
#' Parameters for \code{'method'}:
#' * \code{sdd_fixed} - Applies a fixed uniform dispersal distance per modelling timestep.
#' * \code{sdd_nexpkernel} - Applies a dispersal distance using a negative exponential kernel from its origin.
#' #' The negative exponential kernel is defined as:
#' \deqn{f(x) = \frac{1}{2 \pi a^2} e^{-\frac{x}{a}}}{fx = 1 / (2 * pi * a^2) * exp(-x / a)}
#' where \eqn{a} is the mean dispersal distance (in m) divided by 2.
#' * \code{kissmig} - Applies the kissmig stochastic dispersal model. Requires \code{`kissmig`} package. Applied at each modelling time step.
#' * \code{migclim} - Applies the dispersal algorithm MigClim to the modelled objects. Requires \code{"MigClim"} package.
#'
#' A comprehensive overview of the benefits of including dispersal constrains in
#' species distribution models can be found in Bateman et al. (2013).
#'
#' The following additional parameters can bet set:
#' * \code{pext}: [`numeric`] indicator for \code{`kissmig`} of the probability a
#' colonized cell becomes uncolonised, i.e., the species gets locally extinct
#' (Default: \code{0.1}).
#' * \code{pcor}: [`numeric`] probability that corner cells are considered in the
#' 3x3 neighbourhood (Default: \code{0.2}).
#'
#' @note
#' Unless otherwise stated, the default unit of supplied distance values (e.g. average dispersal
#' distance) should be in \code{"m"}.
#' @references
#' * Bateman, B. L., Murphy, H. T., Reside, A. E., Mokany, K., & VanDerWal, J. (2013).
#' Appropriateness of full‐, partial‐and no‐dispersal scenarios in climate change impact
#'  modelling. Diversity and Distributions, 19(10), 1224-1234.
#' @family constraint
#' @keywords scenario
#'
#'@name add_constraint_dispersal
NULL

#' @rdname add_constraint_dispersal
#' @export
methods::setGeneric("add_constraint_dispersal",
                    signature = methods::signature("mod"),
                    function(mod, method, value = NULL, type = NULL, ...) standardGeneric("add_constraint_dispersal"))

#' @rdname add_constraint_dispersal
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
        if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]', 'yellow', 'Overwriting existing dispersal constraint.')
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

      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]', 'green', 'KISSMIG options: iterations=',value,'| pext=', pext,'| pcor=', pcor)

      cr[['dispersal']] <- list(method = method,
                                params = c("iteration" = value,
                                           "type" = type,
                                           "signed" = FALSE,
                                           "pext" = pext,
                                           "pcor" = pcor
                                           ))

    }
    out <- mod$clone(deep = TRUE)
    if(method == "migclim"){
      # Using the MigClim package for calculating any transitions and
      # This requires prior calculated Thresholds!
      out <- add_constraint_MigClim(mod = out, ...)
    } else {
      # --- #
      out <- out$set_constraints(cr)
    }
    return(out)
  }
)

#' Short-distance fixed dispersal function
#'
#' @param baseline_threshold The [`SpatRaster`] with presence/absence
#' information from a previous year.
#' @param new_suit A new [`SpatRaster`] object.
#' @param value A [`numeric`] value of the fixed dispersal threshold. In unit
#' \code{'meters'}.
#' @param resistance A resistance [`SpatRaster`] object with values to be
#' omitted during distance calculation (Default: \code{NULL}).
#'
#' @noRd
#'
#' @keywords internal
.sdd_fixed <- function(baseline_threshold, new_suit, value, resistance = NULL){
  assertthat::assert_that(
    is.Raster(baseline_threshold), is.Raster(new_suit),
    is_comparable_raster(baseline_threshold, new_suit),
    is.numeric(value),
    is.null(resistance) || is.Raster(resistance)
  )

  # Check for small lon-lat values
  if(terra::is.lonlat(baseline_threshold)){
    if(value < 1){
      message('Very small average dispersal vlaue provided. Check that they are in unit m!')
    }
  }

  # Get original baseline threshold
  ori.tr <- baseline_threshold
  ori.tr[ori.tr>0] <- 1

  # Set resistance layer to 0 if set to zero.
  if(is.Raster(resistance)){
    baseline_threshold[resistance == 1] <- 2
    # Set resistance to the value omitted
    resistance <- 2
    baseline_threshold <- terra::mask(baseline_threshold, resistance)
  }
  # Grow baseline raster by the amount of value at max
  # Furthermore divide by value to get a normalized distance
  dis <- terra::gridDist(baseline_threshold, target = 1)
  ras_dis <- terra::clamp(dis, lower = 0, upper = value) / value
  # Invert
  ras_dis <- abs(ras_dis - 1)

  # Now multiply the net suitability projection with this mask
  # Thus removing any grid cells outside
  out <- new_suit * ras_dis
  # Mask with original so as to retain non-zero values
  out <- terra::mask(out, ori.tr)
  return(out)
}

#' Short-distance negative exponential kernel dispersal function
#'
#' @param baseline_threshold The [`SpatRaster`] with presence/absence information
#' from a previous year.
#' @param new_suit A new [`SpatRaster`] object.
#' @param value A [`numeric`] value of the fixed dispersal threshold. In unit \code{'meters'}.
#' @param normalize Should a normalising constant be used for the exponential
#' dispersal parameter (Default: \code{FALSE}).
#' @param resistance A resistance [`SpatRaster`] object with values to be omitted
#' during distance calculation (Default: \code{NULL}).
#'
#' @noRd
#'
#' @keywords internal
.sdd_nexpkernel <- function(baseline_threshold, new_suit, value, normalize = TRUE, resistance = NULL){
  assertthat::assert_that(
    is.Raster(baseline_threshold), is.Raster(new_suit),
    is_comparable_raster(baseline_threshold, new_suit),
    is.numeric(value),
    is.logical(normalize),
    is.null(resistance) || is.Raster(resistance),
    # Check that baseline threshold raster is binomial
    length(unique(baseline_threshold)[,1])==2
  )

  # Check for small lon-lat values
  if(terra::is.lonlat(baseline_threshold)){
    if(value < 1){
      message('Very small average dispersal vlaue provided. Check that they are in unit m!')
    }
  }

  # Get original baseline threshold
  ori.tr <- baseline_threshold
  ori.tr[ori.tr>0] <- 1

  # Set resistance layer to 0 if set to zero.
  if(is.Raster(resistance)){
    baseline_threshold[resistance == 1] <- 2
    # Set resistance to the value omitted
    resistance <- 2
    baseline_threshold <- terra::mask(baseline_threshold, resistance)
  }

  # Divide alpha values by 2
  alpha <- value/2

  # Grow baseline raster by using an exponentially weighted kernel
  ras_dis <- terra::gridDist(baseline_threshold, target = 1)
  # Normalized (with a constant) negative exponential kernel
  ras_dis <- terra::app(ras_dis, fun = function(x) (1 / (2 * pi * value ^ 2)) * exp(-x / value) )
  # Equivalent to alpha = 1/value and
  # ras_dis <- terra::app(ras_dis, fun = function(x) exp(-alpha * x))
  if(normalize){
    ras_dis <- predictor_transform(ras_dis, option = 'norm')
  }

  # Now multiply the net suitability projection with this mask Thus removing any
  # non-suitable grid cells (0) and changing the value of those within reach
  out <- new_suit * ras_dis
  return(out)
}

#' Keep it simple migration calculation.
#'
#' @param baseline_threshold The [`SpatRaster`] with presence/absence
#' information from a previous year.
#' @param new_suit A new [`SpatRaster`] object.
#' @param params A [vector] or [list] with passed on parameter values.
#' @param resistance A resistance [`SpatRaster`] object with values to be
#' omitted during distance calculation (Default: \code{NULL}).
#'
#' @noRd
#'
#' @keywords internal
.kissmig_dispersal <- function(baseline_threshold, new_suit, params, resistance = NULL){
  assertthat::assert_that(
    is.Raster(baseline_threshold), is.Raster(new_suit),
    is_comparable_raster(baseline_threshold, new_suit),
    is.vector(params) || is.list(params),
    is.null(resistance) || is.logical(resistance) || is.Raster(resistance),
    # Check that baseline threshold raster is binomial
    length(unique(baseline_threshold)[,1])==2
  )

  check_package('kissmig')
  if(!isNamespaceLoaded("kissmig")) { attachNamespace("kissmig");requireNamespace("kissmig") }

  # Set suitability layer to 0 if resistance layer is set
  if(is.Raster(resistance)){
    new_suit[resistance>0] <- 0
  }

  # Simulate kissmig for a given threshold and suitability raster
  km <- kissmig::kissmig(O = terra_to_raster( baseline_threshold ),
                         # Rescale newsuit to 0-1
                         S = predictor_transform(new_suit, 'norm') |>
                           terra_to_raster(),
                         it = as.numeric( params['iteration'] ),
                         type = params['type'],
                         pext = as.numeric(params['pext']),
                         pcor = as.numeric(params['pcor'])
                        )
  # Convert to terra again
  km <- terra::rast(km)
  if(is.factor(km)) km <- terra::as.int(km)

  # Now multiply the net suitability projection with this mask Thus removing any
  # non-suitable grid cells (0) and changing the value of those within reach
  ns <- new_suit * km

  return(
    c(km, ns)
  )
}

# ------------------------ #
#### Connectivity constraints ####

#' Adds a connectivity constraint to a scenario object.
#'
#' @inheritParams add_constraint
#' @param value For many dispersal \code{"constrain"} this is set as [`numeric`]
#' value specifying a fixed constrain or constant in units \code{"m"}
#' (Default: \code{NULL}). For kissmig the value needs to give the number of
#' iteration steps (or within year migration steps). For adaptability
#' constraints this parameter specifies the extent (in units of standard
#' deviation) to which extrapolations should be performed.
#' @param resistance A [`SpatRaster`] object describing a resistance surface or
#' barrier for use in connectivity constrains (Default: \code{NULL}).
#'
#' @details
#' * \code{hardbarrier} - Defines a hard barrier to any dispersal events. By
#' definition this sets all values larger
#' than \code{0} in the barrier layer to \code{0} in the projection. Barrier has
#' to be provided through the \code{"resistance"} parameter.
#' * \code{resistance} - Allows the provision of a static or dynamic layer that is
#' multiplied with the projection at each time step. Can for example be used to
#' reduce the suitability of any given area (using pressures not included in the model).
#' The respective layer(s) have to be provided through the \code{"resistance"} parameter.
#' Provided layers are incorporated as \code{abs(resistance - 1)} and multiplied with
#' the prediction.
#'
#'
#' @family constraint
#' @keywords scenario
#'
#'@name add_constraint_connectivity
NULL

#' @rdname add_constraint_connectivity
#' @export
methods::setGeneric("add_constraint_connectivity",
                    signature = methods::signature("mod"),
                    function(mod, method, value = NULL, resistance = NULL, ...) standardGeneric("add_constraint_connectivity"))

#' @rdname add_constraint_connectivity
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
                        choices = c("hardbarrier", "resistance"), several.ok = FALSE)

    # Check if there is already a dispersal constrain, if yes raise warning
    if(!is.Waiver(mod$get_constraints())){
      # If there are any dispersal constrains in there, raise warning
      if(any( "connectivity" %in% names(mod$get_constraints()) )){
        if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Setup]','yellow','Overwriting existing connectivity constraint')
      }
    }

    # Add processing method #
    # --- #
    co <- list()
    if(method == "hardbarrier"){
      # Assert hard barrier
      assertthat::assert_that(
        is.Raster(resistance),
        !is.null(resistance), msg = "Set a hard barrier via the resistance parameter."
      )
      # Check that resistance layer is a binary mask
      assertthat::assert_that(length(unique(resistance))<=2,
                              terra::global(resistance,'max', na.rm = TRUE)>0,
                              msg = "Resistance layer should be a binary mark with values 0/1.")
      co[['connectivity']] <- list(method = method,
                                params = c("resistance" = resistance))
    } else if(method == "resistance"){
      # Flexible resistance layer
      assertthat::assert_that(
        is.Raster(resistance),
        !is.null(resistance), msg = "The method resistance requires a specified resistance raster."
      )
      # If raster is stack with multiple layers, ensure that time
      if(terra::nlyr(resistance)>1){
        # Check that layers have a z dimension and fall within the timeperiod
        startend <- mod$get_timeperiod()
        assertthat::assert_that( !is.null( terra::time(resistance) ),
                                 all( range(terra::time(resistance))==startend ),
                                 msg = "If a stack of layers is supplied as resistance, it needs a Z value of equal length to the predictors!")
      }
      times <- terra::time(resistance)
      # If resistance layer is bigger than 1, normalize
      if(any(terra::global(resistance, 'max', na.rm = TRUE)>1)){
        if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Setup]','yellow','Resistance values larger than 1. Normalizing...')
        resistance <- predictor_transform(resistance, option = "norm")
      }
      resistance <- abs( resistance - 1 ) # Invert
      if(!is.null(times)) terra::time(resistance) <- times # Reset times again if found

      co[['connectivity']] <- list(method = method,
                                   params = c("resistance" = resistance))
    }
    # --- #
    new <- mod$clone(deep = TRUE)
    new$set_constraints(co)
    return(new)
  }
)

# ------------------------ #
#### Adaptability constraints ####

#' Adds an adaptability constraint to a scenario object
#'
#' @description Adaptability constraints assume that suitable habitat for
#' species in (future) projections might be unsuitable if it is outside the
#' range of conditions currently observed for the species.
#'
#' Currently only `nichelimit` is implemented, which adds a simple constrain on
#' the predictor parameter space, which can be defined through the
#' \code{"value"} parameter. For example by setting it to \code{1} (Default),
#' any projections are constrained to be within the range of at maximum 1
#' standard deviation from the range of covariates used for model training.
#'
#' @inheritParams add_constraint
#' @param names A [`character`] vector with names of the predictors for which an
#'   adaptability threshold should be set (Default: \code{NULL} for all).
#' @param value A [`numeric`] value in units of standard deviation (Default:
#'   \code{1}).
#' @param increment A [`numeric`] constant that is added to value at every time
#'   step (Default: \code{0}). Allows incremental widening of the niche space,
#'   thus opening constraints.
#'
#' @family constraint
#' @keywords scenario
#'
#' @examples
#' \dontrun{
#' scenario(fit) |>
#'  add_constraint_adaptability(value = 1)
#' }
#' @name add_constraint_adaptability
NULL

#' @rdname add_constraint_adaptability
#' @export
methods::setGeneric("add_constraint_adaptability",
                    signature = methods::signature("mod"),
                    function(mod, method = "nichelimit", names = NULL, value = 1, increment = 0, ...) standardGeneric("add_constraint_adaptability"))

#' @rdname add_constraint_adaptability
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
      if(is.null(names)) names <- NA
      co[['adaptability']] <- list(method = method,
                                   params = c("names" = names, "value" = value,
                                              "increment" = increment))
    }
    # --- #
    new <- mod$clone(deep = TRUE)
    new <- new$set_constraints(co)
    return(new)
  }
)

#' Adaptability constrain by applying a limit on extrapolation beyond the niche
#'
#' @param newdata A [`data.frame`] with the information about new data layers.
#' @param model A [`list`] created by the modelling object containing the full
#' predictors and biodiversity predictors.
#' @param names A [`character`] or \code{NULL} of the names of predictors.
#' @param value A [`numeric`] value in units of standard deviation (Default: \code{1}).
#' @param increment A [`numeric`] constant that is added to value at every time
#' step (Default: \code{0}). Allows incremental widening of the niche space,
#' thus opening constraints.
#' @param increment_step A [`numeric`] indicating the number of time increment
#' should be applied.
#'
#' @noRd
#'
#' @keywords internal
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
  # Convert numeric parameters to numeric to be sure
  value <- as.numeric(value)
  increment <- as.numeric(increment)
  increment_step <- as.numeric(increment_step)
  # --- #
  # Now calculate the range across each target predictor and occurrence dataset
  df <- data.frame()
  for(id in names(model$biodiversity)){
    sub <- model$biodiversity[[id]]
    # Which are presence data
    is_presence <- which(factor_to_numeric(sub$observations[['observed']]) > 0)
    df <- rbind(df,
                sub$predictors[is_presence, sub$predictors_names])
  }
  rr <- sapply(df, function(x) range(x, na.rm = TRUE))   # Calculate ranges
  rsd <- sapply(df, function(x) stats::sd(x, na.rm = TRUE))   # Calculate standard deviation

  # Apply value and increment if set
  rsd <- rsd * (value + (increment*increment_step))
  rr[1,] <- rr[1,] - rsd; rr[2,] <- rr[2,] + rsd

  # Now 'clamp' all predictor values beyond these names to 0, e.g. partial out
  nd <- newdata
  for(n in names){
    if(!(n %in% names(rr))) next() # If variable not present in model frame, skip
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
#### Size constraints ####

#' Adds a size constraint on a scenario
#'
#' @description
#' This function applies a minimum size constraint on a `scenario()` created
#' object. The rationale here is that for a given species isolated habitat patches
#' smaller than a given size might not be viable / unrealistic for a species
#' to establish a (long-term) presence.
#'
#' The idea thus is to apply a constraint in that only patches bigger than a
#' certain size are retained between timesteps.
#' It has thus the potential to reduce subsequent colonizations of neighbouring
#' patches.
#'
#' @inheritParams add_constraint
#' @param value A [`numeric`] value describing the minimum amount of area of a
#' given patch
#' @param unit A [`character`] of the unit of area. Options available are
#' \code{km2} (Default) and \code{ha}.
#' @param establishment_step A [`logical`] flag indicating whether a given patch
#' is only to be removed if wasn't small in a previous time step (not yet
#' implemented!)
#'
#' @details
#' Area values in a specific unit need to be supplied.
#'
#' @note
#' *This function requires that a scenario has a set `threshold()`!*
#'
#' @family constraint
#' @keywords scenario
#'
#' @examples
#' \dontrun{
#' scenario(fit) |>
#'  add_predictors(future_covariates) |>
#'  threshold() |>
#'  add_constraint_minsize(value = 1000, unit = "km2") |>
#'  project()
#' }
#'
#' @name add_constraint_minsize
NULL

#' @rdname add_constraint_minsize
#' @export
methods::setGeneric("add_constraint_minsize",
                    signature = methods::signature("mod", "value"),
                    function(mod, value, unit = "km2", establishment_step = FALSE, ...) standardGeneric("add_constraint_minsize"))

#' @rdname add_constraint_minsize
methods::setMethod(
  "add_constraint_minsize",
  methods::signature(mod = "BiodiversityScenario", value = "numeric"),
  function(mod, value, unit = "km2", establishment_step = FALSE, ...){
    assertthat::assert_that(
      inherits(mod, "BiodiversityScenario"),
      !is.Waiver(mod$get_predictors()),
      is.null(value) || is.numeric(value),
      is.character(unit),
      is.logical(establishment_step)
    )
    # Match unit
    unit <- match.arg(arg = unit,
                        choices = c("km2", "ha", "pixel"), several.ok = FALSE)

    if(unit=="pixel"){
      assertthat::assert_that(value>1,
                              msg = "For unit pixel supply values > 1.")
    }
    # Add processing method #
    # --- #
    co <- list()
    co[['min_size']] <- list(method = "min_size",
                             params = c("value" = value,
                                        "unit" = unit,
                                        "establishment_step" = establishment_step))
    # --- #
    new <- mod$clone(deep = TRUE)
    new <- new$set_constraints(co)
    return(new)
  }
)


# ------------------------ #
#### Boundary constraints ####

#' Adds a boundary constraint to a scenario object
#'
#' @description The purpose of boundary constraints is to limit a future
#' projection within a specified area (such as for example a range or
#' ecoregion). This can help to limit unreasonable projections into geographic
#' space.
#'
#' Similar to boundary constraints it is also possible to define a \code{"zone"}
#' for the scenario projections, similar as was done for model training. The
#' difference to a boundary constraint is that the boundary constraint is
#' applied posthoc as a hard cut on any projection, while the zones would allow
#' any projection (and other constraints) to be applied within the zone.
#' **Note: Setting a boundary constraint for future projections effectively potentially suitable areas!**
#'
#' @inheritParams add_constraint
#' @param layer A [`SpatRaster`] or [`sf`] object with the same extent as the
#'   model background. Has to be binary and is used for a posthoc masking of
#'   projected grid cells.
#'
#' @family constraint
#' @keywords scenario
#'
#' @examples
#' \dontrun{
#' # Add scenario constraint
#' scenario(fit) |> add_constraint_boundary(range)
#' }
#'
#' @name add_constraint_boundary
NULL

#' @rdname add_constraint_boundary
#' @export
methods::setGeneric("add_constraint_boundary",
                    signature = methods::signature("mod", "layer"),
                    function(mod, layer, ...) standardGeneric("add_constraint_boundary"))

#' @rdname add_constraint_boundary
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
    bb <- try({ terra::rasterize(layer, ras, 1)}, silent = TRUE)
    if(inherits(bb, "try-error")) stop("Provide a rasterized layer of the boundary constraint!")

    # Call again
    new <- mod$clone(deep = TRUE)
    new <- add_constraint_boundary(new, layer = bb, method = method, ...)
    return( new )
  }
)

#' @rdname add_constraint_boundary
methods::setMethod(
  "add_constraint_boundary",
  methods::signature(mod = "BiodiversityScenario", layer = "ANY"),
  function(mod, layer, method = "boundary", ...){
    assertthat::assert_that(
      inherits(mod, "BiodiversityScenario"),
      is.Raster(layer),
      is.character(method)
    )

    # Check that layer is a single SpatRaster
    if(!inherits(layer, "SpatRaster")){
      assertthat::assert_that(terra::nlyr(layer) == 1)
      layer <- layer[[1]]
    }

    # Add processing method #
    # --- #
    co <- list()
    if(method == "boundary"){
      # Add a constrain on parameter space, e.g. max 1 SD from training data
      # covariates
      assertthat::assert_that(
        length( unique( layer )) <=2
      )
      # If length of values is greater than 1, remove everything else by setting
      # it to NA
      if( length( unique( layer )) >1 ){
        layer[layer<1] <- NA
    }
      co[['boundary']] <- list(method = method,
                                   params = c("layer" = layer))
    }
    # --- #
    new <- mod$clone(deep = TRUE)
    new <- new$set_constraints(co)
    return( new )
  }
)

#' Adds a threshold constraint to a scenario object
#'
#' @description
#' This option adds a [`threshold()`] constraint to a scenario projection,
#' thus effectively applying the threshold as mask to each projection step made
#' during the scenario projection.
#'
#' Applying this constraint thus means that the \code{"suitability"} projection is
#' clipped to the threshold. This method requires
#' the `threshold()` set for a scenario object.
#'
#' It could be in theory possible to re calculate the threshold for each time step
#' based on supplied parameters or even observation records. So far this option has
#' not been necessary to implement.
#'
#' @note
#' Threshold values are taken from the original fitted model.
#'
#' @inheritParams add_constraint
#' @param updatevalue A [`numeric`] indicating to what the masked out values (those outside)
#' the threshold should become (Default: \code{NA}).
#'
#' @family constraint
#' @keywords scenario
#'
#' @examples
#' \dontrun{
#' # Add scenario constraint
#' scenario(fit) |> threshold() |>
#' add_constraint_threshold()
#' }
#'
#' @name add_constraint_threshold
NULL

#' @rdname add_constraint_threshold
#' @export
methods::setGeneric("add_constraint_threshold",
                    signature = methods::signature("mod"),
                    function(mod, updatevalue = NA, ...) standardGeneric("add_constraint_threshold"))

#' @rdname add_constraint_threshold
methods::setMethod(
  "add_constraint_threshold",
  methods::signature(mod = "BiodiversityScenario"),
  function(mod, updatevalue = NA, ...){
    assertthat::assert_that(
      inherits(mod, "BiodiversityScenario"),
      is.na(updatevalue) || is.numeric(updatevalue)
    )

    # Add processing method #
    # --- #
    co <- list()
    co[['threshold']] <- list(method = "threshold",
                             params = c("updatevalue" = updatevalue))
    # --- #
    new <- mod$clone(deep = TRUE)
    new <- new$set_constraints(co)
    return( new )
  }
)
