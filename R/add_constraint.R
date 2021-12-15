#' @include utils.R bdproto-biodiversityscenario.R
NULL

#' Add a constrain to an existing scenario
#'
#' @description This function adds a constrain to a [`BiodiversityScenario-class`] object to
#' constrain (future) projections. These constrains can for instance be constrains on a possible
#' dispersal distance, connectivity between identified patches or limitations on species vital rates.
#' **Most constrains require pre-calculated thresholds to present in the [`BiodiversityScenario-class`] object!**
#' @param mod A [`BiodiversityScenario`] object with specified predictors
#' @param method A [`character`] indicating the type of constrain to be added to the scenario. See details
#' @param value A [`numeric`] value specifying a fixed constrain or constant in units "m". Default: NULL
#' @param resistance A [`RasterLayer`] object describing a resistance surface or barrier for use in connectivity constrains. Default: NULL
#' @param ... passed on parameters
#'
#' @seealso [`add_constrain_dispersal`], [`add_constrain_connectivity`]
#' @details
#' Currently this method functions as a wrapper to support the definition of further modelling constraints.
#' Supported are the options for dispersal and connectivity constrains:
#' [-] sdd_fixed - Applies a fixed uniform dispersal distance per modelling timestep
#' [-] sdd_nexpkernel - Applies a dispersal distance using a negative exponential kernel from its origin,
#' [-] migclim - Applies the dispersal algorithm MigClim to the modelled objects
#' [-] hardbarrier - Defines a hard barrier to any dispersal events
#'
#' A comprehensive overview of the benefits of including dispersal constrains in species distribution models
#' can be found in Bateman et al. (2013).
#' @references Evans, M.E.K., Merow, C., Record, S., McMahon, S.M., Enquist, B.J., 2016. Towards Process-based Range Modeling of Many Species. Trends Ecol. Evol. 31, 860–871. https://doi.org/10.1016/j.tree.2016.08.005
#' @references Bateman, B. L., Murphy, H. T., Reside, A. E., Mokany, K., & VanDerWal, J. (2013). Appropriateness of full‐, partial‐and no‐dispersal scenarios in climate change impact modelling. Diversity and Distributions, 19(10), 1224-1234.
#'
#' @name add_constrain
#' @aliases add_constrain
#' @keywords scenario
#' @exportMethod add_constrain
#' @export
NULL
methods::setGeneric("add_constrain",
                    signature = methods::signature("mod"),
                    function(mod, method,...) standardGeneric("add_constrain"))

#' @name add_constrain
#' @rdname add_constrain
#' @usage \S4method{add_constrain}{BiodiversityScenario, character}(mod, method)
methods::setMethod(
  "add_constrain",
  methods::signature(mod = "BiodiversityScenario"),
  function(mod, method, ...){
    assertthat::assert_that(
      inherits(mod, "BiodiversityScenario"),
      !is.Waiver(mod$get_predictors()),
      is.character(method)
    )
    # Match method
    method <- match.arg(arg = method,
                        choices = c("sdd_fixed", "sdd_nexpkernel", "migclim","hardbarrier"), several.ok = FALSE)

    # Now call the respective functions individually
    o <- switch(method,
                  # Fixed dispersal
                  "sdd_fixed" = add_constrain_dispersal(mod, method = "sdd_fixed", ...),
                  # Short-distance dispersal
                  "sdd_nexpkernel" = add_constrain_dispersal(mod, method = "sdd_nexpkernel", ...),
                  # Using the migclim package
                  "migclim" = add_constrain_dispersal(mod, method = "migclim", ...),
                  # --- #
                 "hardbarrier" = add_constrain_connectivity(mod, method = "hardbarrier", ...)
                  )
    return(o)
  }
)

#' @title Adds a dispersal constrain to a scenario object
#' @name add_constrain_dispersal
#' @aliases add_constrain_dispersal
#' @inheritParams add_constrain
#' @keywords scenario
#' @exportMethod add_constrain_dispersal
#' @export
NULL
methods::setGeneric("add_constrain_dispersal",
                    signature = methods::signature("mod"),
                    function(mod, method, value = NULL, ...) standardGeneric("add_constrain_dispersal"))

#' @name add_constrain_dispersal
#' @rdname add_constrain_dispersal
#' @usage \S4method{add_constrain_dispersal}{BiodiversityScenario, character, numeric}(mod, method, value)
methods::setMethod(
  "add_constrain_dispersal",
  methods::signature(mod = "BiodiversityScenario"),
  function(mod, method, value = NULL, ...){
    assertthat::assert_that(
      inherits(mod, "BiodiversityScenario"),
      !is.Waiver(mod$get_predictors()),
      is.character(method),
      is.null(value) || is.numeric(value)
    )
    # Match method
    method <- match.arg(arg = method,
                        choices = c("sdd_fixed", "sdd_nexpkernel", "migclim"), several.ok = FALSE)

    # Check if there is already a dispersal constrain, if yes raise warning
    if(!is.Waiver(mod$get_constraints())){
      # If there are any dispersal constrains in there, raise warning
      if(any( "dispersal" %in% names(mod$get_constraints()) )){
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','yellow','Overwriting existing dispersal constraint.')
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
      # See if sth. can be taken from here
      # http://doi.wiley.com/10.1111/j.1365-2699.2012.02737.x
      # http://doi.wiley.com/10.1111/gcb.13251
      # http://www.nature.com/articles/nclimate3414
    }
    if(method == "migclim"){
      # Using the MigClim package for calculating any transitions and
      # This requires prior calculated Thresholds!
      out <- add_constrain_MigClim(mod = mod, ...)
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
    is.logical(resistance) || is.Raster(resistance),
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
    is.logical(resistance) || is.Raster(resistance),
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

# ------------------------ #
#### Connectivity constraints ####

#' @title Adds a connectivity constrain to a scenario object
#' @name add_constrain_connectivity
#' @aliases add_constrain_connectivity
#' @inheritParams add_constrain
#' @keywords scenario
#' @exportMethod add_constrain_connectivity
#' @export
NULL
methods::setGeneric("add_constrain_connectivity",
                    signature = methods::signature("mod"),
                    function(mod, method, value = NULL, resistance = NULL, ...) standardGeneric("add_constrain_connectivity"))

#' @name add_constrain_connectivity
#' @rdname add_constrain_connectivity
#' @usage \S4method{add_constrain_connectivity}{BiodiversityScenario, character, numeric, RasterLayer}(mod, method, value, resistance)
methods::setMethod(
  "add_constrain_connectivity",
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
