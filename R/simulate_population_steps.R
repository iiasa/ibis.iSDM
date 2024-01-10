#' @include bdproto-biodiversityscenario.R
NULL

#' Simulate population dynamics following the steps approach
#'
#' @description This function adds a flag to a [`BiodiversityScenario-class`] object
#' to indicate that species abundances are to be simulated based on the expected
#' habitat suitability, as well as demography, density-dependence and dispersal
#' information. The simulation is done using the *steps* package (Visintin et
#' al. 2020) and conducted after a habitat suitability projection has been
#' created. *steps* is a spatially explicit population models coded mostly in R.
#'
#' For a detailed description of *steps* parameters, please see the respective reference
#' and help files. Default assumptions underlying this wrapper are presented in the details
#'
#' @param mod A [`BiodiversityScenario`] object with specified predictors.
#' @param replicates A [`numeric`] vector of the number of replicates (Default: \code{1}).
#' @param vital_rates A symmetrical demographic matrix. Should have column and row
#' names equivalent to the vital stages that are to be estimated.
#' @param carrying_capacity Either [`SpatRaster`] or a [`numeric`] estimate of the
#' maximum carrying capacity, e.g. how many adult individual are likely to occur
#' per grid cell. If set to [`numeric`], then carrying capacity is estimated up
#' to a maximum set (*Note: a more clever way would be to use a species-area relationship
#' for scaling. This is not yet implemented*).
#' @param initial A [`SpatRaster`] giving the initial population size. If not
#' provided, then initial populations are guessed (see details) from the projected
#' suitability rasters (Default: \code{NULL}).
#' @param dispersal A dispersal object defined by the \code{steps} package (Default: \code{NULL}).
#' @param density_dependence Specification of density dependence defined by the
#' \code{steps} package (Default: \code{NULL}).
#' @param include_suitability A [`logical`] flag on whether the projected suitability
#' estimates should be used (Default: \code{TRUE}) or only the initial conditions
#' set to the first time step.
#'
#' @details
#' In order for this function to work the *steps* package has to be installed
#' separately. Instructions to do so can be found on [github](https://github.com/steps-dev/steps).
#'
#' If initial population lifestages are not provided, then they are estimated
#' assuming a linear scaling with suitability, a \code{50:50} split between
#' sexes and a \code{1:3} ratio of adults to juveniles. The provision of
#' different parameters is highly encouraged!
#'
#' @note
#' The *steps* package has multiple options for simulating species population and
#' not all possible options are represented in this wrapper.
#'
#' Furthermore, the package still makes use of the \code{raster} package for much
#' of its internal data processing. Since *ibis.iSDM* switched to [terra] a while
#' ago, there can be efficiency problems as layers need to be translated between
#' packages.
#'
#' @returns Adds flag to a [`BiodiversityScenario`] object to indicate that
#' further simulations are added during projection.
#'
#' @references
#' * Visintin, C., Briscoe, N. J., Woolley, S. N., Lentini, P. E., Tingley, R.,
#'  Wintle, B. A., & Golding, N. (2020). steps: Software for spatially and
#'  temporally explicit population simulations. Methods in Ecology and Evolution,
#'  11(4), 596-603. https://doi.org/10.1111/2041-210X.13354
#'
#' @family constraint
#' @keywords scenario
#'
#' @examples
#' \dontrun{
#' # Define vital rates
#' vt <- matrix(c(0.0,0.5,0.75,
#'                0.5,0.2,0.0,
#'                0.0,0.5,0.9),
#'                nrow = 3, ncol = 3, byrow = TRUE)
#' colnames(vt) <- rownames(vt) <- c('juvenile','subadult','adult')
#'
#' # Assumes that a trained 'model' object exists
#'  mod <- scenario(model) |>
#'   add_predictors(env = predictors, transform = 'scale',
#'                  derivates = "none") |>
#'   # Use Vital rates here, but note the other parameters!
#'   simulate_population_steps(vital_rates = vt) |>
#'   project()
#' }
#'
#' @name simulate_population_steps
NULL

#' @rdname simulate_population_steps
#' @export
methods::setGeneric("simulate_population_steps",
                    signature = methods::signature("mod", "vital_rates"),
                    function(mod, vital_rates, replicates = 1,
                             carrying_capacity = NULL, initial = NULL,
                             dispersal = NULL,
                             density_dependence = NULL,
                             include_suitability = TRUE) standardGeneric("simulate_population_steps"))

#' @rdname simulate_population_steps
methods::setMethod(
  "simulate_population_steps",
  methods::signature(mod = "BiodiversityScenario", vital_rates = "matrix"),
  function(mod, vital_rates, replicates = 1, carrying_capacity = NULL,
           initial = NULL,
           dispersal = NULL, density_dependence = NULL,include_suitability = TRUE ) {
    assertthat::assert_that(
      inherits(mod, "BiodiversityScenario"),
      !is.Waiver(mod$get_predictors()),
      is.matrix(vital_rates) && all(dim(vital_rates)),
      is.numeric(replicates) && replicates >=1,
      is.null(initial) || is.Raster(initial),
      is.null(carrying_capacity) || is.Raster(carrying_capacity),
      is.logical(include_suitability),
      # Check that parameter objects are of correct class
      is.null(dispersal) || inherits(dispersal, "population_dispersal"),
      is.null(density_dependence) || inherits(density_dependence, "population_density_dependence")
    )

    # # Check that carrying capacity is set
    if(is.Raster(initial)){
      assertthat::assert_that(terra::nlyr(initial) == ncol(vital_rates),
                              msg = "Provided initial layer is not matching with the vital rates!")
      assertthat::assert_that(all( names(initial) %in% colnames(vt)),
                              msg = "Initial population raster needs to have the
                              same columns as the vital rates!")
    }

    # Check that the steps package is installed
    check_package("steps")
    if(!isNamespaceLoaded("steps")) { attachNamespace("steps");requireNamespace("steps") }

    sim <- list()

    sim[['simulation']] <- list(method = "steps",
                                params = list(
                                  "replicates" = replicates,
                                  # Vital rates
                                  "vital_rates" = vital_rates,
                                  "population" = initial,
                                  "carrying_capacity" = carrying_capacity,
                                  # Other indicators
                                  "dispersal" = dispersal,
                                  "density_dependence" = density_dependence,
                                  "include_suitability" = include_suitability
                                  )
                                )

    new <- mod$set_simulation(sim)
    return(
      bdproto(NULL, new)
    )
  }
)

# ------------------------ #

#' Internal function for the steps simulations
#'
#' @description
#' This function does the actual computation using the provided objects
#' from the projection.
#'
#' @param proj A [`SpatRaster`] object with multple timeslots
#' @param scenario_simulations A [`list`] with provided settings
#' @returns A [`SpatRaster`] object of the same length as proj.
#'
#' @noRd
#'
#' @keywords internal
.simulate_steps <- function(proj, scenario_simulations){
  assertthat::assert_that(
    is.Raster(proj),
    is.list(scenario_simulations)
  )
  # Check that the steps package is installed
  check_package("steps")
  if(!isNamespaceLoaded("steps")) { attachNamespace("steps");requireNamespace("steps") }

  # For steps to work we need to convert the input to raster
  # And also first normalize them if they somehow exceed 1
  rr <- terra::global(proj, "range",na.rm=T)
  if(!(all( rr[,1] >=0) && all(rr[,2] <=1))){
    proj <- predictor_transform(proj,option = "norm")
  }
  suit <- terra_to_raster(proj);rm(rr)

  # --- #
  # Now get the parameters if found
  params <- scenario_simulations$simulation$params
  # Get the vital rates matrix
  vt <- params$vital_rates

  # Get carrying capacity if set
  if(is.null(params$carrying_capacity)){
    carrying_capacity <- NULL

    # Create a carrying capacity by getting assuming that the total carrying
    # capacity in terms of individuals is roughly a 1/10 of the grid size.
    # FIXME: This is very much arbitrary, but having such a value supports better
    # Calculation below
    # Each size of pixel
    initial <- suit[[1]]
    carrying_capacity <- round( (initial * raster::area(initial)*.1), 0)

  } else {
    carrying_capacity <- params$carrying_capacity
    # Convert to raster for compatability
    if(is.Raster(carrying_capacity)) carrying_capacity <- terra_to_raster(carrying_capacity)
  }

  if(is.null(params$population)){
    # Get the first raster
    initial <- suit[[1]]
    population <- raster::stack(initial,initial,initial)
    names(population) <- colnames(vt)

    # Derive initial population abundances by using stable state age
    # distributions multiplied by the carrying_capacity of the landscape, and then
    # draw from a multinomial distribution to return whole integers
    stable_states <- abs( eigen(vt)$vectors[,1] / sum(eigen(vt)$vectors[,1]))
    popN <- raster::stack(replicate(ncol(vt), carrying_capacity)) * stable_states
    idx <- which(!is.na(raster::getValues(popN[[1]])))
    population <- raster::stack(
      foreach::foreach(i = 1:raster::nlayers(popN)) %do% {
        max_pop <- ceiling(raster::cellStats(popN[[i]], max, na.rm = T))
        pop_values <- popN[[i]][idx]
        popN[[i]][idx] <- stats::rbinom(prob = (pop_values/max_pop),
                                 size = max_pop,
                                 n = length(pop_values))
        popN[[i]]
      })
    names(population) <- colnames(vt)

  } else {
    population <- params[["population"]]
    # Convert to raster for compatability
    if(is.Raster(population)) population <- terra_to_raster(population)
  }

  # Set suitability to NULL if it is not be used
  if(!params$include_suitability) suit <- NULL

  # --- #
  # Define the landscape
  la <- steps::landscape(
    population = population,
    suitability = suit,
    carrying_capacity = carrying_capacity
    )
  assertthat::assert_that(inherits(la, "landscape"))

  # Define the population growth object
  pd <- steps::population_dynamics(
    change = steps::growth(transition_matrix = vt,
                           two_sex = FALSE # NOTE: Ideally capture this from the package
                           ),
    dispersal = params$dispersal,
    density_dependence = params$density_dependence
  )

  # Simulate
  sims <- try({
    steps::simulation(landscape = la,
                      population_dynamics = pd,
                      timesteps = length(terra::time(proj)),
                      replicates = params[['replicates']],
                      verbose = getOption('ibis.setupmessages'))
  },silent = F)
  if(inherits(sims,"try-error")){
    if(getOption('ibis.setupmessages')) myLog('[Scenario]','red','Simulation failed...')
    return(NULL)
  }

  # Extract and reconvert to terra
  new <- raster::stack(
    sapply(1:terra::nlyr(proj),
         function(i) steps::extract_spatial(sims, timestep = i,
                                            landscape_object = "population"))
  )
  # Convert to terra
  new <- terra::rast(new)
  terra::time(new) <- terra::time(proj)
  names(new) <- gsub("suitability", "population", names(proj))

  # A few check
  assertthat::assert_that(is.Raster(new),
                          terra::nlyr(new) == terra::nlyr(proj))
  return(new)
}
