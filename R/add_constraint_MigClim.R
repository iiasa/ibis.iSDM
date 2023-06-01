#' @include utils.R bdproto-biodiversityscenario.R
NULL

#' Add constrains to the modelled distribution projection through the MigClim approach
#'
#' @description This function adds constrain as defined by the MigClimto approach (Engler et al. 2013)
#' to a [`BiodiversityScenario-class`] object to
#' constrain future projections. For a detailed description of MigClim, please the respective reference
#' and the UserGuide. **The default parameters chosen here are suggestions.**
#' @param mod A [`BiodiversityScenario`] object with specified predictors.
#' @param rcThresholdMode A [`character`] of either **binary** or **continuous** value. Default: **continuous**
#' @param dispSteps [`numeric`] parameters on the number of dispersal steps. Dispersal steps are executed for each timestep (prediction layer)
#' and ideally should be aligned with the number of steps for projection. Minimum is \code{1} (Default) and maximum is \code{99}.
#' @param dispKernel A [`vector`] with the number of the dispersal Kernel to be applied.
#' Can be set either to a uniform numeric [vector], e.g. \code{c(1,1,1,1)} or to a proportional decline \code{(1,0.4,0.16,0.06,0.03)} (Default)
#' **Depending on the resolution of the raster, this parameter needs to be adapted**
#' @param barrierType A [character] indicating whether any set barrier should be set as \code{'strong'} or \code{'weak'} barriers.
#' Strong barriers prevent any dispersal across the barrier and weak barriers only do so if the whole [dispKernel] length
#' is covered by the barrier (Default: \code{'strong'}).
#' @param lddFreq [`numeric`] parameter indicating the frequency of long-distance dispersal (LDD) events. Default is \code{0}, so no long-distance dispersal.
#' @param lddRange A [`numeric`] value highlighting the minimum and maximum distance of LDD events.
#' **Note: The units for those distance are in cells, thus the projection units in the raster.**
#' @param iniMatAge Initial maturity age. Used together with `propaguleProd` as a proxy of population growth.
#' Must be set to the cell age in time units which are dispersal steps (Default: \code{1}).
#' @param propaguleProd Probability of a source cell to produce propagules as a function of time since colonization.
#' Set as probability vector that defines the probability of a cell producing propagules.
#' @param replicateNb Number of replicates to be used for the analysis (Default: \code{10}).
#' @param dtmp A [`character`] to a folder where temporary files are to be created.
#' @details The barrier parameter is defined through [add_barrier].
#' @seealso [`MigClim.userGuide()`]
#' @references
#' * Engler R., Hordijk W. and Guisan A. The MIGCLIM R package – seamless integration of
#' dispersal constraints into projections of species distribution models. Ecography,
#' * Robin Engler, Wim Hordijk and Loic Pellissier (2013). MigClim: Implementing dispersal
#' into species distribution models. R package version 1.6.
#' @returns Adds a MigClim onstrain to a [`BiodiversityScenario`] object.
#' @examples
#' \dontrun{
#' # Assumes that a trained 'model' object exists
#'  mod <- scenario(model) |>
#'   add_predictors(env = predictors, transform = 'scale', derivates = "none") |>
#'   add_constrain_MigClim() |>
#'   project()
#' }
#'
#' @name add_constrain_MigClim
#' @aliases add_constrain_MigClim
#' @family constrain
#' @keywords scenario
#' @exportMethod add_constrain_MigClim
#' @export
NULL
methods::setGeneric("add_constrain_MigClim",
                    signature = methods::signature("mod"),
                    function(mod, rcThresholdMode = 'continuous', dispSteps = 1,
                             dispKernel = c(1.0, 0.4, 0.16, 0.06, 0.03), barrierType = "strong",
                             lddFreq = 0, lddRange = c(1000, 10000),
                             iniMatAge = 1, propaguleProdProb = c(0.2, 0.6,0.8, 0.95),
                             replicateNb = 10, dtmp = raster::tmpDir() ) standardGeneric("add_constrain_MigClim"))

#' @name add_constrain_MigClim
#' @rdname add_constrain_MigClim
#' @usage \S4method{add_constrain_MigClim}{BiodiversityScenario, character, numeric, numeric, character, numeric, numeric, numeric, numeric, numeric, character}(mod, rcThresholdMode, dispSteps, dispKernel, barrierType, lddFreq, lddRange, iniMatAge, propaguleProdProb, replicateNb, dtmp)
methods::setMethod(
  "add_constrain_MigClim",
  methods::signature(mod = "BiodiversityScenario"),
  function(mod, rcThresholdMode = 'continuous', dispSteps = 1,
           dispKernel = c(1.0, 0.4, 0.16, 0.06, 0.03), barrierType = "strong",
           lddFreq = 0, lddRange = c(1000, 10000),
           iniMatAge = 1, propaguleProdProb = c(0.2, 0.6, 0.8, 0.95),
           replicateNb = 10, dtmp = raster::tmpDir()) {
    assertthat::assert_that(
      inherits(mod, "BiodiversityScenario"),
      !is.Waiver(mod$get_predictors()),
      is.numeric(dispSteps) && (dispSteps >=1 && dispSteps <=99),
      is.numeric(dispKernel),
      is.character(barrierType),
      is.numeric(lddFreq), is.numeric(lddRange),
      is.numeric(iniMatAge), is.numeric(propaguleProdProb)
    )

    # Parameter list
    # Note this is saved as list and not as mixed vector
    params <- list()

    # Temporary dir
    params[["dtmp"]] <- dtmp

    # Check that package is available
    check_package("MigClim")
    if(!isNamespaceLoaded("MigClim")) { attachNamespace("MigClim");requireNamespace("MigClim") }

    # Get Model
    fit <- mod$get_model()
    assertthat::assert_that(is.Raster(fit$get_data("prediction")), msg = "No prediction found in provided model.")

    # Check that initial threshold has been computed
    tr <- grep('threshold',fit$show_rasters(), ignore.case = TRUE, value = TRUE)
    assertthat::assert_that(length(tr)>0, msg = "No initial threshold has been found in the model object.")

    # Get Initial distribution object
    # iniDist <- data.frame(X = raster::coordinates(fit$get_data(tr))[,1],
    #                       Y = raster::coordinates(fit$get_data(tr))[,2],
    #                       value = values(fit$get_data(tr)))
    # Alternatively save to temporary folder and pass pathname
    dir.create(dtmp, showWarnings = FALSE)
    r <- fit$get_data(tr)
    if(is.factor(r)) r <- raster::deratify(r, complete = TRUE)
    suppressWarnings(
      raster::writeRaster(x = r, filename = file.path(dtmp, paste0(tr, ".tif")),
                          dt = "INT2S", varNA = -9999, prj = TRUE, overwrite = TRUE)
    )
    params[["iniDist"]] <- file.path(dtmp, paste0(tr, ".tif")); rm(r)

    # Get suitability map
    params[["hsMap"]] <- file.path(dtmp,"SuitabilityProjection") # The basename to be used. Each projection will have an incremental number added here

    # Habitat suitability threshold data
    # Set to 0 for continous and to 1-1000 for binary (pick 750 as default)
    params[["rcThreshold"]] <- match.arg(rcThresholdMode, choices = c("binary", "continuous"), several.ok = FALSE)

    # Dispersal steps
    params[["dispSteps"]] <- dispSteps

    # Dispersal kernel definition (expressed as probability)
    assertthat::assert_that(all(dispKernel >=0 && dispKernel <=1))
    params[["dispKernel"]] <- dispKernel

    # Barrier
    params[["barrierType"]] <- match.arg(barrierType, choices = c("weak", "strong"), several.ok = FALSE)

    # LDD params
    params[["lddFreq"]] <- lddFreq
    assertthat::assert_that( is.null(lddRange) || (lddRange[2] > lddRange[1]))
    params[["lddMinDist"]] <- lddRange[1]; params[["lddMaxDist"]] <- lddRange[2]

    # Maturity age and propagule production probability
    params[["iniMatAge"]] <- iniMatAge
    params[["propaguleProdProb"]]  <- propaguleProdProb

    # Simulation name
    params[["simulName"]] <- paste0(raster::tmpDir(), "MigClim_", fit$model$runname )

    # Replicate number
    params[["replicateNb"]] <- replicateNb

    params[["fullOutput"]] <- FALSE # Don't create full outputs to be condense
    params[["overwrite"]] <- TRUE # Always overwrite
    params[["keepTempFiles"]] <- FALSE

    # Add to scenario object
    cr <- list()
    cr[['dispersal']] <- list(method = "MigClim",
                              params = params)

    new <- mod$set_constraints(cr)
    return(
      bdproto(NULL, new)
    )
  }
)
