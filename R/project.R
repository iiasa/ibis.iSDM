#' @include utils.R bdproto-biodiversityscenario.R
NULL

#' Project a fitted model to a new environment and covariates
#'
#' @description
#' Equivalent to [train], this function acts as a
#' wrapper to project the model stored in a [`BiodiversityScenario-class`] object to
#' newly supplied (future) covariates. Supplied predictors are usually spatial-temporal predictors
#' which should be prepared via [`add_predictors()`] (e.g. transformations and derivates) in the same way as they have been during
#' the initial modelling with [`distribution()`].
#' Any constrains specified in the scenario object are applied during the projection.
#'
#' @details
#' In the background the function \code{x$project()} for the respective model object is called, where
#' \code{x} is fitted model object. For specifics on the constrains, see the relevant [constrain] functions,
#' respectively:
#' * [`add_constrain()`] for generic wrapper to add any of the available constrains.
#' * [`add_constrain_dispersal()`] for specifying dispersal constrain on the temporal projections at each step.
#' * [`add_constrain_MigClim()`] Using the \pkg{MigClim} R-package to simulate dispersal in projections.
#' * [`add_constrain_connectivity()`] Apply a connectivity constrain at the projection, for instance by adding
#' a barrier that prevents migration.
#' * [`add_constrain_adaptability()`] Apply an adaptability constrain to the projection, for instance
#' constraining the speed a species is able to adapt to new conditions.
#' * [`add_constrain_boundary()`] To artificially limit the distribution change. Similar as specifying projection limits, but
#' can be used to specifically constrain a projection within a certain area (e.g. a species range or an island).
#'
#' Many constrains also requires thresholds to be calculated. Adding [`threshold()`] to a
#' [`BiodiversityScenario-class`] object enables the computation of thresholds at every step based on the threshold
#' used for the main model (threshold values are taken from there).
#'
#' Finally this function also allows temporal stabilization across prediction steps via enabling
#' the parameter \code{stabilize} and checking the \code{stablize_method} argument. Stabilization can for instance
#' be helpful in situations where environmental variables are quite dynamic, but changes in projected suitability
#' are not expected to abruptly increase or decrease. It is thus a way to smoothen out outliers from the projection.
#' Options are so far for instance \code{'loess'} which fits a [`loess()`] model per pixel and time step. This is conducted at
#' the very of the processing steps and any thresholds will be recalculated afterwards.
#'
#' @seealso [`scenario()`]
#' @param mod A [`BiodiversityScenario`] object with set predictors.
#' Note that some constrains such as [MigClim] can still simulate future change without projections.
#' @param date_interpolation A [`character`] on whether dates should be interpolated. Options
#' include \code{"none"} (Default), \code{"annual"}, \code{"monthly"}, \code{"daily"}.
#' @param stabilize A [`boolean`] value indicating whether the suitability projection should be stabilized (Default: \code{FALSE}).
#' @param stabilize_method [`character`] stating the stabilization method to be applied. Currently supported is \code{`loess`}.
#' @param layer A [`character`] specifying the layer to be projected (Default: \code{"mean"}).
#' @param ... passed on parameters.
#' @returns Saves [`stars`] objects of the obtained predictions in mod.
#'
#' @name project
#' @aliases project
#' @keywords scenarios
#' @exportMethod project
#' @export
NULL
methods::setGeneric("project",
                    signature = methods::signature("mod"),
                    function(mod, date_interpolation = "none", stabilize = FALSE, stabilize_method = "loess",
                             layer = "mean", ...) standardGeneric("project"))

#' @name project
#' @rdname project
#' @usage \S4method{project}{BiodiversityScenario, character, logical, character, character}(mod, date_interpolation, stabilize, stabilize_method, layer)
methods::setMethod(
  "project",
  methods::signature(mod = "BiodiversityScenario"),
  function(mod, date_interpolation = "none", stabilize = FALSE, stabilize_method = "loess",
           layer = "mean", ...){
    # date_interpolation = "none"; stabilize = FALSE; stabilize_method = "loess"; layer="mean"
    assertthat::assert_that(
      inherits(mod, "BiodiversityScenario"),
      !is.Waiver(mod$get_predictors()),
      is.character(date_interpolation),
      is.logical(stabilize),
      is.character(layer)
    )
    # Match methods
    date_interpolation <- match.arg(date_interpolation, c("none", "yearly", "annual", "monthly", "daily"), several.ok = FALSE)
    stabilize_method <- match.arg(stabilize_method, c("loess"), several.ok = FALSE)
    if(!is.Waiver(mod$get_data())) if(getOption('ibis.setupmessages')) myLog('[Scenario]','red','Overwriting existing scenarios...')

    # Get the model object
    fit <- mod$get_model()
    # Check that coefficients and model exist
    assertthat::assert_that(!is.Waiver(fit),
                            nrow(fit$get_coefficients())>0,
                            msg = "No model or coefficients found!")
    # Get predictors
    new_preds <- mod$get_predictors()
    if(is.Waiver(new_preds)) stop('No scenario predictors found.')
    new_crs <- new_preds$get_projection()
    if(is.na(new_crs)) if(getOption('ibis.setupmessages')) myLog('[Scenario]','yellow','Missing projection of future predictors.')

    # Interpolate dates if set
    if(date_interpolation!="none"){
      if(getOption('ibis.setupmessages')) myLog('[Scenario]','green',paste0('Interpolating dates for scenario predictors as: ', date_interpolation))
      o <- approximate_gaps(env = new_preds$get_data(), date_interpolation = date_interpolation)
      # Set new data
      #new_preds$set_data()
    }

    # Get limits if present
    if(!is.null( mod$get_limits() )){
      # FIXME: Scenarios need to be checked that the right layer is taken!!
      # Get prediction
      n <- fit$show_rasters()[grep("threshold",fit$show_rasters())]
      tr <- fit$get_data(n)[[1]]
      tr <- cbind( terra::crds(tr), data.frame(thresh = values(tr)))
      tr[['thresh']] <- ifelse(tr[['thresh']]==0, NA, tr[['thresh']])
      tr <- tr |> (\(.) subset(., stats::complete.cases(thresh)))()

      # Get zones from the limiting area, e.g. those intersecting with input
      suppressMessages(
        suppressWarnings(
          zones <- st_intersection(sf::st_as_sf(tr, coords = c('x','y'), crs = sf::st_crs(fit$model$background)),
                                   mod$get_limits()
          )
        )
      )
      # Limit zones
      zones <- subset(mod$get_limits(), limit %in% unique(zones$limit) )
      # Now clip all provided new predictors and background to this
      new_preds$crop_data(zones)
    }

    # Check that predictor names are all present
    mod_pred_names <- fit$model$predictors_names
    pred_names <- mod$get_predictor_names()
    assertthat::assert_that( all(mod_pred_names %in% pred_names),
                             msg = paste0('Model predictors are missing from the scenario predictor!') )

    # Get constraints, threshold values and other parameters
    scenario_threshold <- mod$get_threshold()
    # Not get the baseline raster
    thresh_reference <- grep('threshold',fit$show_rasters(),value = T)[1] # Use the first one always
    baseline_threshold <- mod$get_model()$get_data(thresh_reference)
    if(!is.Waiver(scenario_threshold)){
      if(is.na(sf::st_crs(baseline_threshold))) sf::st_crs(baseline_threshold) <- sf::st_crs( fit$model$background )
    }

    if(inherits(baseline_threshold, 'SpatRaster')){
      baseline_threshold <- baseline_threshold[[grep(layer, names(baseline_threshold))]]
    }
    scenario_constraints <- mod$get_constraints()

    #  --- Check that everything is there ---
    # Check that thresholds are set for constrains
    if("dispersal" %in% names(scenario_constraints)){
      if(scenario_constraints[["dispersal"]]$method == "MigClim") {
        assertthat::assert_that(is.Raster(baseline_threshold))
      } else if(scenario_constraints[["dispersal"]]$method == "kissmig"){
        assertthat::assert_that( is.Raster(baseline_threshold))
        if(!is.Waiver(scenario_threshold)) {
          if(getOption('ibis.setupmessages')) myLog('[Scenario]','green','Using kissmig to calculate updated distribution thresholds.')
          scenario_threshold <- new_waiver()
        }
      } else {
        assertthat::assert_that(!is.Waiver(scenario_threshold),msg = "Other constrains require threshold option!")
      }
    }
    if("connectivity" %in% names(scenario_constraints) && "dispersal" %notin% names(scenario_constraints)){
      if(getOption('ibis.setupmessages')) myLog('[Scenario]','yellow','Connectivity contraints make most sense with a dispersal constraint.')
    }
    # ----------------------------- #
    #   Start of projection         #
    # ----------------------------- #

    # Now convert to data.frame and subset
    df <- new_preds$get_data(df = TRUE)
    names(df)[1:3] <- tolower(names(df)[1:3]) # Assuming the first three attributes are x,y,t
    assertthat::assert_that(nrow(df)>0,
                            utils::hasName(df,'x'), utils::hasName(df,'y'), utils::hasName(df,'time'),
                            msg = "Error: Projection data and training data are not of equal size and format!")
    df <- subset(df, select = c("x", "y", "time", mod_pred_names) )
    df$time <- to_POSIXct( df$time )
    # Convert all units classes to numeric or character to avoid problems
    df <- units::drop_units(df)

    # ------------------ #
    if(getOption('ibis.setupmessages')) myLog('[Scenario]','green','Starting suitability projections for ', length(unique(df$time)), ' timesteps.')

    # Now for each unique element, loop and project in order
    proj <- terra::rast()
    proj_thresh <- terra::rast()

    pb <- progress::progress_bar$new(format = "Creating projections (:spin) [:bar] :percent",
                                     total = length(unique(df$time)))
    # TODO: Consider doing this in parallel but sequential
    times <- sort(unique(df$time))
    for(step in times){
      # Get data
      nd <- subset(df, time == step)

      # Apply adaptability constrain
      if("adaptability" %in% names(scenario_constraints)){
        if(scenario_constraints[["adaptability"]]$method == "nichelimit") {
          nd <- .nichelimit(newdata = nd, model = mod$get_model()[['model']],
                           names = scenario_constraints[["adaptability"]]$params['names'],
                           value = scenario_constraints[["adaptability"]]$params['value'],
                           increment = scenario_constraints[["adaptability"]]$params['increment'],
                           increment_step = which(step==times) )
        }
      }

      # Project suitability
      out <- fit$project(newdata = nd, layer = layer)
      names(out) <- paste0("suitability", "_", layer, "_", step)
      if(is.na(sf::st_crs(out))) sf::st_crs(out) <- sf::st_crs( fit$model$background )

      # If other constrains are set, apply them posthoc
      if(!is.Waiver(scenario_constraints)){
        # Apply a resistance surface if found
        if("connectivity" %in% names(scenario_constraints)){
          # Get the layer for later
          resistance <- scenario_constraints$connectivity$params$resistance
          # By definition a hard barrier removes all suitable again
          if(any(scenario_constraints$connectivity$method == "resistance")){
            if(terra::nlyr(resistance)>1){
              ind <- which( terra::time(resistance) == as.Date(step) ) # Get specific step
              assertthat::assert_that(is.numeric(ind))
              resistance <- resistance[[ind]]
            }
            out <- out * resistance
          }
        } else {
          resistance <- NULL
        }

        # Calculate dispersal constraint if set
        if("dispersal" %in% names(scenario_constraints) ){
          # MigClim simulations are run posthoc
          if(scenario_constraints$dispersal$method %in% c("sdd_fixed", "sdd_nexpkernel")){
            out <- switch (scenario_constraints$dispersal$method,
                           "sdd_fixed" = .sdd_fixed(baseline_threshold, out,
                                                    value = scenario_constraints$dispersal$params[1],
                                                    resistance = resistance ),
                           "sdd_nexpkernel" = .sdd_nexpkernel(baseline_threshold, out,
                                                              value = scenario_constraints$dispersal$params[1],
                                                              resistance = resistance)
            )
            names(out) <-  paste0('suitability_', step)
          }
          # For kissmig generate threshold and masked suitabilities
          if(scenario_constraints$dispersal$method == "kissmig"){
            out <- .kissmig_dispersal(baseline_threshold,
                                      new_suit = out,
                                      resistance = resistance,
                                      params = scenario_constraints$dispersal$params)
            # Returns a layer of two with both the simulated threshold and the masked suitability raster
            names(out) <- paste0(c('threshold_', 'suitability_'), step)
            # Add threshold to result stack
            suppressWarnings( proj_thresh <- c(proj_thresh, out[[1]] ) )
            baseline_threshold <- out[[1]]
            out <- out[[2]]
          }
        }

        # Connectivity constraints with hard barriers
        if("connectivity" %in% names(scenario_constraints)){
          # By definition a hard barrier removes all suitable again
          if(any(scenario_constraints$connectivity$method == "hardbarrier")){
            out[resistance==1] <- 0
          }
        }

      }

      # Recalculate thresholds if set manually
      if(!is.Waiver(scenario_threshold)){
        # FIXME: Currently this works only for mean thresholds. Think of how the other are to be handled
        scenario_threshold <- scenario_threshold[1]
        out_thresh <- out
        out_thresh[out_thresh < scenario_threshold] <- 0; out_thresh[out_thresh >= scenario_threshold] <- 1
        names(out_thresh) <-  paste0('threshold_', step)
        # If threshold is
        if( terra::global(out_thresh, 'max', na.rm = TRUE) == 0){
          if(getOption('ibis.setupmessages')) myLog('[Scenario]','yellow','Thresholding removed all grid cells. Using last years threshold.')
          out_thresh <- baseline_threshold
        } else { baseline_threshold <- out_thresh }
        # Add to result stack
        suppressWarnings( proj_thresh <- c(proj_thresh, out_thresh) )
      }
      # Add to result stack
      suppressWarnings( proj <- c(proj, out) )
      pb$tick()
    }
    rm(pb)
    terra::time(proj) <- times
    if(terra::nlyr(proj_thresh)>1) terra::time(proj_thresh) <- times

    # Apply MigClim and other post-hoc constraints if set
    # FIXME: Ideally make this whole setup more modular. So create suitability projections first
    if(!is.Waiver(scenario_constraints)){
      # Calculate dispersal constraint if set
      if("dispersal" %in% names(scenario_constraints) ){
        # MigClim simulations are run posthoc
        if(scenario_constraints$dispersal$method == 'MigClim'){
          # Get Parameters
          params <- scenario_constraints$dispersal$params

          pb <- progress::progress_bar$new(total = terra::nlyr(proj))
          for(lyr in 1:terra::nlyr(proj)){
            pb$tick()
            # Normalize the projected suitability rasters to be in range 0-1000 and save
            hsMap <- predictor_transform(env = proj[[lyr]], option = "norm") * 1000
            # Write as filename in the destined folder
            suppressWarnings(
              terra::writeRaster(x = hsMap, filename = paste0( params[["hsMap"]],lyr,".tif"),
                                  datatype = "INT2S", NAflag = -9999, overwrite = TRUE)
              )
            rm(hsMap)
          };rm(pb)
          # For threshold, define based on type
          tr <- ifelse(params[["rcThreshold"]] == "continuous", 0, 750) # Set to 75% as suitability threshold

          # Now run MigClim by switching to temporary dir
          dir.ori <- getwd()
          setwd(params[["dtmp"]])
          try({
            m <- MigClim::MigClim.migrate(
              iniDist = basename(params[["iniDist"]]),
              hsMap = basename(params[["hsMap"]]),
              rcThreshold = tr,
              envChgSteps = terra::nlyr(proj), # Use number of projected suitability layers
              dispSteps = params[["dispSteps"]],
              dispKernel = params[["dispKernel"]],
              barrier = "", # TBD. Loaded via another arguement
              barrierType = params[["barrierType"]],
              iniMatAge = params[["iniMatAge"]],
              propaguleProd = params[["propaguleProdProb"]],
              lddFreq = params[["lddFreq"]],
              lddMinDist = params[["lddMinDist"]], lddMaxDist = params[["lddMaxDist"]],
              simulName = basename(params[["simulName"]]),
              replicateNb = params[["replicateNb"]],
              overWrite = params[["overwrite"]],
              testMode = FALSE,
              fullOutput = params[["fullOutput"]],
              keepTempFiles = params[["keepTempFiles"]]
            )
          })
          # Get average stats
          run_stats <- utils::read.table(
            file.path(basename(params[["simulName"]]), paste0(basename(params[["simulName"]]),"_stats.txt")),
                                  header = TRUE
            )
          run_sums <- utils::read.table(
            file.path(basename(params[["simulName"]]), paste0(basename(params[["simulName"]]),"_summary.txt")),
            header = TRUE
          )
          # Get MigClim outputs
          ll <- list.files(basename(params[["simulName"]]),'asc',full.names = TRUE)
          run_sims <- terra::rast(ll); names(run_sims) <- tools::file_path_sans_ext(basename(ll))
          # Condense the simulation runs into one modal prediction
          run_sim <- terra::app(run_sims, terra::modal)
          sf::st_crs(run_sim) <- sf::st_crs(fit$get_data('prediction'))
          all(sapply(list.files(getwd(),".tif"), file.remove)) # Cleanup
          setwd(dir.ori) # Flip back to original directory

          # Wrap all results in a list
          mc <- list(params = params,
                     stats = run_stats,
                     summary = run_sums,
                     raster = run_sim)
        }
      } # End of MigClim processing chain
    }
    # If not found, set a waiver
    if(!exists("mc")) mc <- new_waiver()

    # ---- #
    assertthat::assert_that(
      is.Raster(proj), is.Raster(proj_thresh),
      msg = "Something went wrong with the projection."
    )

    # Apply boundary constraints if set
    if("boundary" %in% names(scenario_constraints)){
      if(!terra::compareGeom(proj, scenario_constraints$boundary$params$layer, stopOnError = FALSE)){
        scenario_constraints$boundary$params$layer <- alignRasters(
          scenario_constraints$boundary$params$layer,
          proj,
          method = "near", func = terra::modal, cl = FALSE
        )
      }
      proj <- terra::mask(proj, scenario_constraints$boundary$params$layer)
      # Get background and ensure that all values outside are set to 0
      proj[is.na(proj)] <- 0
      proj <- terra::mask(proj, fit$model$background )
      # Also for thresholds if existing
      if(terra::nlyr(proj_thresh)>0){
        proj_thresh <- terra::mask(proj_thresh, scenario_constraints$boundary$params$layer)
        proj_thresh[is.na(proj_thresh)] <- 0
        proj_thresh <- terra::mask(proj_thresh, fit$model$background )
      }
    }

    # Should stabilization be applied?
    if(stabilize){
      if(getOption('ibis.setupmessages')) myLog('[Scenario]','green','Applying stabilization.')
      if(stabilize_method == "loess"){
        # FIXME: Could outsource this code
        impute.loess <- function(y, x.length = NULL, s = 0.75,
                                 smooth = TRUE, na.rm, ...) {
          if (is.null(x.length)) {
            x.length = length(y)
          }
          if(length(y[!is.na(y)]) < 8) {
            y <- rep(NA, x.length)
          } else {
            x <- 1:x.length
            p <- suppressWarnings( stats::loess(y ~ x, span = s,
                                                data.frame(x = x, y = y)) )
            if (smooth == TRUE) {
              y <- stats::predict(p, x)
            } else {
              na.idx <- which(is.na(y))
              if (length(na.idx) > 1) {
                y[na.idx] <- stats::predict(p, data.frame(x = na.idx))
              }
            }
          }
          return(y)
        }
        new_proj <- terra::lapp(proj, fun = impute.loess)
        # Rename again
        names(new_proj) <- names(proj)
        terra::time(new_proj) <- times
        proj <- new_proj; rm(new_proj)
        # Were thresholds calculated? If yes, recalculate on the smoothed estimates
        if(terra::nlyr(proj_thresh)>0){
          new_thresh <- proj
          new_thresh[new_thresh < scenario_threshold[1]] <- 0
          new_thresh[new_thresh >= scenario_threshold[1]] <- 1
          names(new_thresh) <- names(proj_thresh)
          thresh <- new_thresh; rm(new_thresh)
        }
      }
    }

    # Finally convert to stars and rename
    proj <- stars::st_as_stars(proj,
                               crs = sf::st_crs(new_crs)
    ); names(proj) <- 'suitability'

    if(terra::nlyr(proj_thresh)>0){
      # Add the thresholded maps as well
      proj_thresh <- stars::st_as_stars(proj_thresh,
                                       crs = sf::st_crs(new_crs)
      ); names(proj_thresh) <- 'threshold'
      # Correct band if different
      if(all(!stars::st_get_dimension_values(proj, 3) != stars::st_get_dimension_values(proj_thresh, 3 ))){
        proj_thresh <- stars::st_set_dimensions(proj_thresh, 3, values = stars::st_get_dimension_values(proj, 3))
      }
      proj <- stars:::c.stars(proj, proj_thresh)
    }

    # Return output by adding it to the scenario object
    bdproto(NULL, mod,
            scenarios = proj,
            scenarios_migclim = mc
            )
  }
)
