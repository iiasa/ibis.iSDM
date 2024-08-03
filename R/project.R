#' @include class-biodiversityscenario.R
NULL

#' Project a fitted model to a new environment and covariates
#'
#' @description Equivalent to [train], this function acts as a wrapper to
#' project the model stored in a [`BiodiversityScenario-class`] object to newly
#' supplied (future) covariates. Supplied predictors are usually
#' spatial-temporal predictors which should be prepared via [`add_predictors()`]
#' (e.g. transformations and derivates) in the same way as they have been during
#' the initial modelling with [`distribution()`]. Any constrains specified in
#' the scenario object are applied during the projection.
#'
#' @param x A [`BiodiversityScenario`] object with set predictors. Note that some
#' constrains such as \code{MigClim} can still simulate future change without projections.
#' @param date_interpolation A [`character`] on whether dates should be interpolated.
#' Options include \code{"none"} (Default), \code{"annual"}, \code{"monthly"}, \code{"daily"}.
#' @param stabilize A [`logical`] value indicating whether the suitability projection
#' should be stabilized (Default: \code{FALSE}).
#' @param stabilize_method [`character`] stating the stabilization method to be
#' applied. Currently supported is \code{`loess`}.
#' @param layer A [`character`] specifying the layer to be projected (Default: \code{"mean"}).
#' @param verbose Setting this [`logical`] value to \code{TRUE} prints out further
#' information during the model fitting (Default: \code{FALSE}).
#' @param ... passed on parameters.
#'
#' @details In the background the function \code{x$project()} for the respective
#' model object is called, where \code{x} is fitted model object. For specifics
#' on the constraints, see the relevant \code{constrain} functions, respectively:
#'
#' * [`add_constraint()`] for generic wrapper to add any of the available constrains.
#' * [`add_constraint_dispersal()`] for specifying dispersal constraint on the temporal projections at each step.
#' * [`add_constraint_MigClim()`] Using the \pkg{MigClim} R-package to simulate dispersal in projections.
#' * [`add_constraint_connectivity()`] Apply a connectivity constraint at the projection, for instance by adding
#' a barrier that prevents migration.
#' * [`add_constraint_minsize()`] Adds a constraint on the minimum area a given
#' thresholded patch should have, assuming that smaller areas are in fact not
#' suitable.
#' * [`add_constraint_adaptability()`] Apply an adaptability constraint to the projection, for instance
#' constraining the speed a species is able to adapt to new conditions.
#' * [`add_constraint_boundary()`] To artificially limit the distribution change. Similar as specifying projection limits, but
#' can be used to specifically constrain a projection within a certain area
#' (e.g. a species range or an island).
#'
#' Many constrains also requires thresholds to be calculated. Adding
#' [`threshold()`] to a [`BiodiversityScenario-class`] object enables the
#' computation of thresholds at every step based on the threshold used for the
#' main model (threshold values are taken from there).
#'
#' It is also possible to make a complementary simulation with the \code{steps}
#' package, which can be provided via [`simulate_population_steps()`] to the
#' [`BiodiversityScenario-class`] object. Similar as with thresholds, estimates
#' values will then be added to the outputs.
#'
#' Finally this function also allows temporal stabilization across prediction
#' steps via enabling the parameter \code{stabilize} and checking the
#' \code{stablize_method} argument. Stabilization can for instance be helpful in
#' situations where environmental variables are quite dynamic, but changes in
#' projected suitability are not expected to abruptly increase or decrease. It
#' is thus a way to smoothen out outliers from the projection. Options are so
#' far for instance \code{'loess'} which fits a [`loess()`] model per pixel and
#' time step. This is conducted at the very of the processing steps and any
#' thresholds will be recalculated afterwards.
#'
#' @returns Saves [`stars`] objects of the obtained predictions in mod.
#'
#' @seealso [`scenario()`]
#' @keywords scenarios
#'
#' @examples
#' \dontrun{
#' # Fit a model
#' fit <- distribution(background) |>
#'         add_biodiversity_poipa(surveydata) |>
#'         add_predictors(env = predictors) |>
#'         engine_breg() |>
#'         train()
#'
#' # Fit a scenario
#' sc <- scenario(fit) |>
#'         add_predictors(env = future_predictors) |>
#'         project()
#' }
#'
#' @import terra
#'
#' @name project
NULL

#' @rdname project
#' @export
project.BiodiversityScenario <- function(x,...) project(x,...)

#' @rdname project
#' @export
methods::setMethod(
  "project",
  methods::signature(x = "BiodiversityScenario"),
  function(x, date_interpolation = "none", stabilize = FALSE, stabilize_method = "loess",
           layer = "mean", verbose = getOption('ibis.setupmessages', default = TRUE), ...){
    # MJ: Workaround to ensure project generic does not conflict with terra::project
    mod <- x
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
    if(!is.Waiver(mod$get_data())) if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Scenario]','red','Overwriting existing scenarios...')

    # Get the model object
    fit <- mod$get_model(copy = TRUE)
    # Get background
    background <- fit$model$background

    assertthat::assert_that(!is.na(terra::crs( background )),
                            msg = "Model background has no CRS which will only cause problems. Rerun and set!")

    # Check that coefficients and model exist
    assertthat::assert_that(!is.Waiver(fit),
                            nrow(fit$get_coefficients())>0,
                            msg = "No model or coefficients found!")
    # Get predictors
    new_preds <- mod$get_predictors()
    if(is.Waiver(new_preds)) stop('No scenario predictors found.')

    new_crs <- new_preds$get_projection()
    if(is.na(new_crs)) if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Scenario]','yellow','Missing projection of future predictors.')

    # Interpolate predictor set if specified
    if(date_interpolation!="none"){
      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Scenario]','green',paste0('Interpolating dates for scenario predictors as: ', date_interpolation))

      new <- interpolate_gaps(new_preds$get_data(),
                              date_interpolation = date_interpolation)
      # Set new data
      new_preds <- new_preds$set_data(new)
    }

    # Check latent factors and add them to future if found
    if(!is.Waiver(mod$latentfactors) || fit$has_latent()){
      # Get predictor names
      pn <- fit$model$predictors_names
      ind <- grep("nearestpoint_|spatialtrend_|kde__coordinates", pn,value = TRUE)
      if(length(mod$latentfactors)==2){
        # Get the layer if found
        if(!is.null(mod$latentfactors$layer) && length(ind)>0){
          layer <- mod$latentfactors$layer
          new_preds <- new_preds$set_data(layer)
        } else if(length(ind)>0){
          # Get the predictors
          sps <- fit$model$predictors_object$get_data() |> terra::subset(ind)
          new_preds <- new_preds$set_data(sps)
        }
        if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Scenario]','green', "Found latent factors and added them to future predictors stack.")
      } else {
        if(length(ind)>0){
          # Get the predictors
          sps <- fit$model$predictors_object$get_data() |> terra::subset(ind)
          new_preds <- new_preds$set_data(sps)
        }
        if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Scenario]','green', "Found latent factors and added them to future predictors stack.")
      }
      suppressWarnings( try({rm(pn,ind,sps)},silent = TRUE) )
    }

    # Get limits if flagged as present in scenario object
    if(!is.null( mod$get_limits() )){
      # Get Limit and settings from model
      limits <- mod$get_limits()

      # Get Settings object from model to obtain
      settings <- fit$settings
      # Get zone limits form the model
      nrs <- settings$get("limits_zones_categories") |> as.numeric()
      if(length(nrs)>0 && utils::hasName(limits, "limit")){
        limits <- limits |> dplyr::filter(limit %in% nrs)
      }

      # Clip new predictors
      # This also takes account of nearest time zones
      new_preds$crop_data(limits,
                          apply_time = ifelse(any(limits$time == "No set time"), FALSE, TRUE))

      # Filter to first time slot if present for remaining clipping
      if(utils::hasName(limits,"time")){
        if(dplyr::n_distinct(limits$time)>1){
          if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Scenario]','green', "Using first set temporal slot for zoning!")
          first <- sort(limits$time)[1]
          limit <- limits |> dplyr::filter(time == first )
          rm(first)
        }
      }
      # MJ: This code below is necessary for some engine predictions
      # However by default is not done unless limits_clip is set to TRUE.
      # Also adjust the model container
      fit$model$background <- limits
      # Clip the predictor object
      fit$model$predictors_object <- fit$model$predictors_object$clone(deep = TRUE)
      fit$model$predictors_object$crop_data(limits)
      fit$model$predictors_object$mask(limits)
      fit$model$predictors <- fit$model$predictors_object$get_data(df = TRUE, na.rm = FALSE)
      # And offset if found
      if(!is.null(fit$model$offset_object)){
        fit$model$offset_object <- terra::deepcopy(fit$model$offset_object)
        fit$model$offset_object <- terra::crop(fit$model$offset_object, limits)
        fit$model$offset_object <- terra::mask(fit$model$offset_object, limits)
        fit$model$offset <- terra::as.data.frame(fit$model$offset_object, xy = TRUE, na.rm = FALSE)
      }
    }

    # Check that model has any (?) coefficients
    assertthat::assert_that( nrow(fit$get_coefficients())>0,
                             msg = paste0('No coefficients found in model.') )

    # Check that predictor names are all present.
    mod_pred_names <- fit$model$predictors_names
    pred_names <- mod$get_predictor_names()
    # get predictor names of integrated model
    if (!is.Waiver(fit$.internals)) {
      int_pred_names <- lapply(fit$.internals, function(i) i$model$model$predictors_names)
    } else { int_pred_names <- NULL }

    if (is.Waiver(fit$.internals)) {
      assertthat::assert_that(all(mod_pred_names %in% pred_names),
                              msg = paste0('Model predictors are missing from the scenario predictor!'))
    } else {
      # check if missing preds are in .internals list
      assertthat::assert_that(sum(!mod_pred_names %in% pred_names) == 1 ||
                                any(!unlist(int_pred_names, use.names = FALSE) %in% pred_names),
                              msg = "Model predictors are missing from the scenario predictor!")
    }

    # Get constraints, threshold values and other parameters
    scenario_threshold <- mod$get_threshold()
    if(!is.Waiver(scenario_threshold)){
      # Not get the baseline raster
      thresh_reference <- grep('threshold',fit$show_rasters(),value = T)[1] # Use the first one always
      assertthat::assert_that(!is.na(thresh_reference))
      baseline_threshold <- mod$get_model()$get_data(thresh_reference)

      if(is.na(terra::crs(baseline_threshold))) terra::crs(baseline_threshold) <- terra::crs( background )
      # Furthermore apply new limits also to existing predictions (again)
      if(!is.null( mod$get_limits() )){
        # Get Limit and settings from model
        limits <- mod$get_limits()
        # Filter to first time slot if present for remaining clipping
        if(utils::hasName(limits,"time")){
          if(dplyr::n_distinct(limits$time)>1){
            if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Scenario]','green', "Using first set temporal slot for zoning!")
            first <- sort(limits$time)[1]
            limit <- limits |> dplyr::filter(time == first )
            rm(first)
          }
        }
        baseline_threshold <- terra::crop(baseline_threshold, limits)
      }
      if(inherits(baseline_threshold, 'SpatRaster')){
        baseline_threshold <- baseline_threshold[[grep(layer, names(baseline_threshold))]]
      }

      # Set all NA values to 0 (and then mask by background?)
      baseline_threshold[is.na(baseline_threshold)]<-0

    } else {
      baseline_threshold <- new_waiver()
    }

    # Optional constraints or simulations if specified
    scenario_constraints <- mod$get_constraints()
    scenario_simulations <- mod$get_simulation()

    # Create a template for use
    template <- emptyraster( new_preds$get_data() )

    #  --- Check that everything is there ---
    # Check that thresholds are set for constrains
    if("dispersal" %in% names(scenario_constraints)){
      if(scenario_constraints[["dispersal"]]$method == "MigClim") {
        assertthat::assert_that(is.Raster(baseline_threshold))
      } else if(scenario_constraints[["dispersal"]]$method == "kissmig"){
        assertthat::assert_that( is.Raster(baseline_threshold))
        if(!is.Waiver(scenario_threshold)) {
          if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Scenario]','green','Using kissmig to calculate updated distribution thresholds.')
          scenario_threshold <- new_waiver()
        }
      } else {
        assertthat::assert_that(!is.Waiver(scenario_threshold),msg = "Other constrains require threshold option!")
      }
    }

    if(("threshold" %in% names(scenario_constraints)) && is.Waiver(baseline_threshold)){
      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Scenario]','yellow','Threshold constraint found but not threshold set? Apply threshold()!')
    }

    if("connectivity" %in% names(scenario_constraints) && "dispersal" %notin% names(scenario_constraints)){
      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Scenario]','yellow','Connectivity contraints make most sense with a dispersal constraint.')
    }
    # ----------------------------- #
    ####   Start of projection   ####
    # ----------------------------- #

    # Now convert to data.frame and subset
    df <- new_preds$get_data(df = TRUE)

    # cell ids are need to match integrated prediction to future pred
    if (!is.Waiver(fit$.internals)) {
      id_tmp <- terra::extract(x = fit$get_data(), y = df[,1:2], cells = TRUE, layer = 1)
      df$cell <- id_tmp$cell # this should work because order stays the same
    }

    names(df)[1:3] <- tolower(names(df)[1:3]) # Assuming the first three attributes are x,y,t
    assertthat::assert_that(nrow(df)>0,
                            utils::hasName(df,'x'), utils::hasName(df,'y'), utils::hasName(df,'time'),
                            msg = "Error: Projection data and training data are not of equal size and format!")

    df <- dplyr::select(df, dplyr::any_of(c("x", "y", "cell", "time",
                                            unlist(int_pred_names, use.names = FALSE),
                                            mod_pred_names)))

    df$time <- to_POSIXct(df$time)
    # Convert all units classes to numeric or character to avoid problems
    df <- units::drop_units(df)

    # ------------------ #
    if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Scenario]','green','Starting suitability projections for ',
                                              length(unique(df$time)), ' timesteps from ', paste(range(df$time), collapse = " <> "))

    # Now for each unique element, loop and project in order
    proj <- terra::rast()
    proj_thresh <- terra::rast()

    if(getOption('ibis.setupmessages', default = TRUE)){
      pb <- progress::progress_bar$new(format = "Creating projections (:spin) [:bar] :percent",
                                       total = length(unique(df$time)))
    }

    # TODO: Consider doing this in parallel but sequential
    times <- sort(unique(df$time))

    for(step in times){

      # Get data
      nd <- subset(df, time == step)

      # check that timestep has data
      assertthat::assert_that(nrow(nd)>0, !all(is.na(dplyr::select(nd, dplyr::any_of(mod_pred_names)))))

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

      # loop through integrated models
      if (!is.Waiver(fit$.internals)) {

        # loop through all internals
        for (i in 1:length(fit$.internals)) {

          # get current predictors and project
          pred_tmp <- c("x", "y", fit$.internals[[i]]$model$model$predictors_names)
          proj_tmp <- fit$.internals[[i]]$model$project(newdata = dplyr::select(nd, dplyr::any_of(pred_tmp)),
                                                        layer = layer)

          # make sure names match
          names(proj_tmp) <- fit$.internals[[i]]$name

          # add to next nd using cell id (previously issues with xy coords due to numerical diff)
          nd <- dplyr::left_join(x = nd, y = terra::as.data.frame(proj_tmp, cells = TRUE),
                                 by = "cell")

        }
      }

      # Project suitability
      pred_tmp <- c("x", "y", fit$model$predictors_names)
      out <- fit$project(newdata = dplyr::select(nd, dplyr::any_of(pred_tmp)), layer = layer)
      names(out) <- paste0("suitability", "_", layer, "_", as.numeric(step))
      if(is.na(terra::crs(out))) terra::crs(out) <- terra::crs( background )

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
                           "sdd_fixed" = .sdd_fixed(baseline_threshold,
                                                    new_suit = out,
                                                    value = scenario_constraints$dispersal$params[1],
                                                    resistance = resistance ),
                           "sdd_nexpkernel" = .sdd_nexpkernel(baseline_threshold,
                                                              new_suit = out,
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
        scenario_threshold <- scenario_threshold[[1]]
        out_thresh <- out
        out_thresh[out_thresh < scenario_threshold] <- 0; out_thresh[out_thresh >= scenario_threshold] <- 1
        names(out_thresh) <- "threshold"

        # Apply minimum size constraint if set
        if("min_size" %in% names(scenario_constraints)){
          if(scenario_constraints$min_size$method=="min_size"){
            out_thresh <- st_minsize(obj = out_thresh,
                                     value = scenario_constraints$min_size$params["value"] |>
                                       as.numeric(),
                                     unit = scenario_constraints$min_size$params["unit"] |>
                                       as.character(),
                                     establishment_step = scenario_constraints$min_size$params["establishment_step"] |>
                                       as.logical()
            )
          }
        }

        # Apply threshold to suitability projections if added as constraint
        if("threshold" %in% names(scenario_constraints)){
          if(scenario_constraints$threshold$method=="threshold"){
            out <- terra::mask(out, baseline_threshold, maskvalues = 0,
                               updatevalue = scenario_constraints$threshold$params["updatevalue"])
          }
        }

        # If threshold is found, check if it has values, so if projected suitability falls below the threshold
        names(out_thresh) <-  paste0('threshold_', step)
        if( is.na( terra::global(out_thresh, 'max', na.rm = TRUE) )[1] ){
          if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Scenario]','yellow','Threshold not found. Using last years threshold.')
          out_thresh <- baseline_threshold
        }
        if( terra::global(out_thresh, 'max', na.rm = TRUE)[,1] == 0){
          if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Scenario]','yellow','Thresholding removed all grid cells. Using last years threshold.')
          out_thresh <- baseline_threshold
        } else { baseline_threshold <- out_thresh }
        # Add to result stack
        suppressWarnings( proj_thresh <- c(proj_thresh, out_thresh) )
      }
      # Add to result stack
      suppressWarnings( proj <- c(proj, out) )
      if(getOption('ibis.setupmessages', default = TRUE)) pb$tick()
    }
    if(getOption('ibis.setupmessages', default = TRUE)) rm(pb)
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

          if(getOption('ibis.setupmessages', default = TRUE)) pb <- progress::progress_bar$new(total = terra::nlyr(proj))
          for(lyr in 1:terra::nlyr(proj)){
            if(getOption('ibis.setupmessages', default = TRUE)) pb$tick()
            # Normalize the projected suitability rasters to be in range 0-1000 and save
            hsMap <- predictor_transform(env = proj[[lyr]], option = "norm") * 1000
            # Write as filename in the destined folder
            suppressWarnings(
              terra::writeRaster(x = hsMap, filename = paste0( params[["hsMap"]],lyr,".tif"),
                                  datatype = "INT2S", NAflag = -9999, overwrite = TRUE)
              )
            rm(hsMap)
          }
          if(getOption('ibis.setupmessages', default = TRUE)) rm(pb)
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
          method = "ngb", func = terra::modal, cl = FALSE
        )
        scenario_constraints$boundary$params$layer <- terra::extend(scenario_constraints$boundary$params$layer,
                                                                    proj)

      }
      proj <- terra::mask(proj, scenario_constraints$boundary$params$layer)
      # Get background and ensure that all values outside are set to 0
      proj[is.na(proj)] <- 0
      proj <- terra::mask(proj, background )
      # Also for thresholds if existing
      if(terra::nlyr(proj_thresh)>0 && terra::hasValues(proj_thresh) ){
        proj_thresh <- terra::mask(proj_thresh, scenario_constraints$boundary$params$layer)
        proj_thresh[is.na(proj_thresh)] <- 0
        proj_thresh <- terra::mask(proj_thresh, background )
      }
    }

    # Should stabilization be applied?
    if(stabilize){
      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Scenario]','green','Applying stabilization.')
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
        new_proj <- terra::app(proj, fun = impute.loess)
        # Rename again
        names(new_proj) <- names(proj)
        terra::time(new_proj) <- times
        proj <- new_proj; rm(new_proj)
        # Were thresholds calculated? If yes, recalculate on the smoothed estimates
        if(terra::nlyr(proj_thresh)>0 && terra::hasValues(proj_thresh)){
          new_thresh <- proj
          new_thresh[new_thresh < scenario_threshold[1]] <- 0
          new_thresh[new_thresh >= scenario_threshold[1]] <- 1
          names(new_thresh) <- names(proj_thresh)
          thresh <- new_thresh; rm(new_thresh)
        }
      }
    }

    # Optional simulation steps
    if(!is.Waiver(scenario_simulations)){
      if("steps" %in% scenario_simulations$simulation$method ){
        pops <- .simulate_steps(proj, scenario_simulations)
        if(!is.null(pops)){
          assertthat::assert_that(is.Raster(pops))
          pops <- stars::st_as_stars(pops,
                                     crs = sf::st_crs(new_crs)
          ); names(pops) <- 'population'
        }
      }
    } else {
      pops <- NULL
    }

    # --------------------------------- #
    # # # # # # # # # # # # # # # # # # #
    # Finally convert to stars and rename
    if(terra::nlyr(proj)==1){
      # For layers with a single time step, use conversion function
      proj <- raster_to_stars(proj)
      names(proj) <- 'suitability'
    } else {
      proj <- stars::st_as_stars(proj,
                                 crs = sf::st_crs(new_crs)
      ); names(proj) <- 'suitability'
    }

    if(terra::hasValues(proj_thresh)){
      if(terra::nlyr(proj_thresh)==1){
        if(is.na(terra::time(proj_thresh))) terra::time(proj_thresh) <- terra::time(proj)
        proj_thresh <- raster_to_stars(proj_thresh)
        names(proj_thresh) <- 'threshold'
      } else {
        # Small check as there are bugs with categorical data
        if(terra::is.factor(proj_thresh[[1]])) {
          proj_thresh <- terra::as.int(proj_thresh)
          proj_thresh[proj_thresh<0] <- 0 # Negative values (sometimes occur from factor conversion) should be 0
        }
        # Add the thresholded maps as well
        proj_thresh <- stars::st_as_stars(proj_thresh,
                                          crs = sf::st_crs(new_crs)
        ); names(proj_thresh) <- 'threshold'
      }
      # Correct band if different
      if(all(!stars::st_get_dimension_values(proj, 3) != stars::st_get_dimension_values(proj_thresh, 3 ))){
        proj_thresh <- stars::st_set_dimensions(proj_thresh, 3, values = stars::st_get_dimension_values(proj, 3))
      }
      # Try
      new <- try({ c(proj, proj_thresh) },silent = TRUE)
      if(inherits(new, "try-error")){
        # Replace dimensions and check again.
        # This error likely originates from some part of the modelling routing
        # reprojecting input layers
        stars::st_dimensions(proj_thresh) <- stars::st_dimensions(proj)
        new <- try({ c(proj, proj_thresh) },silent = TRUE)
      }
      # Add the threshold as an attribute for later
      attr(new, "threshold") <- scenario_threshold
      proj <- new; rm(new)
    }
    if(!is.null(pops)){
      stars::st_dimensions(pops) <- stars::st_dimensions(proj)

      proj <- try({ c(proj, pops) },silent = TRUE)
    }
    # # # # # # # # # # # # # # # # # # #
    # --------------------------------- #

    # Return output by adding it to the scenario object
    out <- mod$clone(deep = TRUE)
    out$scenarios <- proj
    out$scenarios_migclim <- mc
    return(out)
  }
)
