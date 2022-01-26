#' @include utils.R bdproto-biodiversityscenario.R
NULL

#' Project a fitted model to new covariates
#'
#' @description Wrapper functions to project a [`BiodiversityScenario-class`] object to
#' new (future) covariates
#' @param mod A [`BiodiversityScenario`] object with set predictors
#' @param no_projection A [`logical`] flag whether future projection should be created. (Default: TRUE)
#' Note that some constrains such as [MigClim] can still simulate future change without projections.
#' @param ... passed on parameters
#' @returns Saves [`stars`] objects of the obtained predictions in mod.
#'
#' @name project
#' @aliases project
#' @keywords scenario
#' @exportMethod project
#' @export
NULL
methods::setGeneric("project",
                    signature = methods::signature("mod"),
                    function(mod, no_projection = TRUE,...) standardGeneric("project"))

#' @name project
#' @rdname project
#' @usage \S4method{project}{BiodiversityScenario, logical}(mod, no_projection)
methods::setMethod(
  "project",
  methods::signature(mod = "BiodiversityScenario"),
  function(mod, no_projection = TRUE, ...){
    assertthat::assert_that(
      inherits(mod, "BiodiversityScenario"),
      !is.Waiver(mod$get_predictors()),
      is.logical(no_projection)
    )
    if(!is.Waiver(mod$get_scenarios())) if(getOption('ibis.setupmessages')) myLog('[Scenario]','red','Overwriting existing scenarios...')

    # Get the model object
    fit <- mod$get_model()
    assertthat::assert_that(!is.Waiver(fit), msg = "No model found!")
    # Get predictors
    new_preds <- mod$get_predictors()
    if(is.Waiver(new_preds)) stop('No future predictors found.')
    new_crs <- new_preds$get_projection()
    if(is.na(new_crs)) if(getOption('ibis.setupmessages')) myLog('[Scenario]','yellow','Missing projection of future predictors.')

    # Get limits if present
    if(!is.null( mod$get_limits() )){
      # Get prediction
      n <- fit$show_rasters()[grep("threshold",fit$show_rasters())]
      tr <- fit$get_data(n)[[1]]
      tr <- cbind( raster::coordinates(tr), data.frame(thresh = values(tr)))
      tr[['thresh']] <- ifelse(tr[['thresh']]==0, NA, tr[['thresh']])
      tr <- tr %>% subset(., complete.cases(thresh))

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
                             msg = 'Model predictors are missing from the scenario predictor!')

    # Get constrains and other parameters
    scenario_threshold <- mod$get_threshold()
    # Not get the baseline raster
    thresh_reference <- grep('threshold',fit$show_rasters(),value = T)[1] # Use the first one (mean)
    baseline_threshold <- mod$get_model()$get_data(thresh_reference)
    if(inherits(baseline_threshold, 'RasterStack') || inherits(baseline_threshold, 'RasterBrick')){
      baseline_threshold <- baseline_threshold[[grep("mean",names(baseline_threshold))]] # FIXME: Potentially have an option for this
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
          if(getOption('ibis.setupmessages')) myLog('[Scenario]','yellow','Using kissmig to calculate updated distribution thresholds.')
          scenario_threshold <- new_waiver()
        }
      } else {
        assertthat::assert_that(!is.Waiver(scenario_threshold),msg = "Other constrains require threshold option!")
      }
    }
    if(is.Waiver(baseline_threshold) && !is.Waiver(scenario_constraints)) stop("No baseline threshold layer found!")
    if("connectivity" %in% names(scenario_constraints) && "dispersal" %notin% names(scenario_constraints)){
      if(getOption('ibis.setupmessages')) myLog('[Scenario]','red','Connectivity contraints need a set dispersal constraint.')
    }
    # ----------------------------- #
    # Start of projection           #
    # ----------------------------- #

    # Now convert to data.frame and subset
    df <- new_preds$get_data(df = TRUE)
    names(df)[1:3] <- tolower(names(df)[1:3])
    assertthat::assert_that(nrow(df)>0,
                            hasName(df,'x'), hasName(df,'y'), hasName(df,'time'))
    df <- subset(df, select = c("x", "y", "time", mod_pred_names) )
    # convert time dimension to Posixct
    df$time <- as.POSIXct( df$time )
    # Convert all units classes to numeric to avoid problems
    df <- units::drop_units(df)

    # ------------------ #
    if(getOption('ibis.setupmessages')) myLog('[Scenario]','green','Starting suitability predictions...')

    # Now for each unique element, loop and project in order
    proj <- raster::stack()
    proj_thresh <- raster::stack()

    #TODO: Add no_projection step
    pb <- progress::progress_bar$new(total = length(unique(df$time)))
    # TODO: Do this in parallel
    for(times in sort(unique(df$time))){
      nd <- subset(df, time == times)
      # Project suitability
      # FIXME: Adapt for uncertainty projections as well!
      out <- fit$project(newdata = nd)[[1]] # First prediction being the mean
      names(out) <- paste0("suitability", "_", times)
      if(is.na(raster::projection(out))) raster::projection(out) <- raster::projection( fit$model$background )

      # If constrains are set, apply them
      if(!is.Waiver(scenario_constraints)){
        # Calculate dispersal constraint if set
        if("dispersal" %in% names(scenario_constraints) ){
          # MigClim simulations are run posthoc
          if(scenario_constraints$dispersal$method %in% c("sdd_fixed", "sdd_nexpkernel")){
            out <- switch (scenario_constraints$dispersal$method,
                           "sdd_fixed" = .sdd_fixed(baseline_threshold, out,
                                                    value = scenario_constraints$dispersal$params[1],
                                                    resistance = scenario_constraints$connectivity$params$resistance ),
                           "sdd_nexpkernel" = .sdd_nexpkernel(baseline_threshold, out,
                                                              value = scenario_constraints$dispersal$params[1],
                                                              resistance = scenario_constraints$connectivity$params$resistance)
            )
            names(out) <-  paste0('suitability_', times)
          }
        }
        # For kissmig generate threshold and masked suitabilities
        if(scenario_constraints$dispersal$method %in% c("kissmig")){
          out <- .kissmig_dispersal(baseline_threshold,
                                    new_suit = out,
                                    resistance = scenario_constraints$connectivity$params$resistance,
                                    params = scenario_constraints$dispersal$params)
          # Returns a layer of two with both the simulated threshold and the masked suitability raster
          names(out) <- paste0(c('threshold_', 'suitability_'), times)
          # Add threshold to result stack
          proj_thresh <- raster::addLayer(proj_thresh, out[[1]] )
          baseline_threshold <- out[[1]]
          out <- out[[2]]
        }

        # # Connectivity constraints
        if("connectivity" %in% names(scenario_constraints)){
          # By definition a hard barrier removes all suitable
          if(scenario_constraints$connectivity$method == "hardbarrier"){
            out[scenario_constraints$connectivity$params$resistance] <- 0
          }
        }

      }
      # Calculate thresholds if set manually
      if(!is.Waiver(scenario_threshold)){
        # FIXME: Currently this works only for mean thresholds. Think of how the other are to be handled
        scenario_threshold <- scenario_threshold[1]
        out_thresh <- out
        out_thresh[out_thresh < scenario_threshold] <- 0; out_thresh[out_thresh >= scenario_threshold] <- 1
        names(out_thresh) <-  paste0('threshold_', times)
        # If threshold is
        if( cellStats(out_thresh, 'max') == 0){
          if(getOption('ibis.setupmessages')) myLog('[Scenario]','yellow','Thresholding removed all grid cells. Using last years threshold.')
          out_thresh <- baseline_threshold
        } else { baseline_threshold <- out_thresh }
        # Add to result stack
        proj_thresh <- raster::addLayer(proj_thresh, out_thresh)
      }
      # Add to result stack
      proj <- raster::addLayer(proj, out)
      pb$tick()
    }
    rm(pb)
    proj <- raster::setZ(proj, as.Date(sort(unique(df$time))) )
    if(raster::nlayers(proj_thresh)>1) proj_thresh <- raster::setZ(proj_thresh, as.Date(sort(unique(df$time))) )

    # Apply MigClim if set
    # FIXME: Ideally make this whole setup more modular. So create suitability projections first
    if(!is.Waiver(scenario_constraints)){
      # Calculate dispersal constraint if set
      if("dispersal" %in% names(scenario_constraints) ){
        # MigClim simulations are run posthoc
        if(scenario_constraints$dispersal$method == 'MigClim'){
          # Get Parameters
          params <- scenario_constraints$dispersal$params

          pb <- progress::progress_bar$new(total = raster::nlayers(proj))
          for(lyr in 1:raster::nlayers(proj)){
            pb$tick()
            # Normalize the projected suitability rasters to be in range 0-1000 and save
            hsMap <- predictor_transform(env = proj[[lyr]], option = "norm") * 1000
            # Write as filename in the destined folder
            suppressWarnings(
              raster::writeRaster(x = hsMap, filename = paste0( params[["hsMap"]],lyr,".tif"),
                                  dt = "INT2S", varNA = -9999, prj = TRUE, overwrite = TRUE)
              )
            rm(hsMap)
          };rm(pb)
          # For threshold, define based on type
          tr <- ifelse(params[["rcThreshold"]] == "continuous",
                       0,
                       750) # Set to 75% as suitability threshold

          # Now run MigClim by switching to temporary dir
          dir.ori <- getwd()
          setwd(params[["dtmp"]])
          try({
            m <- MigClim::MigClim.migrate(
              iniDist = basename(params[["iniDist"]]),
              hsMap = basename(params[["hsMap"]]),
              rcThreshold = tr,
              envChgSteps = raster::nlayers(proj), # Use number of projected suitability layers
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
          run_stats <- read.table(
            file.path(basename(params[["simulName"]]), paste0(basename(params[["simulName"]]),"_stats.txt")),
                                  header = TRUE
            )
          run_sums <- read.table(
            file.path(basename(params[["simulName"]]), paste0(basename(params[["simulName"]]),"_summary.txt")),
            header = TRUE
          )
          # Get MigClim outputs
          ll <- list.files(basename(params[["simulName"]]),'asc',full.names = TRUE)
          run_sims <- raster::stack(ll); names(run_sims) <- tools::file_path_sans_ext(basename(ll))
          # Condense the simulation runs into one modal prediction
          run_sim <- raster::calc(run_sims, raster::modal)
          raster::projection(run_sim) <- raster::projection(fit$get_data('prediction'))
          all(sapply(list.files(getwd(),".tif"), file.remove)) # Cleanup
          setwd(dir.ori) # Flip back to original directory

          # Wrap all results in a list
          mc <- list(params = params,
                     stats = run_stats,
                     summary = run_sums,
                     raster = run_sim)
        }
      }
    } # End of MigClim processing chain
    # If not found, set a waiver
    if(!exists("mc")) mc <- new_waiver()
    # ---- #
    assertthat::assert_that(
      is.Raster(proj), is.Raster(proj_thresh),
      msg = "Something went wrong with the projection."
    )

    # Finally convert to stars and rename
    proj <- stars::st_as_stars(proj,
                               crs = sf::st_crs(new_crs)
    ); names(proj) <- 'suitability'

    if(raster::nlayers(proj_thresh)>0){
      # Add the thresholded maps as well
      proj_thresh <- stars::st_as_stars(proj_thresh,
                                       crs = sf::st_crs(new_crs)
      ); names(proj_thresh) <- 'threshold'
      proj <- stars:::c.stars(proj, proj_thresh)
    }

    # Return output by adding it to the scenario object
    bdproto(NULL, mod,
            scenarios = proj,
            scenarios_migclim = mc
            )
  }
)
