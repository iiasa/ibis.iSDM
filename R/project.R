#' @include utils.R bdproto-biodiversityscenario.R
NULL

#' Project a fitted model to new covariates
#'
#' @description Wrapper functions to project a [`BiodiversityScenario-class`] object to
#' new (future) covariates
#' @param mod A [`BiodiversityScenario`] object with set predictors
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
                    function(mod,...) standardGeneric("project"))

#' @name project
#' @rdname project
#' @usage \S4method{project}{BiodiversityScenario}(mod)
methods::setMethod(
  "project",
  methods::signature(mod = "BiodiversityScenario"),
  function(mod, ...){
    assertthat::assert_that(
      inherits(mod, "BiodiversityScenario"),
      !is.Waiver(mod$get_predictors())
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
    assertthat::assert_that(scenario_threshold > 0,msg = "Threshold has to be larger than 0")
    # Not get the baseline raster
    thresh_reference <- grep('threshold',fit$show_rasters(),value = T)[1] # Use the first one (mean)
    baseline_threshold <- mod$get_model()$get_data(thresh_reference)
    if(inherits(baseline_threshold, 'RasterStack') || inherits(baseline_threshold, 'RasterBrick')){
      baseline_threshold <- baseline_threshold[[grep("mean",names(baseline_threshold))]] # FIXME: Potentially have an option for this
    }
    scenario_constraints <- mod$get_constraints()

    # Check that thresholds are set for constrains
    if(is.Waiver(scenario_threshold) && !is.Waiver(scenario_constraints)) stop("Constrains require calculated thresholds!")
    if(is.Waiver(baseline_threshold) && !is.Waiver(scenario_constraints)) stop("No baseline threshold layer found!")
    if("connectivity" %in% names(scenario_constraints) && "dispersal" %notin% names(scenario_constraints)){
      if(getOption('ibis.setupmessages')) myLog('[Scenario]','red','Connectivity contraints need a set dispersal constraint.')
    }

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
    if(getOption('ibis.setupmessages')) myLog('[Scenario]','green','Starting scenario predictions...')

    # Now for each unique element, loop and project in order
    proj <- raster::stack()
    proj_thresh <- raster::stack()

    pb <- progress::progress_bar$new(total = length(unique(df$time)))
    # TODO: Do this in parallel
    for(times in sort(unique(df$time))){
      nd <- subset(df, time == times)
      # Project suitability
      # FIXME: Adapt for uncertainty projections as well!
      out <- fit$project(newdata = nd)[[1]] # First prediction being the mean
      names(out) <- paste0("suitability", "_", times)

      # If constrains are set, apply them
      if(!is.Waiver(scenario_constraints)){
        # Calculate dispersal constraint if set
        if("dispersal" %in% names(scenario_constraints) ){
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
        # # Connectivity constraints
        if("connectivity" %in% names(scenario_constraints)){
          # By definition a hard barrier removes all suitable
          if(scenario_constraints$connectivity$method == "hardbarrier"){
            out[scenario_constraints$connectivity$params$resistance] <- 0
          }
        }

      }
      # Calculate thresholds if set
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
    # ---- #
    assertthat::assert_that(
      is.Raster(proj), is.Raster(proj_thresh),
      compareRaster(proj, proj_thresh),msg = "Something went wrong with the projection."
    )

    # Finally convert to stars and rename
    proj <- stars::st_as_stars(proj,
                               crs = sf::st_crs(new_crs)
    ); names(proj) <- 'suitability'

    if((!is.Waiver(scenario_threshold))){
      # Add the thresholded maps as well
      proj_thresh <- stars::st_as_stars(proj_thresh,
                                       crs = sf::st_crs(new_crs)
      ); names(proj_thresh) <- 'threshold'
      proj <- stars:::c.stars(proj, proj_thresh)
    }

    # Return output by adding it to the scenario object
    bdproto(NULL, mod, scenarios = proj)
  }
)
