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
#' @usage \S4method{mod}{BiodiversityScenario}(mod)
methods::setMethod(
  "project",
  methods::signature(mod = "BiodiversityScenario"),
  function(mod, ...){
    assertthat::assert_that(
      inherits(mod, "BiodiversityScenario"),
      !is.Waiver(mod$get_predictors())
    )

    # Get the model object
    fit <- mod$get_model()
    # Get predictors
    new_preds <- mod$get_predictors()
    if(is.Waiver(new_preds)) stop('No future predictors found.')
    new_crs <- new_preds$get_projection()
    if(is.na(new_crs)) if(getOption('ibis.setupmessages')) myLog('[Scenario]','yellow','Missing projection of future predictors.')

    # Check that predictor names are all present
    mod_pred_names <- fit$model$predictors_names
    pred_names <- mod$get_predictor_names()
    assertthat::assert_that( all(mod_pred_names %in% pred_names),
                             msg = 'Model predictors are missing from the scenario predictor!')

    # Get constrains and other parameters
    scenario_threshold <- mod$get_threshold()
    scenario_constrains <- mod$get_constrains()

    # Now convert to data.frame and subset
    df <- new_preds$get_data(df = TRUE)
    assertthat::assert_that(nrow(df)>0,
                            hasName(df,'x'), hasName(df,'y'), hasName(df,'Time'))
    df <- subset(df, select = c("x", "y", "Time", mod_pred_names) )
    # convert time dimension to Posixct
    df$Time <- as.POSIXct( df$Time )
    # Convert all units classes to numeric to avoid problems
    df <- units::drop_units(df)

    if(!is.Waiver(mod$get_scenarios())) if(getOption('ibis.setupmessages')) myLog('[Scenario]','red','Overwriting existing scenarios...')

    # ------------------ #
    if(getOption('ibis.setupmessages')) myLog('[Scenario]','green','Starting scenario predictions...')

    # Now for each unique element, loop and project in order
    proj <- raster::stack()
    pb <- progress::progress_bar$new(total = length(unique(df$Time)))
    # TODO: Do this in parallel
    for(times in sort(unique(df$Time))){
      nd <- subset(df, Time == times)
      out <- fit$project(newdata = nd)
      names(out) <- paste0('suitability_', times)
      proj <- raster::addLayer(proj, out)
      pb$tick()
    }
    rm(pb)
    proj <- raster::setZ(proj, as.Date(sort(unique(df$Time))) )
    raster::projection(proj) <- new_crs

    # ---- #
    # Calculate thresholds if set
    if(!is.Waiver(scenario_threshold)){
      # FIXME: Currently this works only for mean thresholds. Think of how the other are to be handled
      scenario_threshold <- scenario_threshold[1]
      if(getOption('ibis.setupmessages')) myLog('[Scenario]','green','Applying thresholds...')
      proj_thres <- proj
      proj_thres[proj_thres < scenario_threshold] <- 0; proj_thres[proj_thres >= scenario_threshold] <- 1
      names(proj_thres) <- gsub(pattern = 'suitability_',replacement = 'threshold_',names(proj_thres))
      proj_thres <- raster::setZ(proj_thres, as.Date(sort(unique(df$Time))) )
    }

    # ---- #
    # Apply constrains if specified
    if(!is.Waiver(scenario_constrains)){
      if(getOption('ibis.setupmessages')) myLog('[Scenario]','green','Applying constrains...')

    }

    # Finally convert to stars and rename
    proj <- stars::st_as_stars(proj,
                               crs = sf::st_crs(new_crs)
    ); names(proj) <- 'suitability'

    if((!is.Waiver(scenario_threshold))){
      # Add the thresholded maps as well
      proj_thres <- stars::st_as_stars(proj_thres,
                                       crs = sf::st_crs(new_crs)
      ); names(proj_thres) <- 'threshold'
      proj <- c(proj, proj_thres)
    }

    # Return output by adding it to the scenario object
    bdproto(NULL, mod, scenarios = proj )
  }
)
