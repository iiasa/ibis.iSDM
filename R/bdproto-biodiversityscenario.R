#' @include utils.R bdproto.R
NULL

#' @export
if (!methods::isClass("BiodiversityScenario")) methods::setOldClass("BiodiversityScenario")
NULL

#' Prototype for a biodiversity scenario from a trained model
#'
#' Base [`proto`] class for any biodiversity scenario objects.
#' Serves as container that supplies data and functions to
#' other [`proto`] classes.
#'
#' @name BiodiversityScenario-class
#' @aliases BiodiversityScenario
#' @family bdproto
#' @keywords bdproto
NULL

#' @export
BiodiversityScenario <- bdproto(
  "BiodiversityScenario",
  modelobject = new_waiver(), # The id of the model
  modelid = new_waiver(),
  limits = new_waiver(),
  predictors = new_waiver(),
  constraints = new_waiver(),
  scenarios = new_waiver(),
  # Print message with summary of model
  print = function(self) {
    # Check that model exists
    fit <- self$get_model()
    timeperiod <- self$get_timeperiod()
    # Get set predictors and time period
    pn = ifelse(is.Waiver(self$get_predictor_names()),'None',name_atomic(self$get_predictor_names(), "predictors"))
    tp = ifelse(is.Waiver(timeperiod),'None',
                paste0(
                  paste0( timeperiod,collapse = ' -- '),
                  ' (',round(as.numeric(difftime(self$get_timeperiod()[2],self$get_timeperiod()[1],unit="weeks"))/52.25,1),' years)'
                )
    )
    # Constrains
    cs <- self$get_constraints()
    if(!is.Waiver(cs)) cs <- vapply(cs, function(x) x$method, character(1))
    # Thresholds
    tr <- self$get_threshold()

    message(paste0(
      ifelse(is.Waiver(self$limits),"Spatial-temporal scenario:","Spatial-temporal scenario (limited):"),
                   '\n  Used model: ',ifelse(is.Waiver(fit) || isFALSE(fit), text_red('None'), class(fit)[1] ),
                   "\n --------- ",
                   "\n  Predictors:     ", pn,
                   "\n  Time period:    ", tp,
                   ifelse(!is.Waiver(cs)||!is.Waiver(tr), "\n --------- ", ""),
                   ifelse(is.Waiver(cs),"", paste0("\n  Constraints:    ", text_green(paste(paste0(names(cs),' (',cs,')'),collapse = ', ')) ) ),
                   ifelse(is.Waiver(tr),"", paste0("\n  Threshold:      ", round(tr[1], 3),' (',names(tr[1]),')') ),
                   "\n --------- ",
                   "\n  Scenarios fitted: ", ifelse(is.Waiver(self$scenarios),text_yellow('None'), text_green('Yes'))
      )
    )
  },
  # Verify that set Model exist and check self-validity
  verify = function(self){
    assertthat::validate_that( !is.Waiver(self$modelobject),
                               !is.Waiver(self$modelid)
                               )
    # Get Model object and check that ID is correct
    if(inherits(self$modelobject, "DistributionModel")){
      x <- self$modelobject
    } else {
      x <- get(self$modelobject)
    }
    assertthat::validate_that(x$id == self$modelid)
    # Check that objects are correctly set or found
    assertthat::assert_that(is.Waiver(self$get_predictors()) || inherits(self$get_predictors(), "PredictorDataset"))
    assertthat::assert_that(is.Waiver(self$get_data()) || (inherits(self$get_data(), "stars") || is.Raster(self$get_data())) )
    assertthat::assert_that(is.Waiver(self$get_constraints()) || is.list(self$get_constraints()))
    invisible()
  },
  # Show the name of the Model
  show = function(self) {
    if(is.character(self$modelobject)){
      return( self$modelobject )
    } else {
      return( fit$model$runname )
    }
  },
  # Get projection
  get_projection = function(self){
    return( self$predictors$get_projection() )
  },
  # Get resolution
  get_resolution = function(self){
    return( self$predictors$get_resolution() )
  },
  # Get Model
  get_model = function(self){
    if(is.Waiver(self$modelobject)) return( new_waiver() )
      else {
        if(inherits(self$modelobject, "DistributionModel")){
          return( self$modelobject )
        } else {
          if(!exists(self$modelobject)) return( FALSE ) else {
            return( get(self$modelobject) )
            }
        }
      }
  },
  # Get provided limits
  get_limits = function(self){
    if(is.Waiver(self$limits)) return(NULL)
    return(self$limits)
  },
  # Get Model predictors
  get_predictor_names = function(self) {
    if(is.Waiver(self$predictors)) return(self$predictors)
    if(inherits(self$predictors, "PredictorDataset")) {
      self$predictors$get_names()
    } else {
      stop("Feature data is of an unrecognized class")
    }
  },
  # Get time period of projection
  get_timeperiod = function(self, what = "range"){
    if(is.Waiver(self$predictors)) return(self$predictors)
    if(inherits(self$predictors, "PredictorDataset")) {
      if(what == "range"){
        return(
          c( min(as.Date(self$predictors$timeperiod)), max(as.Date(self$predictors$timeperiod)) )
        )
      } else {
        return(
          sort( self$predictors$timeperiod )
        )
      }
    }
  },
  # Get constrains for model
  get_constraints = function(self){
    return( self$constraints )
  },
  # Get thresholds if specified
  get_threshold = function(self){
    if('threshold' %notin% names(self)) return( new_waiver() )
    return( self$threshold )
  },
  # Apply specific threshold
  apply_threshold = function(self, tr = new_waiver()){
    # Assertions
    if(is.Waiver(tr)) assertthat::assert_that( is.numeric(self$threshold), msg = 'No threshold value found.')
    assertthat::assert_that( !is.Waiver(self$scenarios), msg = 'No scenarios found.')
    # Get prediction and threshold
    sc <- self$get_data()
    if(!is.Waiver(tr)) tr <- self$threshold
    # reclassify to binary
    sc[sc < tr] <- 0; sc[sc >= tr] <- 1
    names(sc) <- 'presence'
    return(sc)
  },
  # Set Predictors
  set_predictors = function(self, x){
    assertthat::assert_that(inherits(x, "PredictorDataset"))
    bdproto(NULL, self, predictors = x)
  },
  # Set constrains
  set_constraints = function(self, x){
    if(!is.Waiver(self$get_constraints())){
      cr <- self$get_constraints()
      # FIXME: Remove duplicates
      bdproto(NULL, self, constraints = c(cr, x))
    } else {
      bdproto(NULL, self, constraints = x)
    }
  },
  # Get Predictors
  get_predictors = function(self){
    return(self$predictors)
  },
  # Remove predictors
  rm_predictors = function(self, names){
    if(is.Waiver(self$predictors) || is.null(self$predictors)) return(NULL)
    if(missing(names)){
      names <- self$get_predictor_names() # Assume all names
    }
    assertthat::assert_that(
      is.character(names) || assertthat::is.scalar(names) || is.vector(names)
    )
    # Get predictor collection
    prcol <- bdproto(NULL, self)
    # Set the object
    prcol$predictors$rm_data(names)
    if(length(prcol$get_predictor_names())==0) prcol$predictors <- new_waiver()
    return(prcol)
  },
  # Get scenario predictions
  get_data = function(self, what = "scenarios"){
    return(self[[what]])
  },
  # Plot the prediction
  plot = function(self, what = "suitability", which = NULL, ...){
    if(is.Waiver(self$get_data())){
      if(getOption('ibis.setupmessages')) myLog('[Scenario]','red','No scenarios found')
      invisible()
    } else {
      # Get unique number of data values. Surely there must be an easier val
      vals <- self$get_data()[what] %>% stars:::pull.stars() %>% as.vector() %>% unique()
      vals <- length(na.omit(vals))
      if(vals>2) col <- ibis_colours$sdm_colour else col <- c('grey25','coral')
      if(is.null(which)){
        stars:::plot.stars( self$get_data()[what], breaks = "equal", col = col )
      } else {
        # Assert that which is actually within the dimensions
        assertthat::assert_that(which <= dim(self$get_data())[3],
                                msg = "Band selection out of bounds.")
        obj <- self$get_data()[what,,,which]
        stars:::plot.stars( obj, breaks = "equal", col = col,
                            main = paste0(what," for ", stars::st_get_dimension_values(obj, "band") )  )
      }
    }
  },
  # Convenience function to plot thresholds if set
  plot_threshold = function(self, which = NULL){
    # Check that baseline and scenario thresholds are all there
    if(!( 'threshold' %in% attributes(self$get_data())$names )) return(new_waiver())

    self$plot(what = "threshold", which = which)
  },
  # Plot Migclim results if existing
  plot_migclim = function(self){
    # Get scenarios
    mc <- self$get_data("scenarios_migclim")
    if(is.Waiver(mc)) return(mc)

    # Otherwise plot the raster
    ras <- mc$raster

    # Colour coding from MigClim::MigClim.plot
    rstVals <- sort(raster::unique(ras))
    negativeNb <- length(which(rstVals < 0))
    positiveNb <- length(which(rstVals > 1 & rstVals < 30000))
    zeroExists <- any(rstVals == 0)
    oneExists <- any(rstVals == 1)
    unilimtedExists <- any(rstVals == 30000)
    Colors <- rep("yellow", negativeNb)
    if(zeroExists) Colors <- c(Colors, "grey94")
    if (oneExists) Colors <- c(Colors, "black")
    Colors <- c(Colors, rainbow(positiveNb, start = 0, end = 0.4))
    if (unilimtedExists) Colors <- c(Colors, "pink")

    # Plot
    # 0 - Cells that have never been occupied and are unsuitable habitat at the end of the simulation
    # 1 - Cells that belong to the species' initial distribution and that have remained occupied during the entire simulation.
    # 1 < value < 30 000 - determine the dispersal step during which it was colonized. E.g. 101 is first dispersal even in first step
    # 30 0000 - Potentially suitable cells that remained uncolonized
    # <0 - Negative values indicate cells that were once occupied but have become decolonized. Code as for colonization
    dev.new(width = 7, height = 7 * ((ymax(ras) - ymin(ras))/(xmax(ras) - xmin(ras))))
    plot(ras, col = Colors, breaks = c(min(rstVals) - 1, rstVals), legend = FALSE,
         main = "Newly colonized and stable habitats")
  },
  # Plot animation of scenarios
  plot_animation = function(self, what = "suitability", fname = NULL){
    assertthat::assert_that(!is.Waiver(self$get_data()) )
    check_package('gganimate')
    # Get scenarios
    obj <- self$get_data()[what]

    # Make the animation plot
    g <- ggplot2::ggplot() +
      stars::geom_stars(data = obj, downsample = c(1,1,0)) +
      ggplot2::coord_equal() +
      ggplot2::theme_bw(base_size = 20) +
      ggplot2::scale_x_discrete(expand=c(0,0)) +
      ggplot2::scale_y_discrete(expand=c(0,0)) +
      ggplot2::scale_fill_gradientn(colours = ibis_colours$sdm_colour, na.value = NA) +
      # Animation options
      gganimate::transition_time(band) +
      gganimate::ease_aes() +
      ggplot2::labs(x = '', y ='', title = "{frame_time}")

    if(is.null(fname)){
      gganimate::anim_save(animation = g, fname)
    } else { g }
  },
  #Plot relative change between baseline and projected thresholds
  plot_relative_change = function(self, position = NULL, variable = 'mean', plot = TRUE){
    # Default position is the last one
    assertthat::assert_that(is.null(position) || is.numeric(position) || is.character(position),
                            is.character(variable))
    # Threshold
    obj <- self$get_model()
    thresh_reference <- grep('threshold',obj$show_rasters(),value = T)
    # If there is more than one threshold only use the one from variable
    if(length(thresh_reference)>1) {
      warning('More than one baseline threshold. Using the first one.')
      thresh_reference <- grep(variable, thresh_reference,value = T)[1]
    }
    # Check that baseline and scenarios are all there
    assertthat::assert_that(
      !is.Waiver(self$get_data()),
      'threshold' %in% attributes(self$get_data())$names,
      length(thresh_reference) >0 & is.character(thresh_reference),
      is.Raster( self$get_model()$get_data('prediction') ),
      msg = "Threshold not found!"
    )

    # Not get the baseline raster
    baseline <- self$get_model()$get_data(thresh_reference)
    # And the last scenario prediction
    scenario <- self$get_data()['threshold']
    time <- stars::st_get_dimension_values(scenario, which = 3) # 3 assumed to be time band
    if(is.numeric(position)) position <- time[position]
    if(is.null(position)) position <- time[length(time)]
    final <- scenario %>%
      stars:::filter.stars(band == position) %>%
      as('Raster')
    raster::projection(final) <- raster::projection(baseline)
    # -- #
    if(!inherits(final, 'RasterLayer')) final <- final[[1]] # In case it is a rasterbrick or similar
    if(!compareRaster(baseline, final,stopiffalse = FALSE)) final <- alignRasters(final, baseline, cl = FALSE) # In case they somehow differ?

    # Calculate overlays
    diff_f <- raster::overlay(baseline, final, fun = function(x, y){x + y * 2})
    diff_f <- raster::ratify(diff_f)
    # 0 = Unsuitable | 1 = Loss | 2 = Gain | 3 = stable
    rat <- levels(diff_f)[[1]]
    rat <- merge.data.frame(rat, data.frame(ID = seq(0,3), diff = c("Unsuitable", "Loss", "Gain", "Stable")),all = TRUE)
    levels(diff_f) <- rat
    diff_f <- raster::mask(diff_f, baseline)
    rm(baseline, final)

    # Plot
    if(plot){
      # Colours
      cols <- c("Unsuitable" = "gray92", "Loss" = "#DE646A", "Gain" = "cyan3", "Stable" = "gray60")

      # Convert to raster
      diff_ff <- as.data.frame(diff_f, xy = TRUE)
      names(diff_ff)[3] <- "Change"
      diff_ff$Change <- factor(diff_ff$Change, levels = names(cols))

      # Use ggplot
      g <- ggplot2::ggplot() +
        ggplot2::coord_equal() +
        ggplot2::geom_raster(data = diff_ff, ggplot2::aes(x = x, y = y, fill = Change)) +
        ggplot2::theme_light(base_size = 18) +
        ggplot2::scale_x_discrete(expand=c(0,0)) +
        ggplot2::scale_y_discrete(expand=c(0,0)) +
        ggplot2::scale_fill_manual(values = cols, na.value = 'transparent') +
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::labs(x = "", y = "", title = paste0('Change between baseline and ', position))
      return(g)
    } else {
      # Return
      return(diff_f)
    }
  },
  # Summarize the change in layers between timesteps
  summary = function(self, layer = "threshold", plot = FALSE, relative = FALSE){
    # Check that baseline and scenario thresholds are all there
    assertthat::assert_that(
      !is.Waiver(self$get_data()),
      is.logical(plot), is.logical(relative)
    )
    if( layer == "threshold" & 'threshold' %in% attributes(self$get_data())$names ){
      # TODO: Try and get rid of dplyr dependency. Currently too much work to not use it
      check_package("dplyr")
      # Get the scenario predictions and from there the thresholds
      scenario <- self$get_data()['threshold']
      time <- stars::st_get_dimension_values(scenario,which = 'band')
      assertthat::assert_that(!is.na(sf::st_crs(scenario)), msg = "Scenario not correctly projected.")
      # HACK: Add area to stars
      ar <- stars:::st_area.stars(scenario)
      # Get the unit
      ar_unit <- units::deparse_unit(ar$area)
      new <- as(scenario,"Raster") * as(ar, "Raster")
      new <- raster::setZ(new, time)
      # Convert to scenarios to data.frame
      df <- stars:::as.data.frame.stars(stars:::st_as_stars(new)) %>% subset(., complete.cases(.))
      names(df)[4] <- "area"
      # --- #
      # Now calculate from this data.frame several metrics related to the area and change in area
      df <- df %>% dplyr::group_by(x,y) %>% dplyr::mutate(id = dplyr::cur_group_id()) %>%
        dplyr::ungroup() %>% dplyr::select(-x,-y) %>%
        dplyr::mutate(area = dplyr::if_else(is.na(area), 0, area)) %>% # Convert missing data to 0
        dplyr::arrange(id, band)
      df$area <- units::as_units(df$area, units::as_units(ar_unit))  # Set Units
      # Convert to km2 and remove units as this causes issues with dplyr
      df$area <- units::set_units(df$area, "km2") %>% units::drop_units()

      # Total amount of area occupied for a given time step
      out <- df %>% dplyr::group_by(band) %>% dplyr::summarise(area_km2 = sum(area, na.rm = TRUE))
      out$totarea <- raster::cellStats((new[[1]]>=0) * as(ar, "Raster"), "sum")
      if(units::deparse_unit(units::as_units(ar_unit)) == "m2") {
        out$totarea <- out$totarea / 1e6
        out <- dplyr::rename(out, totarea_km2 = totarea)
      }

      # Total amount of area lost / gained / stable since previous time step
      totchange_occ <- df %>%
          dplyr::group_by(id) %>%
          dplyr::mutate(change = (area - dplyr::lag(area)) ) %>% dplyr::ungroup() %>%
          subset(., complete.cases(.))
      o <- totchange_occ %>% dplyr::group_by(band) %>%
        dplyr::summarise(totchange_stable_km2 = sum(area[change == 0]),
                         totchange_gain_km2 = sum(change[change > 0]),
                         totchange_loss_km2 = sum(change[change < 0]))
      out <- out %>% dplyr::left_join(o, by = "band")

      if(relative == TRUE){
        # Finally calculate relative change to baseline (first entry) for all entries where this is possible
        relChange <- function(v, fac = 100) (((v- v[1]) / v[1]) * fac)
        out <- subset(out, select = c("band", "area_km2", "totarea_km2"))
        out[,c("area_km2")] <- apply( out[,c("area_km2")], 2, relChange)
      }
    } else {
      # Get the scenario predictions and from there the thresholds
      scenario <- self$get_data()['suitability']
      times <- stars::st_get_dimension_values(scenario, which = 'band')
      # Get area
      ar <- stars:::st_area.stars(scenario)
      # Get the unit
      ar_unit <- units::deparse_unit(ar$area)

      # TODO: Check whether one could not simply multiply with area (poisson > density, binomial > suitable area)
      mod <- self$get_model()
      scenario <- scenario * ar
      out <- summarise_projection(scenario, fun = "mean", relative = relative)
    }

    if(plot){
      if( 'threshold' %in% attributes(self$get_data())$names ){
        if(has_name(out,"band")) out <- dplyr::rename(out, "time" = "band")
        ggplot2::ggplot(out,
                        ggplot2::aes(x = time, y = as.numeric(area_km2))) +
          ggplot2::theme_classic(base_size = 18) +
          ggplot2::geom_line(size = 2) +
          ggplot2::labs(x = "Time", y = expression(Area(km^2)))
      } else {
        ggplot2::ggplot(out,
                        ggplot2::aes(x = time,
                                     y = suitability_q50,
                                     ymin = suitability_q25,
                                     ymax = suitability_q75)) +
          ggplot2::theme_classic(base_size = 18) +
          ggplot2::geom_ribbon(fill = "grey90") +
          ggplot2::geom_line(size = 2) +
          ggplot2::labs(x = "Time", y = expression(Area(km^2)), title = "Relative suitable habitat")
      }
    }
    return(out)
  },
  # Summarize beforeafter
  summary_beforeafter = function(self){
    # Check that baseline and scenario thresholds are all there
    assertthat::assert_that(
      !is.Waiver(self$get_data()),
      ( 'threshold' %in% attributes(self$get_data())$names ),
      msg = "This function only works with added thresholds."
    )
    scenario <- self$get_data()['threshold']

    # Get runname
    runname <- self$get_model()[["model"]]$runname

    return(
      tibble::add_column( summarise_change(scenario), runname = runname, .before = 1)
    )
  },
  # Calculate slopes
  calc_scenarios_slope = function(self, what = 'suitability', plot = TRUE, oftype = "stars"){
    if(is.Waiver(self$get_data())) return( new_waiver() )
    assertthat::assert_that(what %in% attributes(self$get_data())$names )

    oo <- self$get_data()[what]
    tt <- as.numeric( stars::st_get_dimension_values(self$scenarios, 3) )
    # Calc pixel-wise linear slope
    out <- stars::st_apply(
      oo,
      1:2,
      function(x) {
        if (anyNA(x))
          NA_real_
        else
          lm.fit(cbind(1, tt), x)$coefficients[2]
      }
    )
    names(out) <- 'linear_coefficient'
    if(oftype == "stars"){
      if(plot) stars:::plot.stars(out, breaks = "fisher", col = c(ibis_colours$divg_bluered[1:10],"grey90",ibis_colours$divg_bluered[11:20]))
    } else {
      out <- as(out, "Raster")
      if(plot) plot(out, col = c(ibis_colours$divg_bluered[1:10],"grey90",ibis_colours$divg_bluered[11:20]))
    }
    return(out)
  },
  # Save object
  save = function(self, fname, type = 'tif', dt = 'FLT4S'){
    assertthat::assert_that(
      !missing(fname),
      is.character(fname),
      is.character(type),
      !is.Waiver(self$get_data()),
      is.character(dt)
    )
    # Match input types
    type <- match.arg(type, c('gtif','gtiff','tif','nc','ncdf', 'feather'), several.ok = FALSE)
    dt <- match.arg(dt, c('LOG1S','INT1S','INT1U','INT2S','INT2U','INT4S','INT4U','FLT4S','FLT8S'), several.ok = FALSE )

    if(file.exists(fname)) warning('Overwritting existing file...')

    # Respecify type if output filename has already been set
    if(gsub('\\.','',raster::extension(fname)) != type) type <- gsub('\\.','',raster::extension(fname))

    # Change output type for stars
    dtstars <- switch(dt,
                    "LOG1S" = "Byte",
                    "INT1U" = "UInt16",
                    "INT1S" = "Int16",
                    "INT2U" = "UInt16",
                    "INT2S" = "Int16",
                    "INT4S" = "Int32",
                    "INT4U" = "UInt32",
                    "FLT4S" = "Float32",
                    "FLT8S" = "Float64"
    )

    # Get scenario object
    ras <- self$get_data()
    # If Migclim has been computed, save as well
    if(!is.Waiver(self$scenarios_migclim)) ras_migclim <- self$get_data("scenarios_migclim")

    if(type %in% c('gtif','gtiff','tif')){
      # Write stars output for every band
      for(i in 1:length(ras)){
        # Band specific output
        fname2 <- paste0( tools::file_path_sans_ext(fname), "__", names(ras)[i], raster::extension(fname))
        stars::write_stars(
          obj = ras,
          layer = i,
          dsn = fname2,
          options = c("COMPRESS=DEFLATE"),
          type = ifelse(is.factor(ras[[1]]), "Byte", dtstars),
          NA_value = NA_real_,
          update = ifelse(file.exists(fname2), TRUE, FALSE),
          normalize_path = TRUE,
          progress = TRUE
        )
      }
      if(!is.Waiver(self$scenarios_migclim)){
        fname2 <- paste0( tools::file_path_sans_ext(fname), "__migclim", raster::extension(fname))
        writeGeoTiff(ras_migclim, fname = fname, dt = dt)
      }
    } else if(type %in% c('nc','ncdf')) {
      # Save as netcdf, for now in individual files
      for(i in 1:length(ras)){
        # Band specific output
        fname2 <- paste0( tools::file_path_sans_ext(fname), "__", names(ras)[i], raster::extension(fname))
        stars::write_stars(
          obj = ras,
          layer = 1:length(ras),
          dsn = fname2,
          type = ifelse(is.factor(ras[[1]]), "Byte", dtstars),
          NA_value = NA_real_,
          update = ifelse(file.exists(fname2), TRUE, FALSE),
          normalize_path = TRUE,
          progress = TRUE
        )
      }
      if(!is.Waiver(self$scenarios_migclim)){
        fname2 <- paste0( tools::file_path_sans_ext(fname), "__migclim", raster::extension(fname))
        writeNetCDF(ras_migclim, fname = fname, varName = "MigCLIM output", dt = dt)
      }
    } else if(type %in% 'feather'){
      assertthat::assert_that('feather' %in% installed.packages()[,1],
                              msg = 'Feather package not installed!')
      fname <- paste0( tools::file_path_sans_ext(fname), "__migclim", ".feather")
      feather::write_feather(ras, path = fname)
      if(!is.Waiver(self$scenarios_migclim)){
        fname2 <- paste0( tools::file_path_sans_ext(fname), "__migclim", raster::extension(fname))
        feather::write_feather(ras, path = fname)
      }
    }
    invisible()
  }
)
