if (!methods::isClass("BiodiversityScenario")) methods::setOldClass("BiodiversityScenario")

#' Class for a biodiversity scenario from a trained model
#'
#' @description
#' Base [`R6`] class for any biodiversity scenario objects. Serves as
#' container that supplies data and functions to other [`R6`] classes and functions.
#'
#' @keywords classes
#'
#' @name BiodiversityScenario-class
NULL

#' @rdname BiodiversityScenario-class
#' @export
BiodiversityScenario <- R6::R6Class(
  "BiodiversityScenario",
  public = list(
    #' @field modelobject A name of the model for projection.
    #' @field modelid An id of the model used for projection.
    #' @field limits A [`sf`] object used to constraint the prediction.
    #' @field predictors A predictor object for projection.
    #' @field constraints Any constraints set for projection.
    #' @field scenarios The resulting [`stars`] objects.
    modelobject = new_waiver(),
    modelid = new_waiver(),
    limits = new_waiver(),
    predictors = new_waiver(),
    constraints = new_waiver(),
    scenarios = new_waiver(),

    #' @description
    #' Initializes the object and creates an empty list
    #' @return NULL
    initialize = function(){
    },

    #' @description
    #' Print the names and properties of all scenarios.
    #' @return A message on screen
    print = function() {
      # Check that model exists
      fit <- self$get_model()
      timeperiod <- self$get_timeperiod()
      if(!is.Waiver(timeperiod)){
        dur <- round(as.numeric(difftime(timeperiod[2], timeperiod[1], unit="weeks"))/52.25,1)
        if(dur == 0) dur <- "< 1"
      }
      # Get set predictors and time period
      pn = ifelse(is.Waiver(self$get_predictor_names()),'None',name_atomic(self$get_predictor_names(), "predictors"))
      tp = ifelse(is.Waiver(timeperiod),'None',
                  paste0(
                    paste0( timeperiod,collapse = ' -- '),
                    ' (',dur,' years)'
                  )
      )
      # Constrains
      cs <- self$get_constraints()
      if(!is.Waiver(cs)) cs <- vapply(cs, function(x) x$method, character(1))
      # Thresholds
      tr <- self$get_threshold()

      # Any other simulation outputs modules
      simmods <- self$get_simulation()

      message(paste0(
        ifelse(is.Waiver(self$limits),"Spatial-temporal scenario:","Spatial-temporal scenario (limited):"),
        '\n  Used model: ',ifelse(is.Waiver(fit) || isFALSE(fit), text_red('None'), fit$get_name() ),
        "\n --------- ",
        "\n  Predictors:     ", pn,
        "\n  Time period:    ", tp,
        ifelse(!is.Waiver(cs)||!is.Waiver(tr), "\n --------- ", ""),
        ifelse(is.Waiver(cs),"", paste0("\n  Constraints:    ", text_green(paste(paste0(names(cs),' (',cs,')'),collapse = ', ')) ) ),
        ifelse(is.Waiver(tr),"", paste0("\n  Threshold:      ", round(tr[1], 3),' (',names(tr[1]),')') ),
        ifelse(is.Waiver(simmods),"", paste0("\n  Simulations:    ", text_green(paste(paste0(names(simmods),' (',simmods[[1]][[1]],')'),collapse = ', ')) ) ),
        "\n --------- ",
        "\n  Scenarios fitted: ", ifelse(is.Waiver(self$scenarios),text_yellow('None'), text_green('Yes'))
      )
      )
    },

    #' @description
    #' Verify that set Model exist and check self-validity
    #' @return Invisible
    verify = function(){
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
      invisible(self)
    },

    #' @description
    #' Show the name of the Model
    #' @return Model objectname
    show = function() {
      if(is.character(self$modelobject)){
        return( self$modelobject )
      } else {
        return( fit$model$runname )
      }
    },

    #' @description
    #' Get projection of the projection.
    #' @return A [`sf`] object with the geographic projection
    get_projection = function(){
      return( self$predictors$get_projection() )
    },

    #' @description
    #' Get resultion of the projection.
    #' @return A [`numeric`] indication of the resolution.
    get_resolution = function(){
      return( self$predictors$get_resolution() )
    },

    #' @description
    #' Get the actual model used for projection
    #' @param copy A [`logical`] flag on whether a deep copy should be created.
    #' @return A DistributionModel object.
    get_model = function(copy = FALSE){
      if(is.Waiver(self$modelobject)) return( new_waiver() )
      else {
        if(inherits(self$modelobject, "DistributionModel")){
          obj <- self$modelobject
        } else {
          if(!exists(self$modelobject)) return( FALSE ) else {
            obj <- get(self$modelobject)
          }
        }
      }
      if(copy) obj <- obj$clone(deep = TRUE)
      return( obj )
    },

    #' @description
    #' Get provided projection limits if set.
    #' @return A [`sf`] object or NULL.
    get_limits = function(){
      if(is.Waiver(self$limits)) return(NULL)
      return(self$limits)
    },

    #' @description
    #' Get names of predictors for scenario object.
    #' @return A [`character`] vector with the names.
    get_predictor_names = function() {
      if(is.Waiver(self$predictors)) return(self$predictors)
      if(inherits(self$predictors, "PredictorDataset")) {
        self$predictors$get_names()
      } else {
        stop("Feature data is of an unrecognized class")
      }
    },

    #' @description
    #' Get time period of projection.
    #' @param what [`character`] on whether full time period or just the range is to be returned.
    #' @return A time period from start to end.
    get_timeperiod = function(what = "range"){
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

    #' @description
    #' Get constrains for model
    #' @return A [`list`] with the constraints within the scenario.
    get_constraints = function(){
      return( self$constraints )
    },

    #' @description
    #' Get thresholds if specified.
    #' @seealso [threshold()]
    #' @return A [`list`] with method and value for the threshold.
    get_threshold = function(){
      if('threshold' %notin% names(self)) return( new_waiver() )
      return( self$threshold )
    },

    #' @description
    #' Duplicate function for internal consistency to return threshold
    #' @seealso [threshold()]
    #' @return A [`list`] with method and value for the threshold.
    get_thresholdvalue = function(){
      if('threshold' %notin% names(self)) return( new_waiver() )
      return( self$threshold )
    },

    #' @description
    #' Apply a new threshold to the projection.
    #' @note
    #' This sets the threshold method internally to \code{'fixed'}.
    #' @param tr A [`numeric`] value with the new threshold.
    #' @return This object.
    apply_threshold = function(tr = new_waiver()){
      # Assertions
      if(is.Waiver(tr)) assertthat::assert_that( is.numeric(self$threshold),
                                                 msg = 'No threshold value found.')
      assertthat::assert_that( !is.Waiver(self$scenarios),
                               msg = 'No scenarios found.')
      # Get prediction and threshold
      sc <- self$get_data()
      if(is.Waiver(tr)){
        tr <- self$threshold
      } else {
        # Reassign name as method
        names(tr) <- "fixed"
      }
      if(!is.Waiver(tr)){
        # Select only suitability layer
        if("threshold" %in% names(sc)) sc <- sc |> dplyr::select(suitability)
        # reclassify to binary
        new <- sc
        new[new < tr[[1]]] <- 0; new[new >= tr[[1]]] <- 1
        names(new) <- 'threshold'
        # Add to scenario object
        sc <- c(sc, new)
        # Format new threshold object
        self$scenarios <- sc
        self$threshold <- tr
      }
      return(self)
    },

    #' @description
    #' Set new predictors to this object.
    #' @param x [`PredictorDataset-class`] object to be supplied.
    #' @return This object.
    set_predictors = function(x){
      assertthat::assert_that(inherits(x, "PredictorDataset"))
      self$predictors <- x
      return(self)
    },

    #' @description
    #' Set new constrains
    #' @param x A [`SpatRaster`] object to be added as as constraint.
    #' @return This object.
    set_constraints = function(x){
      if(!is.Waiver(self$get_constraints())){
        cr <- self$get_constraints()
        # FIXME: Remove duplicates
        self$constraints <- c(cr, x)
      } else {
        self$constraints <- x
      }
      return(self)
    },

    #' @description
    #' Get simulation options and parameters if gound
    #' @return A [`list`] with the parameters.
    get_simulation = function(){
      if('simulation' %notin% names(self)) return( new_waiver() )
      return( self$simulation )
    },

    #' @description
    #' Set simulation objects.
    #' @param x new simulation entries and options as [`list`] to be set.
    #' @return This object.
    set_simulation = function(x){
      # We only take a single simulation so far
      self$simulation <- x
      return(self)
    },

    #' @description
    #' Get Predictors from the object.
    #' @return A predictor dataset.
    get_predictors = function(){
      return(self$predictors)
    },

    #' @description
    #' Remove predictors from the object.
    #' @param names A [`character`] vector with names
    #' @return This object.
    rm_predictors = function(names){
      if(is.Waiver(self$predictors) || is.null(self$predictors)) return(NULL)
      if(missing(names)){
        names <- self$get_predictor_names() # Assume all names
      }
      assertthat::assert_that(
        is.character(names) || assertthat::is.scalar(names) || is.vector(names)
      )
      # Get predictor collection
      prcol <- self$predictors$clone(deep = TRUE)
      # Set the object
      prcol$rm_data(names)
      if(base::length(prcol$get_predictor_names())==0) prcol <- new_waiver()
      self$predictors <- prcol
      invisible(self)
    },

    #' @description
    #' Get scenario predictions or any other data
    #' @param what A [`character`] vector with names of what
    #' @return This object.
    get_data = function(what = "scenarios"){
      return(self[[what]])
    },

    #' @description
    #' Set new data in object.
    #' @param x A new data object measuing scenarios.
    #' @return This object.
    set_data = function(x){
      # Get projected value
      ff <- self$scenarios
      # Set the object
      ff[["scenarios"]] <- x
      self$scenarios <- ff
      return( self )
    },

    #' @description
    #' Plot the predictions made here.
    #' @param what A [`character`] describing the layers to be plotted.
    #' @param which A [`numeric`] subset to any specific time steps.
    #' @param ... Any other parameters passed on.
    #' @return A graphical representation
    plot = function(what = "suitability", which = NULL, ...){
      assertthat::assert_that(is.character(what))

      if(is.Waiver(self$get_data())){
        if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Scenario]','red','No scenarios found')
        invisible(self)
      } else {
        assertthat::assert_that(what %in% names( self$get_data() ),
                                msg = paste(what, "not found in scenario projections?"))
        # Get unique number of data values. Surely there must be an easier val
        vals <- self$get_data()[what] |> stars:::pull.stars() |> as.vector() |> unique()
        vals <- base::length(stats::na.omit(vals))
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

    #' @description
    #' Convenience function to plot thresholds if set
    #' @param which A [`numeric`] subset to any specific time steps.
    #' @return A graphical representation
    plot_threshold = function(which = NULL){
      # Check that baseline and scenario thresholds are all there
      if(!( 'threshold' %in% attributes(self$get_data())$names )) return(new_waiver())

      self$plot(what = "threshold", which = which)
    },

    #' @description
    #' Plot Migclim results if existing.
    #' @return A graphical representation
    plot_migclim = function(){
      # Get scenarios
      mc <- self$get_data("scenarios_migclim")
      if(is.Waiver(mc)) return(mc)

      # Otherwise plot the raster
      ras <- mc$raster

      # Colour coding from MigClim::MigClim.plot
      rstVals <- sort(terra::unique(ras))
      negativeNb <- base::length(which(rstVals < 0))
      positiveNb <- base::length(which(rstVals > 1 & rstVals < 30000))
      zeroExists <- any(rstVals == 0)
      oneExists <- any(rstVals == 1)
      unilimtedExists <- any(rstVals == 30000)
      Colors <- rep("yellow", negativeNb)
      if(zeroExists) Colors <- c(Colors, "grey94")
      if (oneExists) Colors <- c(Colors, "black")
      Colors <- c(Colors, rainbow(positiveNb, start = 0, end = 0.4))
      if (unilimtedExists) Colors <- c(Colors, "pink")

      # Plot 0 - Cells that have never been occupied and are unsuitable habitat at
      # the end of the simulation
      # 1 - Cells that belong to the species' initial
      # distribution and that have remained occupied during the entire simulation.
      # 1 < value < 30 000 - determine the dispersal step during which it was
      # colonized. E.g. 101 is first dispersal even in first step 30 0000 -
      # Potentially suitable cells that remained uncolonized <0 - Negative values
      # indicate cells that were once occupied but have become decolonized. Code
      # as for colonization
      dev.new(width = 7, height = 7 * ((ymax(ras) - ymin(ras))/(xmax(ras) - xmin(ras))))
      terra::plot(ras, col = Colors, breaks = c(min(rstVals) - 1, rstVals), legend = FALSE,
                  main = "Newly colonized and stable habitats")
    },

    #' @description
    #' Plot animation of scenarios if possible
    #' @note
    #' This requires the \code{"gganimate"} package.
    #' @param what A [`character`] describing the layers to be plotted.
    #' @param fname An optional filename to write the result.
    #' @return A graphical representation
    plot_animation = function(what = "suitability", fname = NULL){
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

    #' @description
    #' Plot relative change between baseline and projected thresholds
    #' @note
    #' This requires a set [`threshold()`] to the scenario object.
    #' @param position Which layer to be plotted
    #' @param variable A [`character`] of the variable to be plotted
    #' @param plot [`logical`] flag on whether to plot the results or return the object.
    #' @return A graphical representation or [`SpatRaster`].
    plot_relative_change = function(position = NULL, variable = 'mean', plot = TRUE){
      # Default position is the last one
      assertthat::assert_that(is.null(position) || is.numeric(position) || is.character(position),
                              is.character(variable))
      # Threshold
      obj <- self$get_model()
      thresh_reference <- grep('threshold',obj$show_rasters(),value = T)
      # If there is more than one threshold only use the one from variable
      if(base::length(thresh_reference)>1) {
        warning('More than one baseline threshold. Using the first one.')
        thresh_reference <- grep(variable, thresh_reference,value = T)[1]
      }
      # Check that baseline and scenarios are all there
      assertthat::assert_that(
        !is.Waiver(self$get_data()),
        'threshold' %in% attributes(self$get_data())$names,
        base::length(thresh_reference) >0 & is.character(thresh_reference),
        is.Raster( self$get_model()$get_data('prediction') ),
        msg = "Threshold not found!"
      )

      # Not get the baseline raster
      baseline <- self$get_model()$get_data(thresh_reference)
      # Get only the variable
      if(terra::nlyr(baseline)>1) baseline <- terra::subset(baseline, grep(variable, names(baseline)))
      # And the last scenario prediction
      scenario <- self$get_data()['threshold']
      time <- stars::st_get_dimension_values(scenario, which = 3) # 3 assumed to be time band
      if(is.numeric(position)) position <- time[position]
      if(is.null(position)) position <- time[base::length(time)]
      final <- scenario |>
        stars:::filter.stars(time == position) |>
        terra::rast()
      if(is.na(terra::crs(final))) terra::crs(final) <- terra::crs(baseline)
      # -- #
      if(!inherits(final, 'SpatRaster')) final <- final[[1]] # In case it is a rasterbrick or similar
      if(!is_comparable_raster(baseline, final)) final <- alignRasters(final, baseline, cl = FALSE) # In case they somehow differ?

      # Calculate overlays
      diff_f <- terra::lapp(c(baseline, final) |> as.numeric(), fun = function(x, y){x + y * 2})
      diff_f <- as.factor(diff_f)
      # 0 = Unsuitable | 1 = Loss | 2 = Gain | 3 = stable
      rat <- levels(diff_f)[[1]]
      rat <- merge.data.frame(rat, data.frame(ID = seq(0,3), diff = c("Unsuitable", "Loss", "Gain", "Stable")),all = TRUE)
      levels(diff_f) <- rat
      diff_f <- terra::mask(diff_f, baseline)
      rm(baseline, final)

      # Plot
      if(plot){
        # Colours
        cols <- c("Unsuitable" = "gray92", "Loss" = "#DE646A", "Gain" = "cyan3", "Stable" = "gray60")

        # Convert to raster
        diff_ff <- terra::as.data.frame(diff_f, xy = TRUE)
        names(diff_ff)[3] <- "Change"
        diff_ff$Change <- factor(diff_ff$Change, levels = names(cols))

        terra::plot(diff_f, col = cols,
                    buffer = TRUE,
                    type = "classes", #legend = "right",
                    levels = names(cols), all_levels = TRUE,
                    main = paste0('Change between\n baseline and ', position))
      } else {
        # Return
        return(diff_f)
      }
    },

    #' @description
    #' Summarize the change in layers between timesteps
    #' @param layer A [`character`] of the variable to be plotted
    #' @param plot [`logical`] flag on whether to plot the results or return the coefficients.
    #' @param relative [`logical`] on coefficients to be converted to relative change.
    #' @return Summarized coefficients as [`data.frame`]
    summary = function(layer = "threshold", plot = FALSE, relative = FALSE){
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
        time <- stars::st_get_dimension_values(scenario, which = 3) # Assuming band 3 is the time dimension
        assertthat::assert_that(!is.na(sf::st_crs(scenario)), msg = "Scenario not correctly projected.")
        # HACK: Add area to stars
        ar <- stars:::st_area.stars(scenario)
        # Get the unit
        ar_unit <- units::deparse_unit(ar$area)
        new <- (scenario |> terra::rast()) * (ar |> terra::rast())
        terra::time(new) <- time
        # Convert to scenarios to data.frame
        df <- stars:::as.data.frame.stars(stars:::st_as_stars(new)) |> (\(.) subset(., stats::complete.cases(.)))()
        # Rename
        names(df)[3:4] <- c("band", "area")
        # --- #
        # Now calculate from this data.frame several metrics related to the area and change in area
        df <- df |> dplyr::group_by(x,y) |> dplyr::mutate(id = dplyr::cur_group_id()) |>
          dplyr::ungroup() |> dplyr::select(-x,-y) |>
          dplyr::mutate(area = dplyr::if_else(is.na(area), 0, area)) |> # Convert missing data to 0
          dplyr::arrange(id, band)
        df$area <- units::as_units(df$area, units::as_units(ar_unit))  # Set Units
        # Convert to km2 and remove units as this causes issues with dplyr
        df$area <- units::set_units(df$area, "km2") |> units::drop_units()

        # Total amount of area occupied for a given time step
        out <- df |> dplyr::group_by(band) |> dplyr::summarise(area_km2 = sum(area, na.rm = TRUE))
        out$totarea <- terra::global((new[[1]]>=0) * (ar |> terra::rast()), "sum", na.rm = TRUE)[,1]
        if(units::deparse_unit(units::as_units(ar_unit)) == "m2") {
          out$totarea <- out$totarea / 1e6
          out <- dplyr::rename(out, totarea_km2 = totarea)
        }

        # Total amount of area lost / gained / stable since previous time step
        totchange_occ <- df |>
          dplyr::group_by(id)  |>
          dplyr::mutate(change = (area - dplyr::lag(area)) ) |> dplyr::ungroup() |>
          (\(.) subset(., stats::complete.cases(.)))()
        o <- totchange_occ |> dplyr::group_by(band) |>
          dplyr::summarise(totchange_stable_km2 = sum(area[change == 0]),
                           totchange_gain_km2 = sum(change[change > 0]),
                           totchange_loss_km2 = sum(change[change < 0]))
        out <- out |> dplyr::left_join(o, by = "band")

        if(relative == TRUE){
          # Finally calculate relative change to baseline (first entry) for all entries where this is possible
          relChange <- function(v, fac = 100) (((v- v[1]) / v[1]) * fac)
          out <- subset(out, select = c("band", "area_km2", "totarea_km2"))
          out[,c("area_km2")] <- apply( out[,c("area_km2")], 2, relChange)
        }
      } else {
        # Get the scenario predictions and from there the thresholds
        scenario <- self$get_data()['suitability']
        # If it has a time dimension also get that
        if( length(stars::st_dimensions(scenario)) >2){
          times <- stars::st_get_dimension_values(scenario, which = 'band')
        } else {
          times <- NULL
          stop("Summary without time dimension not yet implemented!")
        }
        # Get area
        ar <- stars:::st_area.stars(scenario)
        # Get the unit
        ar_unit <- units::deparse_unit(ar$area)

        # TODO: Check whether one could not simply multiply with
        # area (poisson > density, binomial > suitable area)
        mod <- self$get_model()
        scenario <- scenario * ar
        out <- summarise_projection(scenario, fun = "mean", relative = relative)
      }

      if(plot){
        if( 'threshold' %in% attributes(self$get_data())$names ){
          if(utils::hasName(out,"band")) out <- dplyr::rename(out, "time" = "band")
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
            ggplot2::geom_ribbon(fill = "grey85") +
            ggplot2::geom_line(size = 2) +
            ggplot2::labs(x = "Time", y = expression(Area(km^2)), title = "Relative suitable habitat")
        }
      }
      return(out)
    },

    #' @description
    #' Summarize before-after change of first and last layer.
    #' @note
    #' This requires set [`threshold`] prior to projection.
    #' @return Summarized coefficients as [`data.frame`]
    summary_beforeafter = function(){
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

    #' @description
    #' Calculate slopes across the projection
    #' @param what A [`character`] with layer to be plotted (default: \code{"suitability"}).
    #' @param oftype [`character`] of the output type.
    #' @return A plot of the scenario slopes
    plot_scenarios_slope = function(what = 'suitability', oftype = "stars"){
      self$calc_scenarios_slope(what = what, plot = TRUE, oftype = oftype)
      invisible()
    },

    #' @description
    #' Calculate slopes across the projection
    #' @param what A [`character`] with layer to be plotted (default: \code{"suitability"}).
    #' @param plot [`logical`] flag on whether to plot the results or return the coefficients.
    #' @param oftype [`character`] of the output type.
    #' @return A [`SpatRaster`] layer or [`stars`] object.
    calc_scenarios_slope = function(what = 'suitability', plot = TRUE, oftype = "stars"){
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
        if(plot){
          suppressWarnings(
            stars:::plot.stars(out, breaks = "fisher", col = c(ibis_colours$divg_bluered[1:10],"grey85",ibis_colours$divg_bluered[11:20]))
          )
        }
      } else {
        out <- terra::rast(out)
        if(plot) terra::plot(out, col = c(ibis_colours$divg_bluered[1:10],"grey85",ibis_colours$divg_bluered[11:20]))
      }
      return(out)
    },

    #' @description
    #' Convenience function to mask all input projections.
    #' @param mask A \code{SpatRaster} or `sf` object.
    #' @param inverse A `logical` flag if the inverse should be masked instead.
    #' @param ... Any other parameters passed on.
    #' @return Invisible
    mask = function(mask, inverse = FALSE, ...){
      # Check whether prediction has been created
      projection <- self$get_data()
      if(!is.Waiver(projection)){
        # Make valid
        mask <- sf::st_make_valid(mask)
        # If mask is sf, rasterize
        if(!inherits(mask, 'sf')){
          mask <- terra::as.polygons(mask) |> sf::st_as_sf()
        }
        # Check that aligns
        if(sf::st_crs(projection) != sf::st_crs(mask) ){
          mask <- mask |> sf::st_transform(crs = sf::st_crs(projection))
        }
        # If there are multiple, ajoin
        if(nrow(mask)>1){
          mask <- sf::st_combine(mask) |> sf::st_as_sf()
        }

        # Now mask and save
        projection <-
          suppressMessages(
            suppressWarnings(
              projection[mask]
            )
          )

        # Save data
        self[["scenarios"]] <- projection

        invisible(self)
      }
    },

    #' @description
    #' Get centroids of projection layers
    #' @param patch A [`logical`] if centroid should be calculated weighted by values.
    #' @return Returns a [`sf`] object.
    get_centroid = function(patch = FALSE){
      assertthat::assert_that(
        is.logical(patch),
        inherits(self$get_data(), "stars"),
        msg = "This function only works if projections have been calculated!"
      )

      # Check if threshold is present
      rl <- names(self$get_data())
      if(length( grep('threshold',rl,value = TRUE) )>0){
        # Threshold present
        obj <- self$get_data( )[grep('threshold',rl,value = TRUE)]

        # Convert to list of SpatRasters
        ll <- stars_to_raster(obj)
        # Loop through each time step
        cent <- data.frame()
        for(step in names(ll)){
          ras <- ll[[step]]
          # Calculate centroid per time step
          poi <- raster_centroid(ras, patch)
          poi$time <- step
          cent <- rbind(cent, poi)
        }
      } else {
        # Get non-thresholded layer
        obj <- self$get_data( )[ grep('suitability',rl,value = TRUE) ]
        # Convert to list of SpatRasters
        ll <- stars_to_raster(obj)
        # Loop through each time step
        cent <- data.frame()
        for(step in names(ll)){
          ras <- ll[[step]]
          # Calculate centroid per time step,
          # patch to FALSE as this is non-sensical here
          poi <- raster_centroid(ras, patch = FALSE)
          poi$time <- step
          cent <- rbind(cent, poi)
        }
      }
      return(cent)
    },

    #' @description
    #' Save object as output somewhere
    #' @param fname An output filename as [`character`].
    #' @param type A format as [`character`]. Matched against a list of supported formats.
    #' @param dt The datatype used, such as float64
    #' @return Saved spatial prediction on drive.
    save = function(fname, type = 'tif', dt = 'FLT4S'){
      assertthat::assert_that(
        !missing(fname),
        is.character(fname),
        is.character(type),
        !is.Waiver(self$get_data()),
        is.character(dt)
      )
      # Match input types
      type <- match.arg(type, c('gtif','gtiff','tif',
                                'nc','ncdf', 'feather'), several.ok = FALSE)
      dt <- match.arg(dt, c('LOG1S','INT1S','INT1U',
                            'INT2S','INT2U','INT4S',
                            'INT4U','FLT4S','FLT8S'), several.ok = FALSE )

      if(file.exists(fname)) warning('Overwritting existing file...')

      # Respecify type if output filename has already been set
      if(gsub('\\.','', tools::file_ext(fname)) != type) type <- gsub('\\.','', tools::file_ext(fname))

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
          fname2 <- paste0( tools::file_path_sans_ext(fname), "__", names(ras)[i], tools::file_ext(fname))
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
          fname2 <- paste0( tools::file_path_sans_ext(fname), "__migclim", tools::file_ext(fname))
          writeGeoTiff(ras_migclim, fname = fname, dt = dt)
        }
      } else if(type %in% c('nc','ncdf')) {
        # Save as netcdf, for now in individual files
        for(i in 1:length(ras)){
          # Band specific output
          fname2 <- paste0( tools::file_path_sans_ext(fname), "__", names(ras)[i], tools::file_ext(fname))
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
          fname2 <- paste0( tools::file_path_sans_ext(fname), "__migclim", tools::file_ext(fname))
          writeNetCDF(ras_migclim, fname = fname, varName = "MigCLIM output", dt = dt)
        }
      } else if(type %in% 'feather'){
        assertthat::assert_that('feather' %in% utils::installed.packages()[,1],
                                msg = 'Feather package not installed!')
        fname <- paste0( tools::file_path_sans_ext(fname), "__migclim", ".feather")
        feather::write_feather(ras, path = fname)
        if(!is.Waiver(self$scenarios_migclim)){
          fname2 <- paste0( tools::file_path_sans_ext(fname), "__migclim", tools::file_ext(fname))
          feather::write_feather(ras, path = fname)
        }
      }
      invisible(self)
    }
  ),

  # Any private entries
  private = list(
    finalize = function() {
    }
  ),
  # Don't lock object, but support adding new member
  lock_objects = FALSE
)
