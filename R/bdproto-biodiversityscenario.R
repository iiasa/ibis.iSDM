#' @include utils.R bdproto.R
NULL

#' @export
if (!methods::isClass("BiodiversityScenario")) methods::setOldClass("BiodiversityScenario")
NULL

#' Prototype for a biodiversity scenario from a trained model
#'
#' @name BiodiversityScenario-class
#' @aliases BiodiversityScenario
#' @keywords bdproto
NULL

#' @export
BiodiversityScenario <- bdproto(
  "BiodiversityScenario",
  modelobject = new_waiver(), # The id of the model
  modelid = new_waiver(),
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
    cs <- self$get_constrains()
    # Thresholds
    tr <- self$get_threshold()

    message(paste0('Spatial-temporal scenario:',
                   '\n  Used model: ',ifelse(is.Waiver(fit) || isFALSE(fit), text_red('None'), class(fit)[1] ),
                   "\n --------- ",
                   "\n  Predictors:     ", pn,
                   "\n  Time period:    ", tp,
                   ifelse(!is.Waiver(cs)||!is.Waiver(tr), "\n --------- ", ""),
                   ifelse(is.Waiver(cs),"", paste0("\n  Constraints:      ", cs) ),
                   ifelse(is.Waiver(tr),"", paste0("\n  Threshold:      ", round(tr[1], 3),' (',names(tr[1]),')') ),
                   "\n --------- ",
                   "\n  Scenarios fitted: ", ifelse(is.Waiver(self$scenarios),text_yellow('None'), text_green('Yes'))
      )
    )
  },
  # Verify that set Model exist
  verify = function(self){
    assertthat::validate_that( !is.Waiver(self$modelobject),
                               !is.Waiver(self$modelid),
                               exists(self$modelobject) )
    # Get Model object and check that ID is correct
    x <- get(self$modelobject)
    assertthat::validate_that(x$id == self$modelid)
    invisible()
  },
  # Get Model
  get_model = function(self){
    if(is.Waiver(self$modelobject)) return( new_waiver() )
      else
        if(!exists(self$modelobject)) return( FALSE )
          else return( get(self$modelobject) )
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
  get_timeperiod = function(self){
    if(is.Waiver(self$predictors)) return(self$predictors)
    if(inherits(self$predictors, "PredictorDataset")) {
      return(
        c( min(self$predictors$timeperiod), max(self$predictors$timeperiod) )
      )
    }
  },
  # Get constrains for model
  get_constrains = function(self){
    # TODO: To be implemented
    return(new_waiver())
  },
  # Get thresholds if specified
  get_threshold = function(self){
    if('threshold' %notin% names(self)) return( new_waiver() )
    return( self$threshold )
  },
  # Apply specific threshold
  apply_threshold = function(self){
    # Assertions
    assertthat::assert_that( is.numeric(self$threshold), msg = 'No threshold value found.')
    assertthat::assert_that( !is.Waiver(self$scenarios), msg = 'No scenarios found.')
    # Get prediction and threshold
    sc <- self$get_scenarios()
    tr <- self$threshold
    # reclassify to binary
    sc[sc < tr] <- 0; sc[sc >= tr] <- 1
    names(sc) <- 'presence'
    return(sc)
  },
  # Show the name of the Model
  show = function(self) {
    self$modelobject
  },
  # Set Predictors
  set_predictors = function(self, x){
    assertthat::assert_that(inherits(x, "PredictorDataset"))
    bdproto(NULL, self, predictors = x)
  },
  # Get Predictors
  get_predictors = function(self){
    return(self$predictors)
  },
  # Get scenario predictions
  get_scenarios = function(self){
    return(self$scenarios)
  },
  # Calculate slopes
  calc_scenarios_slope = function(self, what = 'suitability', plot = TRUE){
    if(is.Waiver(self$get_scenarios())) return( new_waiver() )
    assertthat::assert_that(what %in% attributes(self$get_scenarios())$names )

    oo <- self$get_scenarios()[what]
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
    if(plot) stars:::plot.stars(out, breaks = "equal", col = ibis_colours$divg_bluered)
    return(out)
  },
  # Plot the prediction
  plot = function(self, what = "suitability",...){
    # FIXME: More plotting options would be good
    if(is.Waiver(self$get_scenarios())){
      if(getOption('ibis.setupmessages')) myLog('[Scenario]','red','No scenarios found')
      invisible()
    } else {
      # Get unique number of data values. Surely there must be an easier val
      vals <- self$get_scenarios()[what] %>% stars:::pull.stars() %>% as.vector() %>% unique() %>% length()
      if(vals>2) col <- ibis_colours$sdm_colour else col <- c('grey25','coral')
      stars:::plot.stars( self$get_scenarios()[what], breaks = "equal", col = col )
    }
  },
  # Plot animation of scenarios
  plot_animation = function(self, what = "suitability", fname = NULL){
    assertthat::assert_that(!is.Waiver(self$get_scenarios()) )
    check_package('gganimate')
    # Get scenarios
    obj <- self$get_scenarios()[what]

    # Make the animation plot
    g <- ggplot2::ggplot() +
      stars::geom_stars(data = obj, downsample = c(1,1,0)) +
      ggplot2::coord_equal() +
      ggplot2::theme_bw() +
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
  plot_relative_change = function(self, position = NULL, variable = 'mean'){
    # Default position is the last one
    assertthat::assert_that(is.null(position) || is.numeric(position) || is.character(position),
                            is.character(variable))
    # Threshold
    thresh_reference <- grep('threshold',modf$show_rasters(),value = T)
    # If there is more than one threshold only use the one from variable
    if(length(thresh_reference)>1) thresh_reference <- grep(variable, thresh_reference,value = T)
    # Check that baseline and scenarios are all there
    assertthat::assert_that(
      !is.Waiver(self$get_scenarios()),
      'threshold' %in% attributes(self$get_scenarios())$names,
      length(thresh_reference) >0 & is.character(thresh_reference),
      is.Raster( self$get_model()$get_data('prediction') )
    )

    # Not get the baseline raster
    baseline <- self$get_model()$get_data(thresh_reference)
    # And the last scenario prediction
    scenario <- self$get_scenarios()['threshold']
    time <- stars::st_get_dimension_values(scenario,which = 'band')
    if(is.numeric(position)) position <- time[position]
    if(is.null(position)) position <- time[length(time)]
    final <- scenario %>%
      stars:::filter.stars(band == position) %>%
      as('Raster')
    raster::projection(final) <- raster::projection(baseline)
    # -- #
    if(!inherits(final, 'RasterLayer')) final <- final[[1]] # In case it is a rasterbrick or similar
    if(!compareRaster(baseline, final,stopiffalse = FALSE)) final <- alignRasters(final, baseline, cl = FALSE) # In case they somehow differ?
    final[final > 0 & !is.na(final)] <- 2
    # Sum up the layers. 1 == Presence earlier | 2 == Presence in future | 3 == Presence in bot
    diff_f <- as.factor( (baseline + final)+1 )
    rat <- levels(diff_f)[[1]]
    rat <- merge.data.frame(rat, data.frame(ID = seq(1,4), diff = c("Absent", "Extinction", "Colonisation", "Stable")))
    levels(diff_f) <- rat
    diff_f <- raster::mask(diff_f, baseline)

    # Plot
    if('rasterVis' %in% installed.packages()[,1]){
      rasterVis::levelplot(diff_f,
                           margin = F,
                           scales = list(draw=TRUE),
                           col.regions = c("grey75","coral","cyan3","grey25"),
                           main = paste0('Change between baseline and ', position)
                           )
    } else {
      # Convert to stars for plotting otherwise
      # FIXME: Stars plotting bugs out if there are fewer than 4 classes
      diff_f <- stars::st_as_stars(diff_f, att = 'diff');names(diff_f) <- 'Change'
      cols <- c("grey75","coral","cyan3","grey25")
      stars:::plot.stars(diff_f, axes = TRUE,key.pos = 1, border = NA, extent = sf::st_bbox(baseline),
                         main = paste0('Change between baseline and ', position),
                         col = cols)
    }
    return(diff_f)

  },
  # Save object
  save = function(self, fname, type = 'gtif', dt = 'FLT4S'){
    assertthat::assert_that(
      !missing(fname),
      is.character(fname),
      type %in% c('gtif','gtiff','tif','nc','ncdf', 'feather'),
      'fits' %in% self$ls(),
      dt %in% c('LOG1S','INT1S','INT1U','INT2S','INT2U','INT4S','INT4U','FLT4S','FLT8S')
    )
    type <- tolower(type)

    # Respecify type if output filename has already been set
    if(gsub('\\.','',raster::extension(fname)) != type) type <- gsub('\\.','',raster::extension(fname))

    # Get raster file in fitted object
    cl <- sapply(self$scenarios, class)
    ras <- self$scenarios[[grep('raster', cl,ignore.case = T)]]

    # Check that no-data value is not present in ras
    assertthat::assert_that(any(!cellStats(ras,min) <= -9999),msg = 'No data value -9999 is potentially in prediction!')

    if(file.exists(fname)) warning('Overwritting existing file...')
    if(type %in% c('gtif','gtiff','tif')){
      # Save as geotiff
      writeGeoTiff(ras, fname = fname, dt = dt)
    } else if(type %in% c('nc','ncdf')) {
      # Save as netcdf
      # TODO: Potentially change the unit descriptions
      writeNetCDF(ras, fname = fname, varName = 'iSDM prediction', varUnit = "",varLong = "")
    } else if(type %in% 'feather'){
      assertthat::assert_that('feather' %in% installed.packages()[,1],
                              msg = 'Feather package not installed!')
      feather::write_feather(ras, path = fname)
    }
    invisible()
  }
)
