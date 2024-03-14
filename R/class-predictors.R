#' @include waiver.R
NULL

if (!methods::isClass("PredictorDataset")) methods::setOldClass("PredictorDataset")

#' PredictorDataset class description
#'
#' @description
#' This class describes the PredictorDataset and is used to store covariates within.
#'
#' @keywords classes
#'
#' @name PredictorDataset-class
NULL

#' @rdname PredictorDataset-class
#' @export
PredictorDataset <- R6::R6Class(
  "PredictorDataset",
  public = list(
    #' @field id The id for this collection as [`character`].
    #' @field data A predictor dataset usually as [`SpatRaster`].
    #' @field name A name for this object.
    #' @field timeperiod A timeperiod field
    id           = character(0),
    data         = new_waiver(),
    name         = character(0),
    timeperiod   = new_waiver(),

    #' @description
    #' Initializes the object and creates an empty list
    #' @param id The id for this collection as [`character`].
    #' @param data A predictor dataset usually as [`SpatRaster`].
    #' @param ... Any other parameters found.
    #' @return NULL
    initialize = function(id, data, ...){
      assertthat::assert_that(
        is.Id(id) || is.character(id)
      )
      self$name <- 'Biodiversity data'
      self$id <- id
      self$data <- data
      # Get Dots and save too
      dots <- list(...)
      for(el in names(dots)){
        self[[el]] <- dots[[el]]
      }
    },

    #' @description
    #' Print the names and properties of all Biodiversity datasets contained within
    #' @param format A [`logical`] flag on whether a message should be printed.
    #' @return A message on screen
    print = function(format = TRUE){
      # Getting names and time periods if set
      nn <- name_atomic(self$get_names(), "predictors")
      # Get Time dimension if existing
      tt <- self$get_time()
      if(!(is.Waiver(tt) || is.null(tt))) tt <- paste0(range(tt),collapse = " <> ") else tt <- NULL
      m <- paste0(self$get_name(),':',
                  '\n Name(s):  ',nn,
                  ifelse(!is.null(tt), paste0("\n Timeperiod:  ", tt), "")
      )
      if(format) message(m) else return( m )
    },

    #' @description
    #' Return name of this object
    #' @return Default [`character`] name.
    get_name = function(){
      'Predictor dataset'
    },

    #' @description
    #' Get Id of this object
    #' @return Default [`character`] name.
    get_id = function(){
      return(self$id)
    },

    #' @description
    #' Get names of data
    #' @return [`character`] names of the data value.
    get_names = function(){
      names(self$get_data())
    },

    #' @description
    #' Alias for get_names
    #' @return [`character`] names of the data value.
    get_predictor_names = function(){
      names(self$get_data())
    },

    #' @description
    #' Get a specific dataset
    #' @param df [`logical`] on whether data is to be returned as [`data.frame`].
    #' @param na.rm [`logical`] if \code{NA} is to be removed from data.frame.
    #' @param ... Any other parameters passed on.
    #' @return A [`SpatRaster`] or [`data.frame`].
    get_data = function(df = FALSE, na.rm = TRUE, ...){
      assertthat::assert_that(
        is.logical(df), is.logical(na.rm)
      )
      if(df) {
        # For SpatRaster
        if(is.Raster(self$data)){
          if(any(is.factor(self$data))){
            # Bugs out for factors, so need
            terra::as.data.frame(self$data, xy = TRUE, na.rm = FALSE, ...)
          } else {
            terra::as.data.frame(self$data, xy = TRUE, na.rm = na.rm, ...)
          }
        } else {
          # For scenario files
          as.data.frame( self$data )
        }
      } else self$data
    },

    #' @description
    #' Get time dimension of object.
    #' @param ... Any other parameters passed on.
    #' @return A [`vector`] with the time dimension of the dataset.
    get_time = function(...){
      # Get data
      d <- self$get_data()
      if(is.Waiver(d)) return(new_waiver())
      if(!inherits(d, 'stars')){
        # Try and get a z dimension from the raster object
        terra::time(d)
      } else {
        # Get dimensions
        o <- stars::st_dimensions(d)
        # Take third entry as the one likely to be the time variable
        return(
          to_POSIXct(
            stars::st_get_dimension_values(d, names(o)[3], center = TRUE)
          )
        )
      }
    },

    #' @description
    #' Get Projection
    #' @return A [`vector`] with the geographical projection of the object.
    get_projection = function(){
      assertthat::assert_that(is.Raster(self$data) || inherits(self$data,'stars'))
      sf::st_crs(self$data)
    },

    #' @description
    #' Get Resolution
    #' @return A [`numeric`] [`vector`] with the spatial resolution of the data.
    get_resolution = function(){
      assertthat::assert_that(is.Raster(self$data) || inherits(self$data,'stars'))
      if(is.Raster(self$data)){
        terra::res(self$data)
      } else {
        stars::st_res(self$data)
      }
    },

    #' @description
    #' Utility function to clip the predictor dataset by another dataset
    #' @param pol A [`sf`] object used for cropping the data
    #' @return Invisibile TRUE
    crop_data = function(pol){
      assertthat::assert_that(is.Raster(self$data) || inherits(self$data,'stars'),
                              inherits(pol, 'sf'),
                              all( unique(sf::st_geometry_type(pol)) %in% c("POLYGON","MULTIPOLYGON") )
      )
      if(is.Raster(self$data)){
        self$data <- terra::crop(self$data, pol)
      } else {
        # Scenario
        suppressWarnings(
          suppressMessages(
            self$data <- sf::st_crop(self$data, pol)
          )
        )
      }
      invisible(self)
    },

    #' @description
    #' Utility function to mask the predictor dataset by another dataset
    #' @param mask A \code{SpatRaster} or `sf` object.
    #' @param inverse A `logical` flag if the inverse should be masked instead.
    #' @param ... Any other parameters passed on to masking.
    #' @return Invisible
    mask = function(mask, inverse = FALSE, ...){
      # Check whether prediction has been created
      prediction <- self$get_data(df = FALSE)
      if(!is.Waiver(prediction)){
        # If mask is sf, rasterize
        if(inherits(mask, 'sf')){
          mask <- terra::rasterize(mask, prediction)
        }
        # Check that mask aligns
        if(!terra::compareGeom(prediction, mask)){
          mask <- terra::resample(mask, prediction, method = "near")
        }
        # Now mask and save
        prediction <- terra::mask(prediction, mask, inverse = inverse, ...)

        # Save data
        self$data <- prediction
        invisible(self)
      }
    },

    #' @description
    #' Add a new Predictor dataset to this collection
    #' @param x [`character`] of the new name to be stored.
    #' @param value A new [`SpatRaster`] object.
    #' @return This object
    set_data = function(x, value){
      assertthat::assert_that(assertthat::is.string(x),
                              is.Raster(value),
                              is_comparable_raster(self$get_data(), value))
      self$data <- suppressWarnings( c(self$get_data(), value) )
      return(self)
    },

    #' @description
    #' Remove a specific Predictor by name
    #' @param x [`character`] of the predictor name to be removed.
    #' @return Invisible
    rm_data = function(x) {
      assertthat::assert_that(is.vector(x) || is.character(x),
                              all(x %in% names(self$get_data()))
      )
      # Match indices
      ind <- base::match(x, self$get_names())
      if(is.Raster(self$get_data() )){
        # Overwrite predictor dataset
        if(base::length(ind) == base::length(self$get_names())){
          self$data <- new_waiver()
        } else {
          self$data <- terra::subset(self$get_data(), ind, negate = TRUE)
        }
      } else {
        if(base::length(ind) == base::length(self$get_names())){
          self$data <- new_waiver()
        } else {
          suppressWarnings(
            self$data <- self$data[-ind]
          )
        }
      }
      invisible(self)
    },

    #' @description
    #' Alias for print method
    #' @return Invisible
    show = function() {
      self$print(format = FALSE)
    },

    #' @description
    #' Collect info statistics with optional decimals
    #' @param digits [`numeric`] Giving the rounding precision
    #' @return A [`data.frame`] summarizing the data.
    summary = function(digits = 2) {
      assertthat::assert_that(
        is.numeric(digits)
      )
      # Get data
      d <- self$get_data()
      if(is.Waiver(d)) return(NULL)
      # Need special handling if there any factors
      if(any(is.factor(self$get_data()))){
        out <- self$get_data()[] |> as.data.frame()
        out[,which(is.factor(self$data))] <- factor( out[,which(is.factor(self$data))] ) # Reformat factors variables
        summary(out, digits = digits)
      } else {
        if(inherits(d, 'stars')){
          return(
            summary( as.data.frame(d) )
          )
        } else {
          # Assume raster
          return(
            terra::summary( d, digits = digits)
          )
        }
      }
      rm(d)
    },

    #' @description
    #' Indication if there are any predictors that are derivates of outers
    #' @seealso [predictor_derivate()]
    #' @return A [`logical`] flag.
    has_derivates = function(){
      return(
        base::length( grep("hinge__|bin__|quad__|thresh__", self$get_names() ) ) > 0
      )
    },

    #' @description
    #' Number of Predictors in object
    #' @return A [`numeric`] estimate
    length = function() {
      if(inherits(self$get_data(),'SpatRaster'))
        terra::nlyr(self$get_data())
      else
        base::length(self$get_data())
    },

    #' @description
    #' Number of cells or values in object
    #' @return A [`numeric`] estimate
    ncell = function() {
      if(inherits(self$get_data(),'SpatRaster'))
        terra::ncell(self$get_data())
      else
        terra::ncell(self$get_data()) |> as.numeric()
    },

    #' @description
    #' Basic Plotting function
    #' @return A graphical interpretation of the predictors in this object.
    plot = function(){
      # Plot the predictors
      par.ori <- graphics::par(no.readonly = TRUE)
      if(is.Raster(self$data)){
        terra::plot( self$data, col = ibis_colours[['viridis_cividis']] )
      } else {
        # Assume stars scenario files
        stars:::plot.stars(self$data, col = ibis_colours[['viridis_cividis']])
      }
      graphics::par(par.ori)
    }
  ),

  # Any private entries
  private = list(
    finalize = function() {
    }
  )
)
