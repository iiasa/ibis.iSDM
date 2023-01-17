#' @include utils.R waiver.R bdproto.R
NULL

#' @export
if (!methods::isClass("PredictorDataset")) methods::setOldClass("PredictorDataset")
NULL

#' PredictorDataset prototype description
#'
#' @name PredictorDataset-class
#' @aliases PredictorDataset
#' @family bdproto
#' @keywords bdproto
NULL

#' @export
PredictorDataset <- bdproto(
  "PredictorDataset",
  id           = character(0),
  data         = new_waiver(),
  # Printing function
  print = function(self){
    # Getting names and time periods if set
    nn <- name_atomic(self$get_names(), "predictors")
    # Get Time dimension if existing
    tt <- self$get_time()
    if(!(is.Waiver(tt) || is.null(tt))) tt <- paste0(range(tt),collapse = " <> ") else tt <- NULL
    message(paste0(self$name(),':',
                   '\n Name(s):  ',nn,
                   ifelse(!is.null(tt), paste0("\n Timeperiod:  ", tt), "")
                   )
    )
  },
  # Return name
  name = function(self){
    'Predictor dataset'
  },
  # Get Id
  id = function(self){
    self$id
  },
  # Get names of data
  get_names = function(self){
    names(self$get_data())
  },
  # Get data
  get_data = function(self, df = FALSE, na.rm = TRUE, ...){
    if(df) {
      if(any(is.factor(self$data))){
        # Bugs for factors, so need
        out <- self$data[] |> as.data.frame()
        out[,which(is.factor(self$data))] <- factor( out[,which(is.factor(self$data))] ) # Reformat factors variables
        cbind(raster::coordinates(self$data), out ) # Attach coordinates and return
      } else {
        raster::as.data.frame(self$data, xy = TRUE, na.rm = na.rm, ...)
      }
    } else self$data
  },
  # Get time dimension
  get_time = function(self, ...){
    # Get data
    d <- self$get_data()
    if(is.Waiver(d)) return(new_waiver())
    if(!inherits(d, 'stars')){
      # Try and get a z dimension from the raster object
      raster::getZ(d)
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
  # Get Projection
  get_projection = function(self){
    assertthat::assert_that(is.Raster(self$data) || inherits(self$data,'stars'))
    sf::st_crs(self$data)
  },
  # Clip the predictor dataset by another dataset
  crop_data = function(self, pol){
    assertthat::assert_that(is.Raster(self$data) || inherits(self$data,'stars'),
                            inherits(pol, 'sf'),
                            all( unique(sf::st_geometry_type(pol)) %in% c("POLYGON","MULTIPOLYGON") )
                            )
    self$data <- raster::crop(self$data, pol)
    invisible()
  },
  # Add a new Predictor dataset to this collection
  set_data = function(self, x, value){
    assertthat::assert_that(assertthat::is.string(x),
                            is.Raster(value),
                            is_comparable_raster(self$get_data(), value))
    bdproto(NULL, self, data = addLayer(self$get_data(), value))
  },
  # Remove a specific Predictor by name
  rm_data = function(self, x) {
    assertthat::assert_that(is.vector(x) || is.character(x),
                            all(x %in% names(self$get_data()))
                            )
    # Match indices
    ind <- match(x, self$get_names())
    # Overwrite predictor dataset
    self$data <- raster::dropLayer(self$get_data(), ind)
    invisible()
  },
  # Print input messages
  show = function(self) {
    self$print()
  },
  # Collect info statistics with optional decimals
  summary = function(self, digits = 2) {
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
          summary(stars:::as.data.frame.stars(d))
        )
      } else {
        # Assume raster
        return(
          round(
            raster::summary( d ), digits = digits
          )
        )
      }
    }
    rm(d)
  },
  # Has derivates?
  has_derivates = function(self){
    if(inherits(self$get_data(),'Raster'))
     return(
       length( grep("hinge__|bin__|quad__|thresh__", self$get_names() ) ) > 0
     )
    else
     return( NULL )
  },
  # Number of Predictors in object
  length = function(self) {
    if(inherits(self$get_data(),'Raster'))
      raster::nlayers(self$get_data())
    else
      ncol(self$get_data)
  },
  # Basic Plotting function
  plot = function(self){
    # Plot the predictors
    par.ori <- par(no.readonly = TRUE)
    raster::plot( self$data, col = ibis_colours[['viridis_cividis']] )
    par(par.ori)
  }
)
