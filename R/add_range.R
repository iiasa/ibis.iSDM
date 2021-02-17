#' Add a range of a species as predictor to a distribution object
#'
#' As options allow specifying including the range either as binary or distance
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param range A [`sf`] object with the range for the target feature
#' @param method [`character`] describing how the range should be included (binary | distance)
#' @param distance_max Numeric threshold on the maximum distance (Default: NULL)
#' @name add_range_predictor
NULL

#' @name add_range_predictor
#' @rdname add_range_predictor
#' @exportMethod add_range_predictor
#' @export
methods::setGeneric(
  "add_range_predictor",
  signature = methods::signature("x", "range", "method"),
  function(x, range, method, distance_max = NULL) standardGeneric("add_range_predictor"))

#' @name add_range_predictor
#' @rdname add_range_predictor
#' @usage \S4method{add_range_predictor}{BiodiversityDistribution,sf, vector}(x, range, method)
methods::setMethod(
  "add_range_predictor",
  methods::signature(x = "BiodiversityDistribution", range = "sf", method = "character"),
  function(x, range, method = 'distance', distance_max = NULL ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.character(method),
                            inherits(range, 'sf'),
                            method %in% c('binary','distance'),
                            is.null(distance_max) || is.numeric(distance_max)
    )
    # Reproject if necessary
    if(st_crs(range) != st_crs(x$background)) range <- st_transform(range, st_crs(x$background))

    # Template raster for background
    if(!is.Waiver(x$predictors)){
      temp <- emptyraster(x$predictors$get_data())
    } else {
      # TODO: Eventually make this work better
      message('CAREFUL - This might not work without predictors already in the model.')
      temp <- raster::raster(extent(x$background),resolution = 1)
    }

    # Rasterize the range
    range$id <- 1:nrow(range) # Assign an id if not already existing
    if( 'fasterize' %in% installed.packages()[,1] ){
      ras_range <- fasterize::fasterize(range, temp, field = 'id')
    } else {
      ras_range <- raster::rasterize(range, temp,field = 'id')
    }

    # -------------- #
    if(methid == 'binary'){
      dis <- ras_range
      names(dis) <- 'binary_range'
    } else if(method == 'distance'){
      # TODO: The below can be much more sophisticated.
      # - For instance adding a exponential decay
      # Calculate the linear distance
      dis <- raster::distance(ras_range)
      dis <- raster::mask(dis, x$background)
      # Set areas not intersecting with range to 0
      dis <- raster::mask(dis,
                          x$background[unlist( st_intersects(st_buffer(range,0), x$background) ),]
      )
      # If max distance is specified
      if(!is.null(distance_max)) dis[dis > distance_max] <- NA # Set values above threshold to NA
      # Convert to relative for better scaling
      dis <- dis / cellStats(dis,'max')
      names(dis) <- 'range_distance'
    }

    # Add as predictor
    if(is.Waiver(x$predictors)){
      x <- add_predictors(x,env = dis,transform = 'none',derivates = 'none')
    } else {
      x$predictors$set_data('range_distance', dis)
    }

    return(x)
  }
)
