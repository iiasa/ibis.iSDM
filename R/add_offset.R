#' Add a range of a species as offset to a distribution object
#'
#' As options allow specifying including the range either as binary or using distance
#' as parameter following Merow et al. (2017)
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param range A [`sf`] object with the range for the target feature
#' @param method [`character`] describing how the range should be included (binary | distance)
#' @param name [`character`] A name of the species if the offset is not specified
#' @param distance_max Numeric threshold on the maximum distance (Default: 150000 [m])
#' @name add_offset_range
NULL

#' @name add_offset_range
#' @rdname add_offset_range
#' @exportMethod add_offset_range
#' @export
methods::setGeneric(
  "add_offset_range",
  signature = methods::signature("x", "range", "method"),
  function(x, range, method = 'distance', name = 'range_distance', distance_max = 150000) standardGeneric("add_offset_range"))

#' Function for when raster is directly supplied (precomputed)
#' @name add_offset_range
#' @rdname add_offset_range
#' @usage \S4method{add_offset_range}{BiodiversityDistribution, raster}(x, range)
methods::setMethod(
  "add_offset_range",
  methods::signature(x = "BiodiversityDistribution", range = "RasterLayer"),
  function(x, range, method = 'range_distance') {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.Raster(range),
                            is.character(method)
    )
    # Messager
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding range offset...')

    # Check that background and range align, otherwise raise error
    if(compareRaster(range, x$background,stopiffalse = FALSE)){
      warning('Supplied range does not align with background! Aligning them now...')
      range <- alignRasters(range, x$background, method = 'bilinear', func = mean, cl = FALSE)
    }
    names(range) <- method

    # Add as a new offset
    x <- x$set_offset(range)
    return(x)
  }
)

#' @name add_offset_range
#' @rdname add_offset_range
#' @usage \S4method{add_offset_range}{BiodiversityDistribution,sf, vector}(x, range, method)
methods::setMethod(
  "add_offset_range",
  methods::signature(x = "BiodiversityDistribution", range = "sf", method = "character"),
  function(x, range, method = 'distance', name = 'range_distance', distance_max = 150000 ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.character(method),
                            inherits(range, 'sf'),
                            method %in% c('binary','distance'),
                            is.null(distance_max) || is.numeric(distance_max)
    )
    # Messager
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding range offset...')

    # Reproject if necessary
    if(st_crs(range) != sf::st_crs(x$background)) range <- sf::st_transform(range, sf::st_crs(x$background))

    # Template raster for background
    if(!is.Waiver(x$predictors)){
      temp <- emptyraster(x$predictors$get_data())
    } else {
      # TODO: Eventually make this work better
      myLog('[Setup]','red','CAREFUL - This might not work without predictors already in the model.')
      temp <- raster::raster(extent(x$background),resolution = 1)
    }

    # Rasterize the range
    range$id <- 1:nrow(range) # Assign an id if not already existing
    if( 'fasterize' %in% installed.packages()[,1] ){
      ras_range <- fasterize::fasterize(range, temp, field = NULL)
    } else {
      ras_range <- raster::rasterize(range, temp,field = NULL)
    }

    # -------------- #
    if(method == 'binary'){
      dis <- ras_range
      names(dis) <- 'binary_range'
    } else if(method == 'distance'){
      # TODO: The below can be much more sophisticated.
      # - For instance adding a exponential decay
      # Calculate the linear distance
      dis <- raster::distance(ras_range)
      dis <- raster::mask(dis, x$background)
      # Set areas not intersecting with range to 0
      suppressMessages( suppressWarnings( range <- sf::st_buffer(range, 0)) )

      suppressMessages(
        dis <- raster::mask(dis,
                        x$background[unlist( st_intersects(range, x$background) ),]
        )
      )
      # If max distance is specified
      if(!is.null(distance_max)) dis[dis > distance_max] <- NA # Set values above threshold to NA
      # Convert to relative for better scaling
      dis <- 1 - predictor_transform(dis,option = 'norm') #1 - (dis / cellStats(dis,'max') )
      names(dis) <- name
    }

    # Add as a new offset
    x <- x$set_offset(dis)
    return(x)
  }
)

