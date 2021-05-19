#' Add a range of a species as predictor to a distribution object
#'
#' As options allow specifying including the range either as binary or distance
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param range A [`sf`] object with the range for the target feature
#' @param method [`character`] describing how the range should be included (binary | distance)
#' @param distance_max Numeric threshold on the maximum distance (Default: NULL)
#' @param priors A [`PriorList-class`] object. Default is set to NULL which uses default prior assumptions
#' @name add_range_predictor
NULL

#' @name add_range_predictor
#' @rdname add_range_predictor
#' @exportMethod add_range_predictor
#' @export
methods::setGeneric(
  "add_range_predictor",
  signature = methods::signature("x", "range", "method"),
  function(x, range, method = 'distance', distance_max = NULL, priors = NULL) standardGeneric("add_range_predictor"))

#' Function for when distance raster is directly supplied (precomputed)
#' @name add_range_predictor
#' @rdname add_range_predictor
#' @usage \S4method{add_range_predictor}{BiodiversityDistribution, raster}(x, range)
methods::setMethod(
  "add_range_predictor",
  methods::signature(x = "BiodiversityDistribution", range = "RasterLayer"),
  function(x, range, method = 'precomputed_range', priors = NULL) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            inherits(range, 'Raster')
    )
    names(range) <- method
    # Add as predictor
    if(is.Waiver(x$predictors)){
      x <- add_predictors(x, env = range,transform = 'none',derivates = 'none', priors)
    } else {
      x$predictors$set_data('range_distance', range)
      if(!is.null(priors)) {
        # FIXME: Ideally attempt to match varnames against supplied predictors vis match.arg or similar
        assertthat::assert_that( all( priors$varnames() %in% names(range) ) )
        x <- x$set_priors(priors)
      }
    }
    return(x)
  }
)

#' @name add_range_predictor
#' @rdname add_range_predictor
#' @usage \S4method{add_range_predictor}{BiodiversityDistribution,sf, vector}(x, range, method)
methods::setMethod(
  "add_range_predictor",
  methods::signature(x = "BiodiversityDistribution", range = "sf", method = "character"),
  function(x, range, method = 'distance', distance_max = NULL, priors = NULL ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.character(method),
                            inherits(range, 'sf'),
                            method %in% c('binary','distance'),
                            is.null(distance_max) || is.numeric(distance_max),
                            is.null(priors) || inherits(priors,'PriorList')
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
      dis <- raster::mask(dis,
                          x$background[unlist( st_intersects(st_buffer(range,0), x$background) ),]
      )
      # If max distance is specified
      if(!is.null(distance_max)) dis[dis > distance_max] <- NA # Set values above threshold to NA
      # Convert to relative for better scaling
      dis <- 1 - (dis / cellStats(dis,'max'))
      names(dis) <- 'range_distance'
    }

    # If priors have been set, save them in the distribution object
    if(!is.null(priors)) {
      # FIXME: Ideally attempt to match varnames against supplied predictors vis match.arg or similar
      assertthat::assert_that( all( priors$varnames() %in% names(dis) ) )
      x <- x$set_priors(priors)
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

#' Add a range of a species as offset to a distribution object
#'
#' As options allow specifying including the range either as binary or using distance
#' as parameter following Merow et al. (2017)
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param range A [`sf`] object with the range for the target feature
#' @param method [`character`] describing how the range should be included (binary | distance)
#' @param distance_max Numeric threshold on the maximum distance (Default: 150000 [m])
#' @name add_range_offset
NULL

#' @name add_range_offset
#' @rdname add_range_offset
#' @exportMethod add_range_offset
#' @export
methods::setGeneric(
  "add_range_offset",
  signature = methods::signature("x", "range", "method"),
  function(x, range, method = 'distance', name = 'range_distance', distance_max = 150000) standardGeneric("add_range_offset"))

#' Function for when raster is directly supplied (precomputed)
#' @name add_range_offset
#' @rdname add_range_offset
#' @usage \S4method{add_range_offset}{BiodiversityDistribution, raster}(x, range)
methods::setMethod(
  "add_range_offset",
  methods::signature(x = "BiodiversityDistribution", range = "RasterLayer"),
  function(x, range, name = 'range_distance') {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            inherits(range, 'Raster')
    )
    names(range) <- name
    # Add as a new offset
    x <- x$set_offset(range)
    return(x)
  }
)

#' @name add_range_offset
#' @rdname add_range_offset
#' @usage \S4method{add_range_offset}{BiodiversityDistribution,sf, vector}(x, range, method)
methods::setMethod(
  "add_range_offset",
  methods::signature(x = "BiodiversityDistribution", range = "sf", method = "character"),
  function(x, range, method = 'distance', name = 'range_distance', distance_max = 150000 ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.character(method),
                            inherits(range, 'sf'),
                            method %in% c('binary','distance'),
                            is.null(distance_max) || is.numeric(distance_max)
    )
    # Reproject if necessary
    if(st_crs(range) != sf::st_crs(x$background)) range <- sf::st_transform(range, sf::st_crs(x$background))

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

