#' Specify a spatial explicit offset
#'
#' @description
#' Including offsets is another option to integrate spatial prior information
#' in linear and additive regression models. Offsets shift the intercept of
#' the regression fit by a certain amount. Although only one offset can be added
#' to a regression model, it is possible to combine several spatial-explicit estimates into
#' one offset by calculating the sum of all spatial-explicit layers.
#'
#' @details
#' This function allows to set any specific offset to a regression model. The offset
#' has to be provided as spatial [`RasterLayer`] object. This function simply adds the layer to
#' a [`distribution()`] object.
#' **Note that any transformation of the offset (such as \code{log}) has do be done externally!**
#'
#' If the layer is range and requires additional formatting, consider using the
#' function [`add_offset_range()`] which has additional functionalities such such distance
#' transformations.
#'
#' @note
#' Since offsets only make sense for linear regressions (and not for instance
#' regression tree based methods such as [engine_bart]), they do not work for all engines.
#' Offsets specified for non-supported engines are ignored during the estimation
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param layer A [`sf`] or [`RasterLayer`] object with the range for the target feature.
#' @param add [`logical`] specifiying whether new offset is to be added. Setting
#' this parameter to \code{FALSE} replaces the current offsets with the new one (Default: \code{TRUE}).
#' @param ... Other parameters or arguments (currently not supported)
#' @references
#' * Merow, C., Allen, J.M., Aiello-Lammens, M., Silander, J.A., 2016. Improving niche and range estimates with Maxent and point process models by integrating spatially explicit information. Glob. Ecol. Biogeogr. 25, 1022–1036. https://doi.org/10.1111/geb.12453
#' @family offset
#' @keywords prior, offset
#' @examples
#' \dontrun{
#'  x <- distribution(background) %>%
#'    add_predictors(covariates) %>%
#'    add_offset(nicheEstimate)
#' }
#' @name add_offset
NULL

#' @name add_offset
#' @rdname add_offset
#' @exportMethod add_offset
#' @export
methods::setGeneric(
  "add_offset",
  signature = methods::signature("x", "layer"),
  function(x, layer, add = TRUE) standardGeneric("add_offset"))

#' @name add_offset
#' @rdname add_offset
#' @usage \S4method{add_offset}{BiodiversityDistribution, raster}(x, layer)
methods::setMethod(
  "add_offset",
  methods::signature(x = "BiodiversityDistribution", layer = "RasterLayer"),
  function(x, layer, add = TRUE) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.Raster(layer),
                            is.logical(add)
    )
    # Messenger
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding spatial explicit offset...')
    ori.name <- names(layer)

    # Check for infinite values
    assertthat::assert_that(
      all( is.finite(cellStats(layer, "range")) ),
      msg = "Infinite values found in the layer (maybe log of 0?)."
    )

    # Check that background and range align, otherwise raise error
    if(compareRaster(layer, x$background,stopiffalse = FALSE)){
      warning('Supplied layer does not align with background! Aligning them now...')
      layer <- alignRasters(layer, x$background, method = 'bilinear', func = mean, cl = FALSE)
      names(layer) <- ori.name
    }

    # Check whether an offset exists already
    if(!is.Waiver(x$offset) && add){
      # Add to current object
      of <- x$offset
      layer <- raster::resample(layer, of, method = 'bilinear', func = mean, cl = FALSE)
      names(layer) <- ori.name # In case the layer name got lost
      of <- raster::stack(of) |> raster::addLayer(layer)
      x <- x$set_offset(of)
    } else {
      # Add as a new offset
      x <- x$set_offset(layer)
    }
    return(x)
  }
)

#' Specify a spatial explicit offset as bias
#'
#' @description
#' Including offsets is another option to integrate spatial prior information
#' in linear and additive regression models. Offsets shift the intercept of
#' the regression fit by a certain amount. Although only one offset can be added
#' to a regression model, it is possible to combine several spatial-explicit estimates into
#' one offset by calculating the sum of all spatial-explicit layers.
#'
#' @details
#' This functions emulates the use of the [`add_offset()`] function, however applies an inverse
#' transformation to remove the provided layer from the overall offset.
#' So if for instance a offset is already specified (such as area), this function
#' removes the provided \code{bias.layer} from it via \code{"offset(log(off.area)-log(bias.layer))"}
#'
#' **Note that any transformation of the offset (such as \code{log}) has do be done externally!**
#'
#' If a generic offset is added, consider using the [`add_offset()`] function. If the layer is a expert-based range and
#' requires additional formatting, consider using the
#' function [`add_offset_range()`].
#'
#' @inheritParams add_offset
#' @references
#' * Merow, C., Allen, J.M., Aiello-Lammens, M., Silander, J.A., 2016. Improving niche and range estimates with Maxent and point process models by integrating spatially explicit information. Glob. Ecol. Biogeogr. 25, 1022–1036. https://doi.org/10.1111/geb.12453
#' @family offset
#' @keywords prior, offset
#' @examples
#' \dontrun{
#'  x <- distribution(background) %>%
#'    add_predictors(covariates) %>%
#'    add_offset_bias(samplingBias)
#' }
#' @name add_offset_bias
NULL

#' @name add_offset_bias
#' @rdname add_offset_bias
#' @exportMethod add_offset_bias
#' @export
methods::setGeneric(
  "add_offset_bias",
  signature = methods::signature("x", "layer"),
  function(x, layer, add = TRUE) standardGeneric("add_offset_bias"))

#' @name add_offset_bias
#' @rdname add_offset_bias
#' @usage \S4method{add_offset_bias}{BiodiversityDistribution, raster}(x, layer)
methods::setMethod(
  "add_offset_bias",
  methods::signature(x = "BiodiversityDistribution", layer = "RasterLayer"),
  function(x, layer, add = TRUE) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.Raster(layer),
                            is.logical(add)
    )
    # Messenger
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding spatial explicit bias offset...')
    ori.name <- names(layer)

    # Check that background and range align, otherwise raise error
    if(compareRaster(layer, x$background,stopiffalse = FALSE)){
      warning('Supplied layer does not align with background! Aligning them now...')
      layer <- alignRasters(layer, x$background, method = 'bilinear', func = mean, cl = FALSE)
      names(layer) <- ori.name
    }
    # Since it is a bias offset and removal is equivalent to simple subtraction, multiply with *-1
    layer <- layer * -1

    # Check for infinite values
    assertthat::assert_that(
      all( is.finite(cellStats(layer, "range")) ),
      msg = "Infinite values found in the layer (maybe log of 0?)."
    )

    # Check whether an offset exists already
    if(!is.Waiver(x$offset) && add){
      # Add to current object
      of <- x$offset
      layer <- raster::resample(layer, of, method = 'bilinear', func = mean, cl = FALSE)
      names(layer) <- ori.name # In case the layer name got lost
      of <- raster::stack(of) |> raster::addLayer(layer)
      x <- x$set_offset(of)
    } else {
      # Add as a new offset
      x <- x$set_offset(layer)
    }
    return(x)
  }
)

#' Specify a expert-based species range as offset
#'
#' @description
#' This function has additional nuance compared to the more generic
#' [`add_offset()`], allowing to specify transformations and modifications.
#' It is customized specifically for expert-based ranges as offsets.
#'
#' @details
#' As options allow specifying including the range either as binary or using distance
#' as parameter following Merow et al. (2017). In this case the existing range needs
#' to be transformed first, for instance using the \pkg{bossMaps} package.
#' @inheritParams add_offset
#' @param method [`character`] describing how the range should be included (Options: \code{"binary"} | \code{"distance"}).
#' @param distance_max Numeric threshold on the maximum distance (Default: \code{150000} [m]).
#' @seealso [bossMaps]
#' @references
#' * Merow, C., Wilson, A.M., Jetz, W., 2017. Integrating occurrence data and expert maps for improved species range predictions. Glob. Ecol. Biogeogr. 26, 243–258. https://doi.org/10.1111/geb.12539
#' * Merow, C., Allen, J.M., Aiello-Lammens, M., Silander, J.A., 2016. Improving niche and range estimates with Maxent and point process models by integrating spatially explicit information. Glob. Ecol. Biogeogr. 25, 1022–1036. https://doi.org/10.1111/geb.12453
#' @keywords prior, offset
#' @family offset
#' @name add_offset_range
NULL

#' @name add_offset_range
#' @rdname add_offset_range
#' @exportMethod add_offset_range
#' @export
methods::setGeneric(
  "add_offset_range",
  signature = methods::signature("x", "layer", "method"),
  function(x, layer, method = 'distance', distance_max = 150000, add = TRUE) standardGeneric("add_offset_range"))

#' Function for when raster is directly supplied (precomputed)
#' @name add_offset_range
#' @rdname add_offset_range
#' @usage \S4method{add_offset_range}{BiodiversityDistribution, raster}(x, layer)
methods::setMethod(
  "add_offset_range",
  methods::signature(x = "BiodiversityDistribution", layer = "RasterLayer"),
  function(x, layer, method = 'range_distance', add = TRUE) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.Raster(layer),
                            is.character(method),
                            is.logical(add)
    )
    # Messager
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding range offset...')

    # Save name
    ori.name <- names(layer)

    # Check for infinite values
    assertthat::assert_that(
      all( is.finite(cellStats(layer, "range")) ),
      msg = "Infinite values found in the layer (maybe log of 0?)."
    )

    # Check that background and range align, otherwise raise error
    if(compareRaster(layer, x$background,stopiffalse = FALSE)){
      warning('Supplied range does not align with background! Aligning them now...')
      layer <- alignRasters(layer, x$background, method = 'bilinear', func = mean, cl = FALSE)
      names(layer) <- ori.name # In case the layer name got lost
    }

    # Check whether an offset exists already
    if(!is.Waiver(x$offset) && add){
      # Add to current object
      of <- x$offset
      layer <- raster::resample(layer, of, method = 'bilinear', func = mean, cl = FALSE)
      names(layer) <- ori.name # In case the layer name got lost
      of <- raster::stack(of) |> raster::addLayer(layer)
      x <- x$set_offset(of)
    } else {
      # Add as a new offset
      x <- x$set_offset(layer)
    }
    return(x)
  }
)

#' @name add_offset_range
#' @rdname add_offset_range
#' @usage \S4method{add_offset_range}{BiodiversityDistribution, sf, vector}(x, layer, method)
methods::setMethod(
  "add_offset_range",
  methods::signature(x = "BiodiversityDistribution", layer = "sf", method = "character"),
  function(x, layer, method = 'distance', distance_max = 150000, add = TRUE ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.character(method),
                            inherits(layer, 'sf'),
                            method %in% c('binary','distance'),
                            is.null(distance_max) || is.numeric(distance_max),
                            is.logical(add)
    )
    # Messager
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding range offset...')

    # Reproject if necessary
    if(st_crs(layer) != sf::st_crs(x$background)) layer <- sf::st_transform(layer, sf::st_crs(x$background))

    # Template raster for background
    if(!is.Waiver(x$predictors)){
      temp <- emptyraster(x$predictors$get_data())
    } else {
      # TODO: Eventually make this work better
      myLog('[Setup]','red','CAREFUL - This might not work without predictors already in the model.')
      temp <- raster::raster(extent(x$background),resolution = 1)
    }

    # Rasterize the range
    layer$id <- 1:nrow(layer) # Assign an id if not already existing
    if( 'fasterize' %in% installed.packages()[,1] ){
      ras_range <- fasterize::fasterize(layer, temp, field = NULL)
    } else {
      ras_range <- raster::rasterize(layer, temp,field = NULL)
    }

    # -------------- #
    if(method == 'binary'){
      dis <- ras_range
      names(dis) <- "range_binary"
    } else if(method == 'distance'){
      # TODO: The below can be much more sophisticated.
      # - For instance adding a exponential decay
      # Calculate the linear distance
      dis <- raster::distance(ras_range)
      dis <- raster::mask(dis, x$background)
      # Set areas not intersecting with range to 0
      suppressMessages( suppressWarnings( layer <- sf::st_buffer(layer, 0)) )

      suppressMessages(
        dis <- raster::mask(dis,
                        x$background[unlist( sf::st_intersects(layer, x$background) ),]
        )
      )
      # If max distance is specified
      if(!is.null(distance_max)) dis[dis > distance_max] <- NA # Set values above threshold to NA
      # Convert to relative for better scaling
      dis <- 1 - predictor_transform(dis, option = 'norm') #1 - (dis / cellStats(dis,'max') )
      names(dis) <- "range_distance"
    }

    # Check whether an offset exists already
    if(!is.Waiver(x$offset) && add){
      # Add to current object
      of <- x$offset
      ori.name <- names(dis)
      dis <- raster::resample(dis, of, method = 'bilinear', func = mean, cl = FALSE)
      names(dis) <- ori.name # In case the layer name got lost
      of <- raster::stack(of) |> raster::addLayer(dis)
      x <- x$set_offset(of)
    } else {
      # Add as a new offset
      x <- x$set_offset(dis)
    }
    return(x)
  }
)
