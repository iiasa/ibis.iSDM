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
#' @param add [`logical`] specifying whether new offset is to be added. Setting
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

#' Function to remove an offset
#'
#' @description
#' This is just a wrapper function for removing specified offsets from a [`BiodiversityDistribution-class`]) object.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param layer A `character` pointing to the specific layer to be removed. If set to \code{NULL}, then
#' all offsets are removed from the object.
#' @family offset
#' @keywords prior, offset, internal
#' @name rm_offset
NULL

#' @name rm_offset
#' @rdname rm_offset
#' @exportMethod rm_offset
#' @export
methods::setGeneric(
  "rm_offset",
  signature = methods::signature("x", "layer"),
  function(x, layer = NULL) standardGeneric("rm_offset"))

#' @name rm_offset
#' @rdname rm_offset
#' @usage \S4method{rm_offset}{BiodiversityDistribution, character}(x, layer)
methods::setMethod(
  "rm_offset",
  methods::signature(x = "BiodiversityDistribution", layer = "character"),
  function(x, layer = NULL) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.character(layer) || is.null(layer)
    )
    # If no offset can be found, just return proto object
    if(is.Waiver(x$offset)){ return(x) }

    # Messenger
    if(getOption('ibis.setupmessages')) myLog('[Setup]','yellow','Removing offsets.')

    offs <- x$get_offset()
    if(!is.null(layer)){
      assertthat::assert_that(layer %in% offs,
                              msg = paste0("Specified offset ", layer, "not found in the offset list."))
    }

    # Now remove the offset
    x$rm_offset()
  }
)

#### Bias offset ----
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
#' requires additional parametrization, consider using the
#' function [`add_offset_range()`] or the \code{bossMaps} R-package.
#'
#' @inheritParams add_offset
#' @param points An optional [`sf`] object with key points. The location of the points are then used to
#' calculate the probability that a cell has been sampled while accounting for area differences.
#' (Default: \code{NULL}).
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
  function(x, layer, add = TRUE, points = NULL) standardGeneric("add_offset_bias"))

#' @name add_offset_bias
#' @rdname add_offset_bias
#' @usage \S4method{add_offset_bias}{BiodiversityDistribution, raster}(x, layer)
methods::setMethod(
  "add_offset_bias",
  methods::signature(x = "BiodiversityDistribution", layer = "RasterLayer"),
  function(x, layer, add = TRUE, points = NULL) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.Raster(layer),
                            is.logical(add),
                            is.null(points) || inherits(points, 'sf')
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
    if(is.null(points)){
      # Since it is a bias offset and removal is equivalent to simple subtraction, multiply with *-1
      layer <- layer * -1
    } else {
      ## Count the number of records per cell
      tab <- raster::cellFromXY(layer, sf::st_coordinates(points))
      r <- emptyraster(layer)
      r[tab] <- layer[tab]
      r <- raster::mask(r, background)

      ## Make zeros a very small number otherwise issues with log(0).
      r[r[]==0] <- 1e-6
      suppressWarnings({ar <- raster::area(r)})

      ## Calculate the probability that a cell has been sampled while accounting for area differences in lat/lon
      off.bias <- (-log(1-exp(-r * ar)) - log(ar))
      names(off.bias) <- "off.bias"
      # Add bias as covariate
      layer <- off.bias
      ## NOTE: if area offset considered, use "+ offset(log(off.area)-log(off.bias))"
    }

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

#### Add a range as offset ----
#' Specify a expert-based species range as offset
#'
#' @description
#' This function has additional options compared to the more generic
#' [`add_offset()`], allowing customized options specifically for expert-based ranges as offsets or spatialized
#' polygon information on species occurrences.
#' If even more control is needed, the user is informed of the \pkg{bossMaps} package Merow et al. (2017). The \pkg{bossMaps} package
#' calculates - based on supplied point information - the probability of occurrences being inside vs outside the
#' range map and can thus be used as a method to 'improve' the mapping of a species range.
#'
#' @details
#' The output created by this function creates a [`RasterLayer`] to be added to a provided distribution object. Offsets
#' in regression models are likelihood specific as they are added directly to the overall estimate of \code{`y^hat`}.
#'
#' Note that all offsets created by this function are by default log-transformed before export. Background values
#' (e.g. beyond [`distance_max`]) are set to a very small constant (\code{1e-10}).
#'
#' @inheritParams add_offset
#' @param distance_max A [`numeric`] threshold on the maximum distance beyond the range that should be considered
#' to have a high likelihood of containing species occurrences (Default: \code{Inf} [m]). Can be set to \code{NULL} or \code{0}
#' to indicate that no distance should be calculated.
#' @param type A [`character`] denoting the type of model to which this offset is to be added. By default
#' it assumes a \code{'poisson'} distributed model and as a result the output created by this function will be log-transformed.
#' If however a \code{'binomial'} distribution is chosen, than the output will be \code{`logit`} transformed.
#' For integrated models leave at default.
#' @param presence_prop [`numeric`] giving the proportion of all records expected to be inside the range. By
#' default this is set to \code{0.9} indicating that 10% of all records are likely outside the range.
#' @param distance_clip [`logical`] as to whether distance should be clipped after the maximum distance (Default: \code{FALSE}).
#' @seealso [`bossMaps`]
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
  signature = methods::signature("x", "layer"),
  function(x, layer, distance_max = Inf, type = "poisson", presence_prop = 0.9, distance_clip = FALSE, add = TRUE) standardGeneric("add_offset_range"))

#' Function for when raster is directly supplied (precomputed)
#' @name add_offset_range
#' @rdname add_offset_range
#' @usage \S4method{add_offset_range}{BiodiversityDistribution, raster}(x, layer)
methods::setMethod(
  "add_offset_range",
  methods::signature(x = "BiodiversityDistribution", layer = "RasterLayer"),
  function(x, layer, add = TRUE) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.Raster(layer),
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
#' @usage \S4method{add_offset_range}{BiodiversityDistribution, sf}(x, layer)
methods::setMethod(
  "add_offset_range",
  methods::signature(x = "BiodiversityDistribution", layer = "sf"),
  function(x, layer, distance_max = Inf, type = "poisson", presence_prop = 0.9, distance_clip = FALSE, add = TRUE ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            inherits(layer, 'sf'),
                            is.null(distance_max) || is.numeric(distance_max) || is.infinite(distance_max),
                            is.numeric(presence_prop),
                            is.logical(distance_clip),
                            is.character(type),
                            is.logical(add)
    )
    # Match the type if set
    type <- match.arg(type, c("poisson", "binomial"), several.ok = FALSE)

    # Messenger
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding range offset...')

    # Reproject if necessary
    if(sf::st_crs(layer) != sf::st_crs(x$background)) layer <- sf::st_transform(layer, sf::st_crs(x$background))

    # If distance max is null, set to 0
    if(is.null(distance_max)) distance_max <- 0

    # Template raster for rasterization background
    if(!is.Waiver(x$predictors)){
      temp <- emptyraster(x$predictors$get_data())
    } else {
      # Try and guess an sensible background raster
      myLog('[Setup]','red',
      'CAREFUL - This might not work without predictors already in the model.
      Add offset after predictors')
      temp <- raster::raster(raster::extent(x$background),
                             resolution = diff(sf::st_bbox(x$background)[c(1,3)]) / 100,
                             crs = sf::st_crs(x$background))
    }

    # Check to make the entries valid
    if( any(!sf::st_is_valid(layer)) ){
      layer <- sf::st_make_valid(layer) # Check whether to make them valid
      if( any(!sf::st_is_valid(layer)) ){
        # If still has errors, combine
        layer <- layer |> sf::st_combine() |> sf::st_as_sf()
      }
    }

    # If layer has multiple entries join them
    if(nrow(layer)>1) layer <- layer |> sf::st_union() |> sf::st_as_sf()

    # Rasterize the range
    if( 'fasterize' %in% installed.packages()[,1] ){
      ras_range <- try({ fasterize::fasterize(layer, temp, field = NULL, background = NA) })
      if(inherits(ras_range,"try-error")){
        myLog('[Setup]','yellow','Fasterize package needs to be re-installed!')
        ras_range <- raster::rasterize(layer, temp, field = 1, background = NA)
      }
    } else {
      ras_range <- raster::rasterize(layer, temp, field = 1, background = NA)
    }
    # Calculate distance if required
    if(distance_max > 0){
      # Calculate a distance raster
      dis <- raster::gridDistance(ras_range, origin = 1)
      # If max distance is specified
      if(distance_clip && is.finite(distance_max)){
        dis[dis > distance_max] <- NA # Set values above threshold to a very small constant
      }
      # Inverse of distance
      if(is.infinite(distance_max)) distance_max <- cellStats(dis,"max")
      # ---- #
      alpha <- 1 / (distance_max / 4 ) # Divide by 4 for a quarter in each direction
      # Grow baseline raster by using an exponentially weighted kernel
      dis <- raster::calc(dis, fun = function(x) exp(-alpha * x))
      # Set the remaining ones to very small constant
      dis[is.na(dis)] <- 1e-10 # Background values
      dis <- raster::mask(dis, x$background)

    } else {
      dis <- ras_range
      dis[is.na(dis)] <- 1e-10 # Background values
      dis <- raster::mask(dis, x$background)
    }

    # Inside I want all X across the entire area for the PPMs,
    # indicating a lambda per area of at least X/A (per unit area) within the range
    suppressWarnings( ar <- raster::area(ras_range) ) # Calculate area
    pres <- 1 + ( ( raster::cellStats(ar * ras_range, "sum") / raster::cellStats(ar, "sum")) * (presence_prop) )
    abs <- 1 + ( ( raster::cellStats(ar * ras_range, "sum") / raster::cellStats(ar, "sum")) * (1-presence_prop) )
    # Now set all values inside the range to pres and outside to abs
    ras_range[ras_range == 1] <- pres
    ras_range[is.na(ras_range)] <- abs
    # Multiply with distance layer
    ras_range <- ras_range * dis
    # Normalize the result by dividing by the sum
    ras_range <- ras_range / raster::cellStats(ras_range, "sum", na.rm = TRUE)

    # -------------- #
    # Log transform
    ras_range  <- log(ras_range)
    # Rescaling does not affect relative differences.
    ras_range <- raster::scale(ras_range, scale = F)
    names(ras_range) <- "range_distance"

    assertthat::assert_that(
      is.finite( raster::cellStats(ras_range, "max") ),
      msg = "Range offset has infinite values. Check parameters!"
    )

    # Check whether an offset exists already
    if(!is.Waiver(x$offset) && add){
      # Add to current object
      of <- x$offset
      ori.name <- names(ras_range)
      ras_range <- raster::resample(ras_range, of, method = 'bilinear', func = mean, cl = FALSE)
      names(ras_range) <- ori.name # In case the layer name got lost
      of <- raster::stack(of) |> raster::addLayer(ras_range)
      x <- x$set_offset(of)
    } else {
      # Add as a new offset
      x <- x$set_offset(ras_range)
    }
    return(x)
  }
)

#' Specify elevational preferences as offset
#'
#' @description
#' This function implements the elevation preferences offset defined in Ellis‐Soto et al. (2021).
#' The code here was adapted from the Supporting materials script.
#' @details
#' Specifically this functions calculates a continuous decay and decreasing probability of a species to occur
#' from elevation limits. It requires a [`RasterLayer`] with elevation information.
#' A generalized logistic transform (aka Richard's curve) is used to calculate decay from the suitable elevational
#' areas, with the [`rate`] parameter allowing to vary the steepness of decline.
#'
#' Note that all offsets created by this function are by default log-transformed before export. In addition
#' this function also mean-centers the output as recommended by Ellis-Soto et al.
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param elev A [`RasterLayer`] with the elevation for a given background.
#' @param pref A [`numeric`] vector of length \code{2} giving the lower and upper bound of known elevational preferences.
#' Can be set to \code{Inf} if unknown.
#' @param rate A [`numeric`] for the rate used in the offset (Default: \code{.0089}). This parameter specifies the
#' decay to near zero probability at elevation above and below the expert limits.
#' @param add [`logical`] specifying whether new offset is to be added. Setting
#' this parameter to \code{FALSE} replaces the current offsets with the new one (Default: \code{TRUE}).
#' @references
#' * Ellis‐Soto, D., Merow, C., Amatulli, G., Parra, J.L., Jetz, W., 2021. Continental‐scale 1 km hummingbird diversity derived from fusing point records with lateral and elevational expert information. Ecography (Cop.). 44, 640–652. https://doi.org/10.1111/ecog.05119
#' * Merow, C., Allen, J.M., Aiello-Lammens, M., Silander, J.A., 2016. Improving niche and range estimates with Maxent and point process models by integrating spatially explicit information. Glob. Ecol. Biogeogr. 25, 1022–1036. https://doi.org/10.1111/geb.12453
#' @keywords prior, offset
#' @family offset
#' @name add_offset_elevation
NULL

#' @name add_offset_elevation
#' @rdname add_offset_elevation
#' @exportMethod add_offset_elevation
#' @export
methods::setGeneric(
  "add_offset_elevation",
  signature = methods::signature("x", "elev", "pref"),
  function(x, elev, pref, rate = .0089, add = TRUE) standardGeneric("add_offset_elevation"))


#' @name add_offset_elevation
#' @rdname add_offset_elevation
#' @usage \S4method{add_offset_elevation}{BiodiversityDistribution, raster, numeric}(x, elev, pref)
methods::setMethod(
  "add_offset_elevation",
  methods::signature(x = "BiodiversityDistribution", elev = "RasterLayer", pref = "numeric"),
  function(x, elev, pref, rate = .0089, add = TRUE) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.Raster(elev),
                            is.numeric(pref),
                            length(pref)==2,
                            is.logical(add)
    )
    # Messenger
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding elevation offset...')

    # Check for infinite values
    assertthat::assert_that(
      all( is.finite(cellStats(elev, "range")) ),
      msg = "Infinite values found in the layer (maybe log of 0?)."
    )

    # Check that background and range align, otherwise raise error
    if(compareRaster(elev, x$background,stopiffalse = FALSE)){
      warning('Supplied range does not align with background! Aligning them now...')
      elev <- alignRasters(elev, x$background, method = 'bilinear', func = mean, cl = FALSE)
    }
    # Generalized logistic transform (aka Richard's curve) function from bossMaps.
    genLogit <- function(x, lower = 0, upper = 1, rate = 0.04, skew = 0.2, shift = 0){
      upper - ((upper - lower)/((1 + exp(-rate * (x - shift)))^(1/skew)))
    }

    # ---- #
    # if(getOption("ibis.runparallel")) raster::beginCluster(n = getOption("ibis.nthread"))
    # Now calculate the elevation offset by projecting the values onto the elevation layer
    # max avail > min expert
    tmp.elev1 = -1 * (elev - pref[1])
    tmp.elev1.1 = raster::calc(tmp.elev1, function(x) genLogit(x, 1,100, rate, .2 ))
    # min avail < max expert
    tmp.elev2 = elev - pref[2]
    tmp.elev2.1 = raster::calc(tmp.elev2, function(x) genLogit(x, 1,100, rate, .2))
    # Combine both and calculate the minimum
    elev.prior = min( raster::stack(tmp.elev1.1, tmp.elev2.1))
    rm(tmp.elev1,tmp.elev1.1,tmp.elev2,tmp.elev2.1) # clean up

    # Normalize the result by dividing by the sum
    elev.prior = elev.prior / raster::cellStats(elev.prior, "sum", na.rm = TRUE)
    # Mean center prior
    elev.prior = log(elev.prior)
    prior.means = raster::cellStats(elev.prior,"mean")
    raster::values(elev.prior) = do.call('cbind',lapply(1:length(prior.means), function(x) raster::values(elev.prior[[x]]) + abs(prior.means[x])))
    names(elev.prior)='elev.prior'

    # if(getOption("ibis.runparallel")) raster::endCluster()
    # ---- #

    # Check whether an offset exists already
    if(!is.Waiver(x$offset) && add){
      # Add to current object
      of <- x$offset
      elev.prior <- raster::resample(elev.prior, of, method = 'bilinear', func = mean, cl = FALSE)
      names(elev.prior) <- 'elev.prior' # In case the layer name got lost
      of <- raster::stack(of) |> raster::addLayer(elev.prior)
      x <- x$set_offset(of)
    } else {
      # Add as a new offset
      x <- x$set_offset(elev.prior)
    }
    return(x)
  }
)
