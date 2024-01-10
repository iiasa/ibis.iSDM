#' Specify a spatial explicit offset
#'
#' @description Including offsets is another option to integrate spatial prior
#' information in linear and additive regression models. Offsets shift the
#' intercept of the regression fit by a certain amount. Although only one offset
#' can be added to a regression model, it is possible to combine several
#' spatial-explicit estimates into one offset by calculating the sum of all
#' spatial-explicit layers.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param layer A [`sf`] or [`SpatRaster`] object with the range for the target
#' feature.
#' @param add [`logical`] specifying whether new offset is to be added. Setting
#' this parameter to \code{FALSE} replaces the current offsets with the new
#' one (Default: \code{TRUE}).
#'
#' @details This function allows to set any specific offset to a regression
#' model. The offset has to be provided as spatial [`SpatRaster`] object. This
#' function simply adds the layer to a [`distribution()`] object.
#' **Note that any transformation of the offset (such as \code{log}) has do be done externally!**
#'
#' If the layer is range and requires additional formatting, consider using the
#' function [`add_offset_range()`] which has additional functionalities such
#' such distance transformations.
#'
#' @note Since offsets only make sense for linear regressions (and not for
#' instance regression tree based methods such as [engine_bart]), they do not
#' work for all engines. Offsets specified for non-supported engines are ignored
#' during the estimation
#'
#' @returns Adds an offset to a [`distribution`] object.
#'
#' @references
#' * Merow, C., Allen, J.M., Aiello-Lammens, M., Silander, J.A., 2016. Improving niche and
#' range estimates with Maxent and point process models by integrating spatially explicit
#' information. Glob. Ecol. Biogeogr. 25, 1022–1036. https://doi.org/10.1111/geb.12453
#'
#' @family offset
#' @keywords prior offset
#'
#' @examples
#' \dontrun{
#'  x <- distribution(background) |>
#'    add_predictors(covariates) |>
#'    add_offset(nicheEstimate)
#' }
#'
#' @name add_offset
NULL

#' @rdname add_offset
#' @export
methods::setGeneric(
  "add_offset",
  signature = methods::signature("x", "layer"),
  function(x, layer, add = TRUE) standardGeneric("add_offset"))

#' @rdname add_offset
methods::setMethod(
  "add_offset",
  methods::signature(x = "BiodiversityDistribution", layer = "SpatRaster"),
  function(x, layer, add = TRUE) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.Raster(layer),
                            is.logical(add)
    )
    # Messenger
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding spatial explicit offset...')

    # Sanitize names if specified
    if(getOption('ibis.cleannames')) names(layer) <- sanitize_names(names(layer))
    ori.name <- names(layer)

    # Check for infinite values
    assertthat::assert_that(
      all( is.finite( terra::global(layer, "range", na.rm = TRUE)[,1]) ),
      msg = "Infinite values found in the layer (maybe log of 0?)."
    )

    # Check that background and range align, otherwise raise error
    if(is_comparable_raster(layer, x$background)){
      layer <- alignRasters(layer, x$background, method = 'bilinear', func = mean, cl = FALSE)
      names(layer) <- ori.name
    }

    # Check whether an offset exists already
    if(!is.Waiver(x$offset) && add){
      # Add to current object
      of <- x$offset
      layer <- terra::resample(layer, of, method = 'bilinear', threads = getOption("ibis.nthread"))
      names(layer) <- ori.name # In case the layer name got lost
      of <- c(of, layer)
      x <- x$set_offset(of)
    } else {
      # Add as a new offset
      x <- x$set_offset(layer)
    }
    return(x)
  }
)

#' @rdname add_offset
methods::setMethod(
  "add_offset",
  methods::signature(x = "BiodiversityDistribution", layer = "sf"),
  function(x, layer, add = TRUE) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            inherits(layer, "sf"),
                            is.logical(add)
    )

    # Template raster for rasterization background
    if(!is.Waiver(x$predictors)){
      temp <- emptyraster(x$predictors$get_data())
    } else {
      # Try and guess an sensible background raster
      myLog('[Setup]','red',
            'CAREFUL - This might not work without predictors already in the model.
      Add offset after predictors')
      temp <- terra::rast( extent = terra::ext(x$background),
                           resolution = diff(sf::st_bbox(x$background)[c(1,3)]) / 100,
                           crs = terra::crs(x$background))
    }

    # Check to make the entries valid
    if( any(!sf::st_is_valid(layer)) ){
      layer <- sf::st_make_valid(layer) # Check whether to make them valid
      if( any(!sf::st_is_valid(layer)) ){
        # If still has errors, combine
        suppressMessages( layer <- layer |> sf::st_combine() |> sf::st_as_sf() )
      }
    }

    # If layer has multiple entries join them
    if(nrow(layer)>1) suppressMessages( layer <- layer |> sf::st_union() |> sf::st_as_sf() )

    # Rasterize the range
    ras_range <- terra::rasterize(layer, temp, field = 1, background = 0)
    ras_range <- terra::mask(ras_range, x$background)
    names(ras_range) <-  "spatial_offset"

    # Call with new SpatRaster object
    x <- add_offset(x, ras_range, add)
    return(x)
  }
)

#' Function to remove an offset
#'
#' @description This is just a wrapper function for removing specified offsets
#' from a [`BiodiversityDistribution-class`]) object.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param layer A `character` pointing to the specific layer to be removed. If
#'   set to \code{NULL}, then all offsets are removed from the object.
#'
#' @returns Removes an offset from a [`distribution`] object.
#'
#' @family offset
#' @keywords prior offset
#'
#' @examples
#' \dontrun{
#'  rm_offset(model) -> model
#' }
#'
#' @name rm_offset
NULL

#' @rdname rm_offset
#' @export
methods::setGeneric(
  "rm_offset",
  signature = methods::signature("x"),
  function(x, layer = NULL) standardGeneric("rm_offset"))

#' @rdname rm_offset
methods::setMethod(
  "rm_offset",
  methods::signature(x = "BiodiversityDistribution"),
  function(x, layer = NULL) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            missing(layer) || is.character(layer) || is.null(layer)
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
#' @description Including offsets is another option to integrate spatial prior
#' information in linear and additive regression models. Offsets shift the
#' intercept of the regression fit by a certain amount. Although only one offset
#' can be added to a regression model, it is possible to combine several
#' spatial-explicit estimates into one offset by calculating the sum of all
#' spatial-explicit layers.
#'
#' @inheritParams add_offset
#' @param points An optional [`sf`] object with key points. The location of the
#' points are then used to calculate the probability that a cell has been
#' sampled while accounting for area differences. (Default: \code{NULL}).
#'
#' @details This functions emulates the use of the [`add_offset()`] function,
#' however applies an inverse transformation to remove the provided layer from
#' the overall offset. So if for instance a offset is already specified (such as
#' area), this function removes the provided \code{bias.layer} from it via
#' \code{"offset(log(off.area)-log(bias.layer))"}
#'
#' **Note that any transformation of the offset (such as \code{log}) has do be done externally!**
#'
#' If a generic offset is added, consider using the [`add_offset()`] function.
#' If the layer is a expert-based range and requires additional parametrization,
#' consider using the function [`add_offset_range()`] or the \code{bossMaps}
#' R-package.
#'
#' @returns Adds a bias offset to a [`distribution`] object.
#'
#' @references
#' * Merow, C., Allen, J.M., Aiello-Lammens, M., Silander, J.A., 2016. Improving
#' niche and range estimates with Maxent and point process models by integrating
#' spatially explicit information. Glob. Ecol. Biogeogr. 25, 1022–1036. https://doi.org/10.1111/geb.12453
#'
#' @family offset
#' @keywords prior offset
#'
#' @examples
#' \dontrun{
#'  x <- distribution(background) |>
#'    add_predictors(covariates) |>
#'    add_offset_bias(samplingBias)
#' }
#'
#' @name add_offset_bias
NULL

#' @rdname add_offset_bias
#' @export
methods::setGeneric(
  "add_offset_bias",
  signature = methods::signature("x", "layer"),
  function(x, layer, add = TRUE, points = NULL) standardGeneric("add_offset_bias"))

#' @rdname add_offset_bias
methods::setMethod(
  "add_offset_bias",
  methods::signature(x = "BiodiversityDistribution", layer = "SpatRaster"),
  function(x, layer, add = TRUE, points = NULL) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.Raster(layer),
                            is.logical(add),
                            is.null(points) || inherits(points, 'sf')
    )
    # Messenger
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding spatial explicit bias offset...')

    # Sanitize names if specified
    if(getOption('ibis.cleannames')) names(layer) <- sanitize_names(names(layer))
    ori.name <- names(layer)

    # Check that background and range align, otherwise raise error
    if(is_comparable_raster(layer, x$background)){
      layer <- alignRasters(layer, x$background, method = 'bilinear', func = mean, cl = FALSE)
      names(layer) <- ori.name
    }
    if(is.null(points)){
      # Since it is a bias offset and removal is equivalent to simple subtraction, multiply with *-1
      layer <- layer * -1
    } else {
      ## Count the number of records per cell
      tab <- terra::cellFromXY(layer, sf::st_coordinates(points))
      r <- emptyraster(layer)
      r[tab] <- layer[tab]
      r <- terra::mask(r, background)

      ## Make zeros a very small number otherwise issues with log(0).
      r[r[]==0] <- 1e-6
      suppressWarnings({ar <- terra::cellSize(r)})

      ## Calculate the probability that a cell has been sampled while accounting
      ## for area differences in lat/lon Direction sign is negative and if area
      ## offset considered, use "+ offset(log(off.area)-log(off.bias))"
      off.bias <- (-log(1-exp(-r * ar)) - log(ar))
      names(off.bias) <- "off.bias"
      # Add bias as covariate
      layer <- off.bias
    }

    # Check for infinite values
    assertthat::assert_that(
      all( is.finite( terra::global(layer, "range", na.rm = TRUE)[,1]) ),
      msg = "Infinite values found in the layer (maybe log of 0?)."
    )

    # Check whether an offset exists already
    if(!is.Waiver(x$offset) && add){
      # Add to current object
      of <- x$offset
      layer <- terra::resample(layer, of, method = 'bilinear', threads = getOption("ibis.nthread"))
      names(layer) <- ori.name # In case the layer name got lost
      suppressWarnings( of <- c( of, layer ) )
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
#' @description This function has additional options compared to the more
#' generic [`add_offset()`], allowing customized options specifically for
#' expert-based ranges as offsets or spatialized polygon information on species
#' occurrences. If even more control is needed, the user is informed of the
#' \code{"bossMaps"} package Merow et al. (2017). Some functionalities of that
#' package emulated through the \code{"distance_function"} set to \code{"log"}.
#' This tries to fit a 5-parameter logistic function to estimate the distance
#' from the range (Merow et al. 2017).
#'
#' @inheritParams add_offset
#' @param distance_max A [`numeric`] threshold on the maximum distance beyond
#' the range that should be considered to have a high likelihood of containing
#' species occurrences (Default: \code{Inf} \code{"m"}). Can be set to
#' \code{NULL} or \code{0} to indicate that no distance should be calculated.
#' @param family A [`character`] denoting the type of model to which this offset
#' is to be added. By default it assumes a \code{'poisson'} distributed model
#' and as a result the output created by this function will be
#' log-transformed. If however a \code{'binomial'} distribution is chosen,
#' than the output will be \code{`logit`} transformed. For integrated models
#' leave at default.
#' @param presence_prop [`numeric`] giving the proportion of all records
#' expected to be inside the range. By default this is set to \code{0.9}
#' indicating that 10% of all records are likely outside the range.
#' @param distance_clip [`logical`] as to whether distance should be clipped
#' after the maximum distance (Default: \code{FALSE}).
#' @param distance_function A [`character`] specifying the distance function to
#' be used. Available are linear (\code{"linear"}), negative exponential kernels (\code{"negexp"},
#' default) and a five parameters logistic curve (code{"logcurve"}) as
#' proposed by Merow et al. 2017.
#' @param point An optional [`sf`] layer with points or [`logical`] argument. In
#' the case of the latter the point data is ignored (Default: \code{FALSE}).
#' @param field_occurrence A [`numeric`] or [`character`] location of
#' biodiversity point records.
#' @param fraction An optional [`SpatRaster`] object that is multiplied with
#' digitized raster layer. Can be used to for example to remove or reduce the
#' expected value (Default: \code{NULL}).
#'
#' @details The output created by this function creates a [`SpatRaster`] to be
#' added to a provided distribution object. Offsets in regression models are
#' likelihood specific as they are added directly to the overall estimate of
#' \code{`y^hat`}.
#'
#' Note that all offsets created by this function are by default log-transformed
#' before export. Background values (e.g. beyond \code{"distance_max"}) are set
#' to a very small constant (\code{1e-10}).
#'
#' @returns Adds a range offset to a [`distribution`] object.
#'
#' @references
#' * Merow, C., Wilson, A.M., Jetz, W., 2017. Integrating occurrence data and expert
#' maps for improved species range predictions. Glob. Ecol. Biogeogr. 26, 243–258.
#' https://doi.org/10.1111/geb.12539
#' * Merow, C., Allen, J.M., Aiello-Lammens, M., Silander, J.A., 2016. Improving
#' niche and range estimates with Maxent and point process models by integrating
#' spatially explicit information. Glob. Ecol. Biogeogr. 25, 1022–1036.
#' https://doi.org/10.1111/geb.12453
#'
#' @seealso \code{"bossMaps"}
#' @family offset
#' @keywords prior offset
#'
#' @examples
#' \dontrun{
#'  # Train a presence-only model with a simple offset
#'  fit <- distribution(background) |>
#'  add_biodiversity_poipo(virtual_points, field_occurrence = "Observed") |>
#'  add_predictors(predictors) |>
#'  add_offset_range(virtual_range, distance_max = 5,distance_function = "logcurve",
#'  distance_clip = TRUE ) |>
#'  engine_glm() |>
#'  train()
#' }
#'
#' @name add_offset_range
NULL

#' @rdname add_offset_range
#' @export
methods::setGeneric(
  "add_offset_range",
  signature = methods::signature("x", "layer"),
  function(x, layer, distance_max = Inf, family = "poisson", presence_prop = 0.9,
           distance_clip = FALSE, distance_function = "negexp",
           field_occurrence = "observed", fraction = NULL,
           point = FALSE, add = TRUE) standardGeneric("add_offset_range"))

#' Function for when raster is directly supplied (precomputed)
#' @rdname add_offset_range
methods::setMethod(
  "add_offset_range",
  methods::signature(x = "BiodiversityDistribution", layer = "SpatRaster"),
  function(x, layer, fraction = NULL, add = TRUE) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.Raster(layer),
                            is.null(fraction) || is.Raster(fraction),
                            is.logical(add)
    )
    # Messenger
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding range offset...')

    # Sanitize names if specified
    if(getOption('ibis.cleannames')) names(layer) <- sanitize_names(names(layer))
    ori.name <- names(layer)

    # Check for infinite values
    assertthat::assert_that(
      all( is.finite(terra::global(layer, "range", na.rm = TRUE)) ),
      msg = "Infinite values found in the layer (maybe log of 0?)."
    )

    # Check that background and range align, otherwise raise error
    if(is_comparable_raster(layer, x$background)){
      warning('Supplied range does not align with background! Aligning them now...')
      layer <- alignRasters(layer, x$background, method = 'bilinear', func = mean, cl = FALSE)
      names(layer) <- ori.name # In case the layer name got lost
    }

    # Multiply with fraction layer if set
    if(!is.null(fraction)){
      # Rescale if necessary and set 0 to a small constant 1e-6
      if(terra::global(fraction, "min")[,1] < 0) fraction <- predictor_transform(fraction, option = "norm")
      fraction[fraction==0] <- 1e-6
      layer <- layer * fraction
    }

    # Check whether an offset exists already
    if(!is.Waiver(x$offset) && add){
      # Add to current object
      of <- x$offset
      layer <- terra::resample(layer, of, method = 'bilinear', threads = getOption("ibis.nthread"))
      names(layer) <- ori.name # In case the layer name got lost
      suppressWarnings( of <- c( of, layer ) )
      x <- x$set_offset(of)
    } else {
      # Add as a new offset
      x <- x$set_offset(layer)
    }
    return(x)
  }
)

#' @rdname add_offset_range
methods::setMethod(
  "add_offset_range",
  methods::signature(x = "BiodiversityDistribution", layer = "sf"),
  function(x, layer, distance_max = Inf, family = "poisson", presence_prop = 0.9,
           distance_clip = FALSE, distance_function = "negexp",
           field_occurrence = "observed", fraction = NULL, point = FALSE, add = TRUE ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            inherits(layer, 'sf'),
                            is.null(distance_max) || is.numeric(distance_max) || is.infinite(distance_max),
                            is.numeric(presence_prop),
                            is.logical(distance_clip),
                            is.character(distance_function),
                            is.null(fraction) || is.Raster(fraction),
                            is.character(family),
                            inherits(point, "sf") || is.logical(point),
                            is.character(field_occurrence),
                            is.logical(add)
    )
    # distance_max = Inf; family = "poisson"; presence_prop = 0.9; distance_clip = FALSE; distance_function = "negexp"; field_occurrence = "observed"; fraction = NULL; add = TRUE;point =NULL
    # Match the type if set
    family <- match.arg(family, c("poisson", "binomial"), several.ok = FALSE)

    # Distance function
    distance_function <- match.arg(distance_function, c("linear","negexp", "logcurve"), several.ok = FALSE)

    # Check that necessary dependency is present for log curve
    if(distance_function=="logcurve"){
      check_package("gnlm")
      if(!("gnlm" %in% loadedNamespaces()) || ('gnlm' %notin% utils::sessionInfo()$otherPkgs) ) {
        try({requireNamespace('gnlm');attachNamespace("gnlm")},silent = TRUE)
      }
    }

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
      temp <- terra::rast( extent = terra::ext(x$background),
                             resolution = diff(sf::st_bbox(x$background)[c(1,3)]) / 100,
                             crs = terra::crs(x$background))
    }

    # Check to make the entries valid
    if( any(!sf::st_is_valid(layer)) ){
      layer <- sf::st_make_valid(layer) # Check whether to make them valid
      if( any(!sf::st_is_valid(layer)) ){
        # If still has errors, combine
        suppressMessages( layer <- layer |> sf::st_combine() |> sf::st_as_sf() )
      }
    }

    # If layer has multiple entries join them
    if(nrow(layer)>1) suppressMessages( layer <- layer |> sf::st_union() |> sf::st_as_sf() )

    # Get Point information if set
    if(isTRUE(point) || is.null(point)){
      assertthat::assert_that(length(x$get_biodiversity_types()) > 0)
      #TODO: Collate point from x
      stop("Automatic point collation not yet implemented. Please supply a sf layer to point!")
    } else if(inherits(point, 'sf')){
      assertthat::assert_that( assertthat::has_name(point, field_occurrence),
                               nrow(point)>1)
      # Transform to be sure
      point <- point |> sf::st_transform(crs = sf::st_crs(layer))
      # If family is poisson distributed, add some pseudo-absence points
      if(family=="poisson"){
        point <- add_pseudoabsence(point, field_occurrence = field_occurrence,
                                   template = terra::init(temp,1),
                                   settings = pseudoabs_settings(nrpoints = 0,
                                                                 min_ratio = 1,
                                                                 layer = layer,
                                                                 method = "range",
                                                                 inside = FALSE))
      }
    }

    # Rasterize the range
    ras_range <- terra::rasterize(layer, temp, field = 1, background = 0)
    ras_range <- terra::mask(ras_range, x$background)

    # Calculate distance if required
    if(distance_max > 0){
      # Calculate a distance raster in km
      dis <- terra::gridDist(ras_range, target = 1, scale = 1000)
      # If max distance is specified
      if(distance_clip && is.finite(distance_max)){
        dis[dis > distance_max] <- NA # Set values above threshold to a very small constant
      }
      # Inverse of distance
      if(is.infinite(distance_max)) distance_max <- terra::global(dis, "max", na.rm = TRUE)[,1]
      suppressWarnings( ar <- terra::cellSize(ras_range, unit = "km") ) # Calculate area in km
      # ---- #
      if(distance_function == "negexp"){
        alpha <- 1 / (distance_max / 4 ) # Divide by 4 for a quarter in each direction
        # Grow baseline raster by using an exponentially weighted kernel
        dis <- terra::app(dis, fun = function(x) exp(-alpha * x))
        # Set the remaining ones to very small constant
        dis[is.na(dis)] <- 1e-10 # Background values
        dis <- terra::mask(dis, x$background)

        # Inside I want all X across the entire area for the PPMs, indicating a
        # lambda per area of at least X/A (per unit area) within the range
        pres <- 1 + ( ( terra::global(ar * ras_range, "sum", na.rm = TRUE)[,1] / terra::global(ar, "sum", na.rm = TRUE)[,1]) * (presence_prop) )
        abs <- 1 + ( ( terra::global(ar * ras_range, "sum", na.rm = TRUE)[,1] / terra::global(ar, "sum", na.rm = TRUE)[,1]) * (1-presence_prop) )
        # Now set all values inside the range to pres and outside to abs
        ras_range[ras_range == 1] <- pres
        ras_range[ras_range == 0] <- abs
        # Multiply with distance layer
        ras_range <- ras_range * dis
        # Normalize the result by dividing by the sum
        ras_range <- ras_range / terra::global(ras_range, "sum", na.rm = TRUE)[,1]
      } else if(distance_function == "logcurve"){
        # Extract the point values from the raster
        ex <- get_rastervalue(coords = point, env = dis)
        obs <- point[[field_occurrence]]
        ex <- ex[,names(dis)]
        # Get only valid observation
        if(any(is.na(ex))){
          obs <- obs[which(is.finite(ex))]
          ex <- ex[which(is.finite(ex))]
        }

        if(family == "binomial"){
          assertthat::assert_that( length(unique(obs)) == 2)
          y <- cbind(obs, 1-obs)
        } else y <- cbind(obs)

        # Grid search for optimal parameters
        co <- .searchLogisticCurve(y = y, x = ex,
                                    family = family,
                                    search = TRUE)

        # Convert output to SpatRaster using logistic Richard curve
        ras_range <- logisticRichard(x = dis,
                               upper = co["upper"],
                               lower = co["lower"],
                               rate = co["rate"],
                               shift = co["shift"],
                               skew = co["skew"])
        attr(ras_range, "logistic_coefficients") <- co

      } else if (distance_function == "linear") {
        # Multiply with distance layer
        ras_range <-  abs( dis / terra::global(ras_range, "sum", na.rm = TRUE)[,1]) * -1
        ras_range[is.na(ras_range)] <- terra::global(ras_range, "min", na.rm = TRUE)[,1]
      } else {
        stop("Distance method not yet implemented.")
      }

    } else {
      dis <- ras_range
      dis[is.na(dis)] <- 1e-10 # Background values
      dis <- terra::mask(dis, x$background)
    }

    # Multiply with fraction layer if set
    if(!is.null(fraction)){
      # Rescale if necessary and set 0 to a small constant 1e-6
      if(terra::global(fraction, "min", na.rm = TRUE)[,1] < 0) fraction <- predictor_transform(fraction, option = "norm")
      fraction[fraction==0] <- 1e-6
      ras_range <- ras_range * fraction
    }

    # -------------- #
    # Log transform for better scaling
    if(family %in% c("negexp", "linear")){
      ras_range <- switch (family,
                           "poisson" = terra::app(ras_range, log),
                           "binomial" = terra::app(ras_range, logistic)
      )
    }
    # Rescaling does not affect relative differences.
    ras_range <- terra::scale(ras_range, scale = F)
    names(ras_range) <- "range_distance"

    assertthat::assert_that(
      is.finite( terra::global(ras_range, "max", na.rm = TRUE)[,1] ),
      msg = "Range offset has infinite values. Check parameters!"
    )

    # Sanitize names if specified
    if(getOption('ibis.cleannames')) names(layer) <- sanitize_names(names(layer))

    # Set some attributes
    attr(ras_range, "distance_function") <- distance_function
    attr(ras_range, "distance_max") <- distance_max

    ras_range <- terra::mask(ras_range, x$background)

    # Check whether an offset exists already
    if(!is.Waiver(x$offset) && add){
      # Add to current object
      of <- x$offset
      ori.name <- names(ras_range)
      ras_range <- terra::resample(ras_range, of, method = 'bilinear', threads = getOption("ibis.nthread") )
      names(ras_range) <- ori.name # In case the layer name got lost
      suppressWarnings( of <- c(of, ras_range) )
      x <- x$set_offset(of)
    } else {
      # Add as a new offset
      x <- x$set_offset(ras_range)
    }
    return(x)
  }
)

#' Function to calculate the best logistic Richard's curve given a distance
#' @description Internal function not to be used outside `add_offset_range`.
#'
#' @param y A [`numeric`] of the response.
#' @param x A [`numeric`] of the variable.
#' @param family A [`character`] with the family. Default is `binomial`.
#' @param search A [`logical`] whether grid search by AIC be conducted (Default: \code{TRUE}).
#' @param iniParam The initial parameters for the logistic curve.
#'
#' @returns A [`numeric`] vector with the coefficients of regression.
#'
#' @noRd
#'
#' @keywords internal
.searchLogisticCurve <- function(y, x, family, search = TRUE,
                                 iniParam = c(upper = 1,
                                              lower = 0,
                                              rate = 0.04,
                                              shift = 1,
                                              skew = 0.2)){
  assertthat::assert_that(
    is.numeric(y),
    is.numeric(x),
    length(y)>2, length(x)>2,
    is.character(family),
    is.logical(search),
    is.vector(iniParam) && length(iniParam)==5
  )
  # Check package
  check_package("gnlm")

  family <- match.arg(family, c("binomial", "poisson"), several.ok = FALSE)

  if(family == "binomial"){
    if(is.vector(y)) y <- cbind(y, 1 - y)
  } else if(family == "poisson"){
    family <- "Poisson"
    y <- y[,1]
  }

  # Define search grid parameters
  if(search){
    pp <- expand.grid(upper = 1,
                      lower = seq(0, .75, .15),
                      rate = seq(0.01, 0.3, 0.03),
                      shift = 1,
                      skew = seq(0.1, 0.3, 0.05) )
  } else {pp <- rbind(data.frame(t(iniParam))) }
  assertthat::assert_that(ncol(pp)==5)

  # Now find the best parameter combinations for the given set
  result <- data.frame(i = 1:nrow(pp))
  result$aic <- NA

  # Progress
  if(getOption('ibis.setupmessages')) pb <- progress::progress_bar$new(total = nrow(pp))
  for(i in 1:nrow(pp)){
    if(getOption('ibis.setupmessages')) pb$tick()
    # Default starting
    # holdEnv <- list(y = y,
    #                 x = x,
    #                 iniParam = pp[i,])
    # suppressMessages( attach(holdEnv) )
    # on.exit(detach("holdEnv"))
    xx <- pp[i,] |> as.list()

    # Fit Model
    if(family == "binomial"){
      suppressMessages(
        suppressWarnings(
          logisticParam <- try({
            gnlm::bnlr(y = y, link = "logit",
                                        mu = ~ (upper - ((upper - lower)/(1 + exp(-rate * ( - shift)))^(1/skew))),
                                        pmu = xx)
          }, silent = TRUE)
        )
      )
    } else {
      # Add estimate here
      xx$x <- x
      suppressMessages(
        suppressWarnings(
          logisticParam <- try({
            gnlm::gnlr(y = y, distribution = family,
                       mu = ~ (upper - ((upper - lower)/(1 + exp(-rate * (x - shift)))^(1/skew))),
                       pmu = xx)
          }, silent = TRUE)
        )
      )
    }
    if(!inherits(logisticParam, "try-error")){
      result[i,"aic"] <- logisticParam$aic
    }
    rm(logisticParam)
  }

  # Now get the best combination and refit
  # holdEnv <- list(y = y,
  #                 x = x,
  #                 iniParam = pp[which.min(result$aic),])
  # suppressMessages( attach(holdEnv) )
  # on.exit(detach("holdEnv"))
  xx <- pp[which.min(result$aic),] |> as.list()
  # Fit Model
  if(family == "binomial"){
    suppressMessages(
      suppressWarnings(
        logisticParam <- try({
          gnlm::bnlr(y = y, link = "logit",
                     mu = ~ (upper - ((upper - lower)/(1 + exp(-rate * ( - shift)))^(1/skew))),
                     pmu = xx)
        }, silent = TRUE)
      )
    )
  } else {
    xx$x <- x
    suppressMessages(
      suppressWarnings(
        logisticParam <- try({
          gnlm::gnlr(y = y, distribution = family,
                     mu = ~ (upper - ((upper - lower)/(1 + exp(-rate * (x - shift)))^(1/skew))),
                     pmu = xx)
        }, silent = TRUE)
      )
    )
  }

  # Get the coefficients of the best model
  if(inherits(logisticParam, "try-error")) stop("Offset calculating failed...")
  co <- logisticParam$coefficients
  names(co) <- c("upper", "lower", "rate", "shift", "skew")
  return(co)
}

#### Elevational offset ####

#' Specify elevational preferences as offset
#'
#' @description This function implements the elevation preferences offset
#' defined in Ellis‐Soto et al. (2021). The code here was adapted from the
#' Supporting materials script.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param elev A [`SpatRaster`] with the elevation for a given background.
#' @param pref A [`numeric`] vector of length \code{2} giving the lower and
#' upper bound of known elevational preferences. Can be set to \code{Inf} if
#' unknown.
#' @param rate A [`numeric`] for the rate used in the offset (Default:
#' \code{.0089}). This parameter specifies the decay to near zero probability
#' at elevation above and below the expert limits.
#' @param add [`logical`] specifying whether new offset is to be added. Setting
#' this parameter to \code{FALSE} replaces the current offsets with the new
#' one (Default: \code{TRUE}).
#'
#' @details Specifically this functions calculates a continuous decay and
#' decreasing probability of a species to occur from elevation limits. It
#' requires a [`SpatRaster`] with elevation information. A generalized logistic
#' transform (aka Richard's curve) is used to calculate decay from the suitable
#' elevational areas, with the \code{"rate"} parameter allowing to vary the
#' steepness of decline.
#'
#' Note that all offsets created by this function are by default log-transformed
#' before export. In addition this function also mean-centers the output as
#' recommended by Ellis-Soto et al.
#'
#' @returns Adds a elevational offset to a [`distribution`] object.
#'
#' @references
#' * Ellis‐Soto, D., Merow, C., Amatulli, G., Parra, J.L., Jetz, W., 2021. Continental‐scale
#' 1 km hummingbird diversity derived from fusing point records with lateral and
#' elevational expert information. Ecography (Cop.). 44, 640–652. https://doi.org/10.1111/ecog.05119
#' * Merow, C., Allen, J.M., Aiello-Lammens, M., Silander, J.A., 2016. Improving
#' niche and range estimates with Maxent and point process models by integrating
#' spatially explicit information. Glob. Ecol. Biogeogr. 25, 1022–1036.
#' https://doi.org/10.1111/geb.12453
#'
#' @keywords prior offset
#' @family offset
#'
#' @examples
#' \dontrun{
#'  # Adds the offset to a distribution object
#'  distribution(background) |> add_offset_elevation(dem, pref = c(400, 1200))
#' }
#'
#' @name add_offset_elevation
NULL

#' @rdname add_offset_elevation
#' @export
methods::setGeneric(
  "add_offset_elevation",
  signature = methods::signature("x", "elev", "pref"),
  function(x, elev, pref, rate = .0089, add = TRUE) standardGeneric("add_offset_elevation"))

#' @rdname add_offset_elevation
methods::setMethod(
  "add_offset_elevation",
  methods::signature(x = "BiodiversityDistribution", elev = "SpatRaster", pref = "numeric"),
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
      all( is.finite( terra::global(elev, "range", na.rm = TRUE)[,1]) ),
      msg = "Infinite values found in the layer (maybe log of 0?)."
    )

    # Check that background and range align, otherwise raise error
    if(is_comparable_raster(elev, x$background)){
      warning('Supplied range does not align with background! Aligning them now...')
      elev <- alignRasters(elev, x$background, method = 'bilinear', func = mean, cl = FALSE)
    }

    # ---- #
    # if(getOption("ibis.runparallel")) raster::beginCluster(n =
    # getOption("ibis.nthread")) Now calculate the elevation offset by
    # projecting the values onto the elevation layer max avail > min expert
    tmp.elev1 <- -1 * (elev - pref[1])
    tmp.elev1.1 <- terra::app(tmp.elev1, function(x) logisticRichard(x, 1,100, rate, .2 ))
    # min avail < max expert
    tmp.elev2 <- elev - pref[2]
    tmp.elev2.1 <- terra::app(tmp.elev2, function(x) logisticRichard(x, 1,100, rate, .2))
    # Combine both and calculate the minimum
    elev.prior <- min( c(tmp.elev1.1, tmp.elev2.1))
    rm(tmp.elev1,tmp.elev1.1,tmp.elev2,tmp.elev2.1) # clean up

    # Normalize the result by dividing by the sum
    elev.prior <- elev.prior / terra::global(elev.prior, "sum", na.rm = TRUE)[,1]
    # Mean center prior
    elev.prior <- log(elev.prior)
    prior.means <- terra::global(elev.prior, "mean", na.rm = TRUE)[,1]
    terra::values(elev.prior) <- do.call('cbind', lapply(1:length(prior.means), function(x) terra::values(elev.prior[[x]]) + abs(prior.means[x])) )
    names(elev.prior) <- 'elev.prior'

    # if(getOption("ibis.runparallel")) raster::endCluster()
    # ---- #

    # Sanitize names if specified
    if(getOption('ibis.cleannames')) names(elev.prior) <- sanitize_names(names(elev.prior))

    # Check whether an offset exists already
    if(!is.Waiver(x$offset) && add){
      # Add to current object
      of <- x$offset
      elev.prior <- terra::resample(elev.prior, of, method = 'bilinear', threads = getOption("ibis.nthread"))
      names(elev.prior) <- 'elev.prior' # In case the layer name got lost
      suppressWarnings( of <- c( of, elev.prior ) )
      x <- x$set_offset(of)
    } else {
      # Add as a new offset
      x <- x$set_offset(elev.prior)
    }
    return(x)
  }
)
