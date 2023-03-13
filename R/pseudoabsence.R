#' @include utils.R
NULL

#' Settings for specifying pseudo-absence points within the model background
#'
#' @description
#' This function defines the settings for pseudo-absence sampling of the background. For many
#' engines such points are necessary to model Poisson (or Binomial) distributed point process data.
#' Specifically we call absence points for Binomial (Bernoulli really) distributed responses \code{'pseudo-absence'}
#' and absence data for Poisson responses \code{'background'} points. For more details read Renner et al. (2015).
#'
#' The function \code{'add_pseudoabsence'} allows to add absence points to any [`sf`] object. See **Details** for
#' additional parameter description and examples on how to 'turn' a presence-only dataset into a presence-(pseudo-)absence.
#'
#' @details
#' There are multiple methods available for sampling a biased background layer.
#' Possible parameters for \code{method} are:
#'
#' * \code{'random'} Absence points are generated randomly over the background (Default),
#' * \code{'buffer'} Absence points are generated only within a buffered distance of existing points.
#' This option requires the specification of the parameter \code{buffer_distance}.
#' * \code{'mcp'} Can be used to only generate absence points within or outside a
#' minimum convex polygon of the presence points. The parameter \code{inside} specifies whether
#' points should be sampled inside or outside (Default) the minimum convex polygon.
#' * \code{'range'} Absence points are created either inside or outside a provided additional
#' layer that indicates for example a range of species (controlled through parameter \code{inside}).
#' * \code{'zones'} A ratified (e.g. of type [factor]) [`RasterLayer`] layer depicting zones from which absence
#' points are to be sampled. This method checks which points fall within which zones and then samples absence
#' points either within or outside these zones exclusively. Both \code{'layer'} and \code{'inside'} have to be set
#' for this option.
#' * \code{'target'} Make use of a target background for sampling absence points. Here a [`RasterLayer`] object
#' has to be provided through the parameter \code{'layer'}. Absence points are then sampled exclusively within
#' the target areas for grid cells with non-zero values.
#'
#' @param background A [`RasterLayer`] or [`sf`] object over which background points can be sampled. Default is
#' \code{NULL} (Default) and the background is then added when the sampling is first called.
#' @param nrpoints A [`numeric`] given the number of absence points to be created. Has to be larger than \code{0} and
#' normally points are not created in excess of the number of cells of the background (Default: \code{10 000}).
#' @param method [`character`] denoting how the sampling should be done. See details for options (Default: \code{"random"}).
#' @param buffer_distance [`numeric`] A distance from the observations in which pseudo-absence points are not to be generated.
#' Note that units follow the units of the projection (e.g. \code{m} or \code{°}). Only used when \code{method = "buffer"}.
#' @param inside A [`logical`] value of whether absence points should be sampled outside (Default) or inside a
#' minimum convex polygon or range provided the respective method is chosen (parameter \code{method = "mcp"} or
#' \code{method = "range"}).
#' @param min_ratio A [`numeric`] with the minimum ratio of background points relative to the presence points.
#' Setting this value to \code{1} generates an equal amount of absence points relative to the presence points.
#' Usually ignored unless the ratio exceeds the \code{nrpoints} parameters (Default: \code{0.25}).
#' @param layer A [`sf`] or [`RasterLayer`] (in the case of method \code{'zones'}) object indicating the range of a species.
#' Only used with \code{method = "range"} or \code{method = "zones"} (Default: \code{NULL}).
#' @param bias A [`RasterLayer`] with the same extent and projection and background. Absence points will
#' be preferentially sampled in areas with higher (!) bias. (Default: \code{NULL}).
#' @param ... Any other settings to be added to the pseudoabs settings.
#' @examples
#' \dontrun{
#' # This setting generates 10000 pseudo-absence points outside the minimum convex polygon of presence points
#' ass1 <- pseudoabs_settings(nrpoints = 10000, method = 'mcp', inside = FALSE)
#'
#' # This setting would match the number of presence-absence points directly.
#' ass2 <- pseudoabs_settings(nrpoints = 0, min_ratio = 1)
#'
#' # These settings can then be used to add pseudo-absence data to a presence-only dataset
#' # this effectively adds these simulated absence points to the resulting model
#' all_my_points <- add_pseudoabsence(df = virtual_points, field_occurrence = 'Observed',
#'                                      template = background, settings = ass1)
#' }
#' @references
#' * Renner IW, Elith J, Baddeley A, Fithian W, Hastie T, Phillips SJ, Popovic G, Warton DI. 2015. Point process models for presence-only analysis. Methods in Ecology and Evolution 6:366–379. DOI: 10.1111/2041-210X.12352.
#' * Renner, I. W., & Warton, D. I. (2013). Equivalence of MAXENT and Poisson point process models for species distribution modeling in ecology. Biometrics, 69(1), 274-281.
#' @name pseudoabs_settings
#' @aliases pseudoabs_settings
#' @keywords train
#' @exportMethod pseudoabs_settings
#' @export
NULL
methods::setGeneric("pseudoabs_settings",
                    signature = methods::signature("background"),
                    function(background = NULL, nrpoints = 10000, min_ratio = 0.25,
                             method = "random", buffer_distance = 10000, inside = FALSE,
                             layer = NULL, bias = NULL, ...) standardGeneric("pseudoabs_settings"))

#' @name pseudoabs_settings
#' @rdname pseudoabs_settings
#' @usage \S4method{pseudoabs_settings}{ANY, numeric, numeric, character, numeric, logical, logical ANY}(background, nrpoints, min_ratio, method, buffer_distance, inside, layer, bias)
methods::setMethod(
  "pseudoabs_settings",
  methods::signature(background = "ANY"),
  function(background = NULL, nrpoints = 10000, min_ratio = 0.25,
           method = "random", buffer_distance = 10000, inside = FALSE,
           layer = NULL, bias = NULL, ...){
    # Check inputs
    assertthat::assert_that(
      is.Raster(background) || inherits(background, 'sf') || is.null(background),
      is.numeric(nrpoints),
      is.numeric(min_ratio),
      is.logical(inside),
      (inherits(layer, 'sf') || is.Raster(layer)) || is.null(layer),
      is.character(method),
      is.numeric(buffer_distance),
      is.Raster(bias) || is.null(bias)
    )
    method <- match.arg(method, c("random", "buffer", "mcp", "range", "zones", "target"), several.ok = FALSE)
    # Create the settings object
    settings <- bdproto(NULL, Settings)
    settings$name <- "Background"
    settings$set('background', background)
    # Set all options
    settings$set('nrpoints', nrpoints)
    settings$set('min_ratio', min_ratio)
    settings$set('inside', inside)
    settings$set('method', method)
    settings$set('buffer_distance', buffer_distance)
    settings$set('layer', layer)
    settings$set('bias', bias)
    # Other settings
    mc <- match.call(expand.dots = FALSE)
    settings$data <- c( settings$data, mc$... )

    return(settings)
  }
)

#' Add pseudo-absence points to a point data set
#'
#' @description
#' For most engines, background or pseudo-absence points are necessary. The distinction
#' lies in how the absence data are handled. For [`poisson`] distributed responses,
#' absence points are considered background points over which the intensity of sampling (\code{lambda})
#' is integrated (in a classical Poisson point-process model).
#'
#' In contrast in [`binomial`] distributed responses, the absence information is assumed to
#' be an adequate representation of the true absences and treated by the model as such...
#' Here it is advised to specify absence points in a way that they represent potential true absence,
#' such as for example through targeted background sampling or by sampling them within/outside a given range.
#'
#' @details
#' A [`pseudoabs_settings()`] object can be added to setup how absence points should be sampled.
#' A \code{bias} parameter can be set to specify a bias layer to sample from, for instance a layer
#' of accessibility.
#' Note that when modelling several datasets, it might make sense to check across all datasets
#' whether certain areas are truly absent.
#' By default, the pseudo-absence points are not sampled in areas in which there are already presence points.
#' @note
#' This method removes all columns from the input \code{df} object other than the \code{field_occurrence} column
#' and the coordinate columns (which will be created if not already present).
#' @param df A [`sf`], [`data.frame`] or [`tibble`] object containing point data.
#' @param field_occurrence A [`character`] name of the column containing the presence information (Default: \code{observed}).
#' @param template A [`RasterLayer`] object that is aligned with the predictors (Default: \code{NULL}). If set to \code{NULL},
#' then \code{background} in the [`pseudoabs_settings()`] has to be a [`RasterLayer`] object.
#' @param settings A [`pseudoabs_settings()`] objects. Absence settings are taken from [ibis_options] otherwise (Default).
#' @references
#' * Stolar, J., & Nielsen, S. E. (2015). Accounting for spatially biased sampling effort in presence‐only species distribution modelling. Diversity and Distributions, 21(5), 595-608.
#' * Bird, T.J., Bates, A.E., Lefcheck, J.S., Hill, N.A., Thomson, R.J., Edgar, G.J., Stuart-Smith, R.D., Wotherspoon, S., Krkosek, M., Stuart-Smith, J.F. and Pecl, G.T., 2014. Statistical solutions for error and bias in global citizen science datasets. Biological Conservation, 173, pp.144-154.
#' @keywords train
#' @returns A [`data.frame`] containing the newly created pseudo absence points.
#' @export
add_pseudoabsence <- function(df, field_occurrence = "observed", template = NULL, settings = getOption("ibis.pseudoabsence")){
  assertthat::assert_that(
    is.data.frame(df) || inherits(df, 'sf') || tibble::is_tibble(df),
    is.Raster(template) || is.null(template),
    is.character(field_occurrence),
    assertthat::has_name(df, field_occurrence),
    inherits(settings, "Settings")
  )
  # Check that no 0 are present, otherwise raise a warning.
  assertthat::see_if( any(df[[field_occurrence]] != 0) )

  # Try and guess the geometry
  if(!inherits(df, 'sf')) df <- guess_sf(df)
  assertthat::assert_that(inherits(df, 'sf'), msg = "Could not convert input to sf. Prepare data first.")
  # Add coordinates if not present
  if(!assertthat::has_name(df, 'x') || !assertthat::has_name(df, 'y')) {
    df$x <- sf::st_coordinates(df[attr(df, "sf_column")])[,1]
    df$y <- sf::st_coordinates(df[attr(df, "sf_column")])[,2]
  }
  # Select relevant columns and assign type
  df$type <- "Presence"

  # Check whether the background is set and if not, use the bbox
  if(is.Waiver(settings$get("background"))){
    # Check whether temlate wasn't provided
    if(!is.null(template)){
      background <- template
    } else {
      background <- sf::st_as_sf(
        sf::st_as_sfc(sf::st_bbox(df))
      );background$bg <- 1
    }
  } else {
    background <- settings$get("background")
  }
  # Check that background is a raster, otherwise rasterize with identical resolution
  if(!is.Raster(background)){
    assertthat::assert_that(is.Raster(template),
                            msg = "No suitable RasterLayer was provided through Settings or as template!")
    if("fasterize" %in% utils::installed.packages()[,1]){
      background <- fasterize::fasterize(sf = background, raster = emptyraster(template), field = NULL)
    } else {
      background <- raster::rasterize(background, emptyraster(template), field = 1)
    }
    assertthat::assert_that(is.Raster(background))
  }

  # --- #
  # Now depending on the settings create absence points
  # Get number of points to sample and ratio
  nrpoints <- settings$get('nrpoints')
  min_ratio <- settings$get('min_ratio')

  method <- settings$get('method')
  buffer_distance <- settings$get('buffer_distance')

  # If the nr of points is 0, set it equal to the number of min_ratio or presented presence points
  nrpoints <- max(nrpoints, round( nrow(df) * min_ratio ))
  if(nrpoints > raster::ncell(background)) nrpoints <- raster::ncell(background)

  if(!is.Waiver(settings$get("bias"))){
    bias <- settings$get("bias")
    if(!compareRaster(bias, background, stopiffalse = FALSE)){
      # Resample to ensure same coverage
      bias <- raster::crop(bias, raster::extent(background))
      bias <- raster::resample(bias, background, method = "bilinear")
    }
    # Normalize if not already set
    if(raster::cellStats(bias, 'max') > 1 || raster::cellStats(bias, 'min') < 0 ){
      bias <- predictor_transform(bias, option = "norm")
    }
  } else { bias <- NULL }

  # Rasterize the presence estimates
  bg1 <- raster::rasterize(df[,c('x','y')] %>% sf::st_drop_geometry(),
                           background, fun = 'count', background = 0)
  bg1 <- raster::mask(bg1, background)

  assertthat::assert_that(
    is.finite(raster::cellStats(bg1,'max',na.rm = T)[1])
  )

  # Generate pseudo absence data
  if(method == "random"){
    # Now sample from all cells not occupied
    if(!is.null(bias)){
      # Get probability values for cells where no sampling has been conducted
      prob_bias <- bias[which(bg1[]==0)]
      if(any(is.na(prob_bias))) prob_bias[is.na(prob_bias)] <- 0
      abs <- sample(which(bg1[]==0), size = nrpoints, replace = TRUE, prob = prob_bias)
    } else {
      # raster::sampleStratified(bg1, nrpoints)[,1]
      abs <- sample(which(bg1[]==0), size = nrpoints, replace = TRUE)
    }
  } else if(method == "buffer"){
    assertthat::assert_that(is.numeric(buffer_distance),msg = "Buffer distance parameter not numeric!")
    # Get units of projection and print for
    un <- sf:::crs_parameters(sf::st_crs(df))$ud_unit
    if(getOption('ibis.setupmessages')) myLog('[Export]','yellow', paste0('Calculating pseudo-absence outside a ', buffer_distance ,units::deparse_unit(un),' buffer'))
    # Calculate buffer
    buf <- sf::st_buffer(x = df, dist = buffer_distance)
    if("fasterize" %in% utils::installed.packages()[,1]){
      buf <- fasterize::fasterize(sf = buf, raster = emptyraster(template), field = NULL)
    } else {
      buf <- raster::rasterize(buf, emptyraster(template), field = 1)
    }
    bg2 <- raster::mask(template, buf, inverse = FALSE, updatevalue = 1)
    assertthat::assert_that(cellStats(bg2, "max")>0,msg = "Considered buffer distance too big!")
    # Now sample from all cells not occupied
    if(!is.null(bias)){
      # Get probability values for cells where no sampling has been conducted
      prob_bias <- bias[which(bg2[]==1)]
      if(any(is.na(prob_bias))) prob_bias[is.na(prob_bias)] <- 0
      abs <- sample(which(bg2[]==1), size = nrpoints, replace = TRUE, prob = prob_bias)
    } else {
      # raster::sampleStratified(bg1, nrpoints)[,1]
      abs <- sample(which(bg2[]==1), size = nrpoints, replace = TRUE)
    }
    rm(bg2, buf)
  } else if(method == "mcp"){
    # Idea is to draw a MCP polygon around all points and sample background points only outside of it.
    pol <- sf::st_as_sf( sf::st_convex_hull(sf::st_union(df)) )
    # Now mask out this area from the background
    inside <- settings$get("inside")
    assertthat::assert_that(is.logical(inside), msg = "MCP inside / outside parameter has to be set.")
    bg2 <- raster::mask(bg1, mask = pol, inverse = !inside)
    if(!is.null(bias)){
      # Get probability values for cells where no sampling has been conducted
      prob_bias <- bias[which(bg2[]==0)]
      if(any(is.na(prob_bias))) prob_bias[is.na(prob_bias)] <- 0
      abs <- sample(which(bg2[]==0), size = nrpoints, replace = TRUE, prob = prob_bias)
    } else {
      abs <- sample(which(bg2[]==0), size = nrpoints, replace = TRUE)
    }
    rm(bg2)
  } else if(method == "range"){
    assertthat::assert_that(!is.null(settings$get("layer")),
                            inherits(settings$get("layer"), "sf"),
                            msg = "For method range a layer in sf format has to be specified!")
    # Get the layer
    layer <- settings$get("layer")
    if(sf::st_crs(layer) != sf::st_crs(df)){
      layer <- sf::st_transform(layer, crs = sf::st_crs(df))
    }
    # Get parameter on location
    inside <- settings$get("inside")

    # Now mask out from the sampling background the area from which to sample
    bg2 <- raster::mask(bg1, mask = layer, inverse = !inside)
    if(!is.null(bias)){
      # Get probability values for cells where no sampling has been conducted
      prob_bias <- bias[which(bg2[]==0)]
      if(any(is.na(prob_bias))) prob_bias[is.na(prob_bias)] <- 0
      abs <- sample(which(bg2[]==0), size = nrpoints, replace = TRUE, prob = prob_bias)
    } else {
      abs <- sample(which(bg2[]==0), size = nrpoints, replace = TRUE)
    }
    rm(bg2)
  } else if(method == "zones"){
    assertthat::assert_that(is.Raster(settings$get("layer")),
                            is.factor(settings$get("layer")),
                            msg = "For method zones a factorized RasterLayer has to be specified!")
    # Get the layer
    layer <- settings$get("layer")
    if(!raster::compareCRS(layer, bg1)){
      layer <- raster::projectRaster(layer, bg1,method = "ngb")
    }
    # Get parameter on location
    inside <- settings$get("inside")

    # Cross-tabulate with presence raster
    tab <- raster::crosstab(layer, bg1, long = TRUE)
    tab <- tab[tab[,2]>0,] # Get only those zones where there are presence points
    # Remove any 0 class if there
    if(any(tab[,1] == 0)) tab <- tab[tab[,1]!=0,]

    # Now make layer of the zones that have only the values include or non-include
    if(inside){
      zones <- layer %in% unique(tab[,1])
    } else {
      zones <- !layer %in% unique(tab[,1])
    }
    # Mask again to be sure
    zones <- raster::mask(zones, bg1)
    zones[zones==0] <- NA

    # Now mask out from the sampling background the area from which to sample
    bg2 <- raster::mask(bg1, mask = zones)
    if(!is.null(bias)){
      # Get probability values for cells where no sampling has been conducted
      prob_bias <- bias[which(bg2[]==0)]
      if(any(is.na(prob_bias))) prob_bias[is.na(prob_bias)] <- 0
      abs <- sample(which(bg2[]==0), size = nrpoints, replace = TRUE, prob = prob_bias)
    } else {
      abs <- sample(which(bg2[]==0), size = nrpoints, replace = TRUE)
    }
    rm(bg2)
  } else if(method == "target"){
    assertthat::assert_that(is.Raster(settings$get("layer")),
                            msg = "For method target a rasterized occurrence layer has to be specified!")
    # Get the layer
    layer <- settings$get("layer")
    if(!raster::compareCRS(layer, bg1)){
      layer <- raster::projectRaster(layer, bg1,method = "ngb")
    }
    # Now mask out from the sampling background the area from which to sample
    bg2 <- raster::mask(bg1, mask = layer)
    if(!is.null(bias)){
      # Get probability values for cells where no sampling has been conducted
      prob_bias <- bias[which(bg2[]==0)]
      if(any(is.na(prob_bias))) prob_bias[is.na(prob_bias)] <- 0
      abs <- sample(which(bg2[]==0), size = nrpoints, replace = TRUE, prob = prob_bias)
    } else {
      abs <- sample(which(bg2[]==0), size = nrpoints, replace = TRUE)
    }
    rm(bg2)
  } else {
    stop("Method not found. Check validity of absence settings!")
  }

  # Append to the presence information.
  abs <- sf::st_as_sf(data.frame(raster::xyFromCell(bg1, abs)),
                      coords = c("x", "y"), crs = sf::st_crs(df) )
  abs$x <- sf::st_coordinates(abs)[,1]; abs$y <- sf::st_coordinates(abs)[,2]
  abs$type <- "Pseudo-absence"; abs[[field_occurrence]] <- 0
  sf::st_geometry(abs) <- attr(df, "sf_column") # Rename geom column to be the same as for df
  assertthat::assert_that( nrow(abs) > 0,
                           all(names(abs) %in% names(df)))
  # Unique to remove any duplicate values (otherwise double counted cells)
  # FIXME: Ignoring this as one might want to stress contrast to biases cells
  # abs <- unique(abs)
  # Combine with presence information and return
  out <- rbind.data.frame(subset(df, select = names(abs)),
                          abs)
  return(out)
}
