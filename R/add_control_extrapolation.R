#' Add a control to a BiodiversityModel object to control extrapolation
#'
#' @description One of the main aims of species distribution models (SDMs) is to project
#' in space and time. For projections a common issue is extrapolation as - unconstrained -
#' SDMs can indicate areas as suitable which are unlikely to be occupied by species
#' or habitats (often due to historic or biotic factors). To some extent this can
#' be related to an insufficient quantification of the niche (e.g. niche truncation
#' by considering only a subset of observations within the actual distribution),
#' in other cases there can also be general barriers or constraints that limit
#' any projections (e.g. islands). This control method adds some of those options
#' to a model distribution object. Currently supported methods are:
#'
#' [*] \code{"zones"} - This is a wrapper to allow the addition of zones to a
#' distribution model object, similar to what is also possible via [distribution()].
#' Required is a spatial layer that describes a environmental zoning.
#'
#' [*] \code{"mcp"} - Rather than using an external or additional layer, this option constraints
#' predictions by a certain distance of points in its vicinity. Buffer distances
#' have to be in the unit of the projection used and can be configured via
#' \code{"mcp_buffer"}.
#'
#' [*] \code{"nt2"} - Constraints the predictions using the multivariate combination novelty index (NT2)
#' following Mesgaran et al. (2014). This method is also available in the [similarity()]
#' function.
#'
#' [*] \code{"shape"} - This is an implementation of the 'shape' method introduced
#' by Velazco et al. (2023). Through a user defined threshold it effectively limits
#' model extrapolation so that no projections are made beyond the extent judged as
#' defensible and informed by the training observations.
#'
#' See also details for further explanations.
#'
#' @details
#' For method \code{"zones"} a zoning layer can be supplied which is then used to intersect
#' the provided training points with. Any projections made with the model can
#' then be constrained so as to not project into areas that do not consider any
#' training points and are unlikely to have any. Examples for zones are for the
#' separation of islands and mainlands, biomes, or lithological soil conditions.
#'
#' If no layer is available, it is also possible to constraint predictions by the
#' distance to a minimum convex polygon surrounding the training points with
#' method \code{"mcp"} (optionally buffered). This can make sense particular for
#' rare species or those fully sampled across their niche.
#'
#' For the \code{"NT2"} and \code{"MESS"} index it is possible to constrain
#' the prediction to conditions within (\code{novel = "within"}) or also include
#' outside (\code{novel = "outside"}) conditions.
#'
#' @note
#' The method \code{"zones"} is also possible directly within [distribution()].
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param layer A [`SpatRaster`] or [`sf`] object that limits the prediction
#'   surface when intersected with input data (Default: \code{NULL}).
#' @param method A [`character`] vector describing the method used for controlling
#'  extrapolation. Available options are \code{"zones"}, \code{"mcp"} (Default),
#'  or \code{"nt2"} or \code{"shape"}.
#' @param mcp_buffer A [`numeric`] distance to buffer the mcp (Default
#'   \code{0}). Only used if \code{"mcp"} is used.
#' @param novel Which conditions are to be masked out respectively, either the
#' novel conditions within only \code{"within"} (Default) or also including outside
#' reference conditions \code{"outside"}. Only use for \code{method = "nt2"}, for
#' \code{method = "mess"} this variable is always \code{"within"}.
#' @param limits_clip [`logical`] Should the limits clip all predictors before
#'   fitting a model (\code{TRUE}) or just the prediction (\code{FALSE},
#'   default).
#'
#' @references
#' * Randin, C. F., Dirnböck, T., Dullinger, S., Zimmermann, N. E., Zappa, M., & Guisan, A. (2006). Are niche‐based species distribution models transferable in space?. Journal of biogeography, 33(10), 1689-1703. https://doi.org/10.1111/j.1365-2699.2006.01466.x
#' * Chevalier, M., Broennimann, O., Cornuault, J., & Guisan, A. (2021). Data integration methods to account for spatial niche truncation effects in regional projections of species distribution. Ecological Applications, 31(7), e02427. https://doi.org/10.1002/eap.2427
#' * Velazco, S. J. E., Brooke, M. R., De Marco Jr., P., Regan, H. M., & Franklin, J. (2023). How far can I extrapolate my species distribution model? Exploring Shape, a novel method. Ecography, 11, e06992. https://doi.org/10.1111/ecog.06992
#' * Mesgaran, M. B., R. D. Cousens, B. L. Webber, and J. Franklin. (2014) Here be dragons: a tool for quantifying novelty due to covariate range and correlation change when projecting species distribution models. Diversity and Distributions 20:1147-1159.
#' @returns Adds extrapolation control option to a [`distribution`] object.
#' @keywords control
#' @aliases add_control_extrapolation
#' @examples
#' \dontrun{
#'  # To add a zone layer for extrapolation constraints.
#'  x <- distribution(background) |>
#'    add_predictors(covariates) |>
#'    add_control_extrapolation(method = "zones", layer = zones)
#' }
#' @name add_control_extrapolation
NULL

#' @name add_control_extrapolation
#' @rdname add_control_extrapolation
#' @export
methods::setGeneric(
  "add_control_extrapolation",
  signature = methods::signature("x"),
  function(x, layer, method = "mcp", mcp_buffer = 0,
           novel = "within", limits_clip = FALSE) standardGeneric("add_control_extrapolation"))

#' @name add_control_extrapolation
#' @rdname add_control_extrapolation
methods::setMethod(
  "add_control_extrapolation",
  methods::signature(x = "BiodiversityDistribution"),
  function(x, layer, method = "mcp", mcp_buffer = 0, novel = "within", limits_clip = FALSE) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            missing(layer) || (is.Raster(layer) || inherits(layer, "sf")),
                            (is.numeric(mcp_buffer) && mcp_buffer >=0),
                            is.logical(limits_clip),
                            is.character(novel),
                            is.character(method)
    )
    # Match method
    method <- match.arg(method, c("zones", "mess", "nt2", "mcp", "shape"), several.ok = FALSE)
    novel <- match.arg(novel, c("within", "outside"), several.ok = FALSE)

    # Apply method specific settings
    if(method == "zones"){
      assertthat::assert_that((is.Raster(layer) || inherits(layer, "sf")),
                              msg = "No zone layer specified!")

      if(inherits(layer,'SpatRaster')){
        assertthat::assert_that(terra::is.factor(layer),
                                msg = 'Provided limit raster needs to be ratified (categorical)!')
        layer <- sf::st_as_sf( terra::as.polygons(layer, dissolve = TRUE) ) |> sf::st_cast("MULTIPOLYGON")
      }
      assertthat::assert_that(inherits(layer, "sf"),
                              unique(sf::st_geometry_type(layer)) %in% c('MULTIPOLYGON','POLYGON'),
                              msg = "Limits need to be of polygon type."
      )

      # Get background
      background <- x$background

      # Ensure that limits has the same projection as background
      if(sf::st_crs(layer) != sf::st_crs(background)) layer <- sf::st_transform(layer, background)

      # Ensure that limits is intersecting the background
      if(is.Raster(background)){
        if(suppressMessages(length( sf::st_intersects(layer, terra::as.polygons(background) |> sf::st_as_sf()) )) == 0 ) { layer <- NULL; warning('Provided zones do not intersect the background!') }
      } else {
        if(suppressMessages(length( sf::st_intersects(layer, background |> sf::st_as_sf()) )) == 0 ) { layer <- NULL; warning('Provided zones do not intersect the background!') }
      }

      # Get first column for zone description and rename
      layer <- layer[,1]; names(layer) <- c('limit','geometry')
      limits <- list(layer = layer, "limits_method" = method,
                     "mcp_buffer" = mcp_buffer, "limits_clip" = limits_clip)
      x <- x$set_limits(x = limits)
    } else if(method == "mcp"){
      # Specify the option to calculate a mcp based on the added data.
      # This is done directly in train.
      limits <- list("layer" = NULL, "limits_method" = method,
                     "mcp_buffer" = mcp_buffer, "limits_clip" = limits_clip)
      x <- x$set_limits(x = limits)
    } else if(method == "nt2"){
      # Specify that the multivariate combination novelty index (NT2) is
      # to be applied
      limits <- list("layer" = NULL, "limits_method" = method,
                     "mcp_buffer" = mcp_buffer, "novel" = novel,
                     "limits_clip" = limits_clip)
      x <- x$set_limits(x = limits)
    } else if(method == "mess"){
      # Specify that the multivariate combination novelty index (NT2) is
      # to be applied
      limits <- list("layer" = NULL, "limits_method" = method,
                     "mcp_buffer" = mcp_buffer, "novel" = novel,
                     "limits_clip" = limits_clip)
      x <- x$set_limits(x = limits)
    } else if(method == "shape"){
      stop("Method not yet implemented")
    }

    # Return the altered object
    return(x)
  }
)
