#' Are SpatRasters comparable?
#'
#' This function checks if two `SpatRaster-class` objects
#' are comparable.
#'
#' @param x [`SpatRaster-class`] object.
#' @param y [`SpatRaster-class`] or [`sf`] object.
#' @keywords internal, utils
#' @return [`logical`] indicating if the two [`SpatRaster-class`] objects have the same
#'   resolution, extent, dimensionality, and coordinate system.
#' @noRd
is_comparable_raster <- function(x, y) {
  if(inherits(y, "sf")){
    assertthat::assert_that(
      terra::same.crs(x, y) &&
        sf::st_crs(x) == sf::st_crs(y)
    )
  } else {
    assertthat::assert_that(
      (is.Raster(x) && is.Raster(y)),
      terra::same.crs(x, y) &&
        terra::compareGeom(x, y, stopOnError = FALSE)
    )
  }
}

#' Easy conversion function
#'
#' @description
#' As a consequence of switching to [terra] from [raster], there might be situations
#' where it is necessary to convert between them.
#' This function does the job.
#'
#' @param input A [`SpatRaster`] object to convert to raster.
#' @keywords internal, utils
#' @noRd
terra_to_raster <- function(input){
  assertthat::assert_that(
    is.Raster(input)
  )
  # Check that package is available
  check_package("raster")
  if(!isNamespaceLoaded("raster")) { attachNamespace("raster");requireNamespace("raster") }

  out <- terra::as.raster(input)
  if(raster::nlayers(out)>1) out <- raster::stack(out)
  return(out)
}

#' Do extents intersect?
#'
#' Verify if the extents of two spatial objects intersect or not.
#'
#' @param x [`SpatRaster-class`], [`Spatial-class`] or [`sf::sf()`] object.
#' @param y [`SpatRaster-class`], [`Spatial-class`] or [`sf::sf()`] object.
#' @keywords internal, utils
#' @return [`logical`].
#' @noRd
intersecting_extents <- function(x, y) {
  assertthat::assert_that(
    inherits(x, c("SpatRaster", "Spatial", "sf")),
    inherits(y, c("SpatRaster", "Spatial", "sf")))
  isTRUE(sf::st_intersects(
    terra::ext(x) |> vect() |> sf::st_as_sf(),
    terra::ext(y) |> vect() |> sf::st_as_sf(),
    sparse = FALSE)[[1]])
}

#' Extract polygon data from intersecting point data
#' @param poly A [sf] object.
#' @param points A [`Spatial`] or [sf] object.
#' @param coords A [vector] pointing to the coordinate columns. (Default: \code{c("x", "y")})
#' @keywords utils, internal
#' @return An object with the spatial intersection
#' @noRd
point_in_polygon <- function(poly, points, coords = c('x','y')){
  assertthat::assert_that(
    inherits(poly,'sf'),
    inherits(points,'sf') || inherits(points,'Spatial') || inherits(points,'data.frame'),
    length(coords)>0
  )
  # Convert to sf
  points <- sf::st_as_sf(points, coords = coords, crs = sf::st_crs(poly))
  assertthat::assert_that(
    sf::st_crs(poly) == sf::st_crs(points)
  )

  # Parallize if number of points large and allowed
  if(getOption("ibis.runparallel") && nrow(points) > 5000){
    # Within test
    # Paralise any simple features analysis.
    # @source https://www.spatialanalytics.co.nz/post/2018/04/01/fixing-st-par/
    st_parallel <- function(sf_df, sf_func, n_cores, ...){

      # Create a vector to split the data set up by.
      split_vector <- rep(1:n_cores, each = nrow(sf_df) / n_cores, length.out = nrow(sf_df))

      # FIXME:
      # MC.cores does not work properly on windows. To be replaced with future
      if(Sys.info()['sysname']=="Windows") n_cores <- 1
      # Perform GIS analysis
      split_results <- split(sf_df, split_vector)  |>
        parallel::mclapply(function(x) sf_func(x, ...), mc.cores = n_cores)

      # Define the output_class. If length is greater than two, then grab the second variable.
      output_class <- class(split_results[[1]])
      if (length(output_class) == 2){
        output_class <- output_class[2]
      }

      # Combine results back together. Method of combining depends on the output from the function.
      if (output_class == "matrix"){
        result <- do.call("rbind", split_results)
        names(result) <- NULL
      } else if (output_class == "sfc") {
        result <- do.call("c", split_results)
        result <- sf_func(result) # do.call combines the list but there are still n_cores of the geometry which had been split up. Running st_union or st_collect gathers them up into one, as is the expected output of these two functions.
      } else if (output_class %in% c('list', 'sgbp') ){
        result <- do.call("c", split_results)
        names(result) <- NULL
      } else if (output_class == "data.frame" ){
        result <- do.call("rbind", split_results)
      } else {
        stop("Unknown class. st_parallel only accepts the following outputs at present: sfc, list, sf, matrix, sgbp.")
      }

      # Return result
      return(result)
    }
    ov <- st_parallel(points, function(x) sf::st_join(x, poly, join = st_within), getOption("ibis.nthread"))
  } else {
    # Within test
    ov <- suppressMessages( sf::st_join(points, poly, join = st_within) )
  }
  return(ov)
}

#' Create mask based on a zonal layer
#'
#' @description
#' This function has options to create a mask based on provided point data. It is identical in functionality to
#' the parameter \code{'limit'} in `train()`. Currently it has two available options:
#'
#' [*] It is either possible to provide a categorical zonal [`SpatRaster`] layer takes available point data and intersects it with a
#' zonal layer. The output is a [`SpatRaster`] object with only those classes in which a point occurrence fell.
#' Typical example for instance is layer of the distribution of Biomes, Ecoregions or Climatic Zones.
#'
#' [*] Buffer, in which case a buffer in the units of the geographic projection are created. Buffer width have to
#' be supplied as non-NULL parameter to \code{'buffer_width'}. The output mask thus allows to limit the prediction
#' to a spatial buffer of provided extent within any geovalid occurrence record.
#'
#' @param df A [`sf`] object with point information.
#' @param zones A [`sf`] or [`SpatRaster`] object with polygons of the zones to be used for occurrence masking.
#' @param buffer_width A [`numeric`] value specifying the buffer width. Ignored if a Zones layer is provided.
#' @param column A [`character`] giving the column in which zonal ids are found. Only used when zones is of
#' type [`sf`] (Default: \code{"limits"}).
#' @param template An optional [`SpatRaster`] object on which which the zones should be rasterized (Default: \code{NULL}).
#' @returns A [`sf`] or [`SpatRaster`] object.
#' @keywords utils, internal
#' @noRd
create_zonaloccurrence_mask <- function(df, zones = NULL, buffer_width = NULL, column = "limits", template = NULL){
  assertthat::assert_that(
    inherits(df, "sf"),
    unique(sf::st_geometry_type(df)) %in% "POINT",
    is.character(column),
    is.null(zones) || (inherits(zones, "sf") || is.Raster(zones)),
    is.null(buffer_width) || is.numeric(buffer_width),
    is.null(template) || is.Raster(template),
    # Can't have both set
    !(is.null(zones) && is.null(buffer_width))
  )
  # Make zones mask
  if(!is.null(zones)){
    # If zones is sf, check that it is of type polygon
    if(inherits(zones, "sf")) assertthat::assert_that( all( unique(sf::st_geometry_type(zones)) %in% c("POLYGON", "MULTIPOLYGON") ) )

    if(inherits(zones, "sf")){
      if(sf::st_crs(df)!=sf::st_crs(zones)){
        zones <- zones |> sf::st_transform(crs = sf::st_crs(df))
      }

      # Get zones from the limiting area, e.g. those intersecting with input
      suppressMessages(
        suppressWarnings(
          zones <- sf::st_intersection(df, zones)
        )
      )
      # Extract values from zonal raster layer
      limit <- terra::extract(zones, df) |> unique()

      # Limit zones
      zones <- subset(zones, limit %in% unique(zones[[column]]) )

      # Finally rasterize if template is set
      if(!is.null(template)) zones <- terra::rasterize(zones, template, field = column)
    } else {
      # Extract values from zonal raster layer
      ex <- terra::extract(zones, df) |> unique()
      # Remove NA if found
      if(anyNA(ex)) ex <- ex[-which(is.na(ex))]

      # Now create copy of zonal raster and set all values other than ex to NA
      new <- emptyraster(zones)
      new[zones %in% ex] <- 1
      zones <- new
      # Align with template if set
      if(!is.null(template)){
        if(terra::compareGeom(zones, template, stopOnError = FALSE)){
          zones <- terra::resample(zones, template, method = "near", threads = getOption("ibis.nthread"))
        }
      }
    }
  } else {
    assertthat::assert_that(
      is.Raster(template),msg = "A background layer has to be provided for this function to work!"
    )
    # Buffer points width provided layer
    suppressWarnings(
      buf <- sf::st_buffer(x = df, dist = buffer_width, nQuadSegs = 50)
    )
    # Rasterize
    zones <- terra::rasterize(buf, template, field = 1, background = 0)
    zones <- terra::mask(zones, template)
    # Ratify
    zones <- terra::droplevels(zones)
  }
  return(zones)
}

#' Converts a bounding box to a Well Known Text (WKT) polygon
#'
#' @param minx Minimum x value, or the most western longitude
#' @param miny Minimum y value, or the most southern latitude
#' @param maxx Maximum x value, or the most eastern longitude
#' @param maxy Maximum y value, or the most northern latitude
#' @param all A [`vector`] of length 4, with the elements: minx, miny, maxx, maxy
#' @return An object of class [`character`], a Well Known Text (WKT) string of the form
#' 'POLYGON((minx miny, maxx miny, maxx maxy, minx maxy, minx miny))'
#' @keywords internal, utils
#' @noRd
bbox2wkt <- function(minx=NA, miny=NA, maxx=NA, maxy=NA, bbox=NULL){
  if(is.null(bbox)) bbox <- c(minx, miny, maxx, maxy)
  assertthat::assert_that(length(bbox)==4) #check for 4 digits
  assertthat::assert_that(assertthat::noNA(bbox)) #check for NAs
  assertthat::assert_that(is.numeric(as.numeric(bbox))) #check for numeric-ness

  paste('POLYGON((',
        sprintf('%s %s',bbox[1],bbox[2]), ',', sprintf('%s %s',bbox[3],bbox[2]), ',',
        sprintf('%s %s',bbox[3],bbox[4]), ',', sprintf('%s %s',bbox[1],bbox[4]), ',',
        sprintf('%s %s',bbox[1],bbox[2]),
        '))', sep="")
}

#' Expand an extent by a certain number
#'
#' @param e An [`extent`] object.
#' @param f [`numeric`] value to increase the extent (Default: \code{0.1}).
#' @keywords utils, internal
#' @return Returns the unified total [`extent`] object.
#' @noRd
extent_expand <- function(e,f=0.1){
  assertthat::assert_that(inherits(e,'SpatExtent'),
                          is.numeric(f))
  # Convert to vector
  e <- as.vector(e)
  xi <- (e['xmax']-e['xmin'])*(f/2)
  yi <- (e['ymax']-e['ymin'])*(f/2)

  xmin <- e['xmin']-xi
  xmax <- e['xmax']+xi
  ymin <- e['ymin']-yi
  ymax <- e['ymax']+yi

  return(terra::ext(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)))
}

#' Helper function rename the geometry of a provided
#'
#' @param g A [`sf`] object containing some data.
#' @param name A [`character`] with the new name for the geometry.
#' @source https://gis.stackexchange.com/questions/386584/sf-geometry-column-naming-differences-r
#' @keywords internal, utils
#' @noRd
rename_geometry <- function(g, name){
  assertthat::assert_that(
    inherits(g, "sf"),
    is.character(name)
  )
  current = attr(g, "sf_column")
  if(current == name) return(g)
  names(g)[names(g)==current] = name
  sf::st_geometry(g)=name
  g
}

#' Convert a data.frame or tibble to simple features
#'
#' @description This function tries to guess the coordinate field and converts a data.frame
#' to a simple feature.
#'
#' @param df A [`data.frame`], [`tibble`] or [`sf`] object.
#' @param geom_name A [`character`] indicating the name of the geometry column (Default: \code{'geometry'}).
#' @keywords internal, utils
#' @noRd
guess_sf <- function(df, geom_name = 'geometry'){
 assertthat::assert_that(
   inherits(df,'data.frame') || inherits(df, 'sf') || inherits(df, 'tibble')
 )
 # If sf, return immediately
 if(inherits(df, 'sf')) return(df)
 # If there is an attribute, but for some reason the file is not sf, use that one
 if(!is.null(attr(df, "sf_column"))) {
   df <-  sf::st_as_sf(df)
   if(attr(df, "sf_column") != geom_name){
     names(df)[which(names(df) == attr(df, "sf_column"))] <- geom_name
     sf::st_geometry(df) <- geom_name
   }
   return(df)
 }
 # Commonly used column names
 nx = c("x","X","lon","longitude")
 ny = c("y", "Y", "lat", "latitude")
 ng = c("geom", "geometry", "geometry")

 # Check if geom is present
 if(any( ng %in% names(df) )){
   attr(df, "sf_column") <- ng[which(ng %in% names(df))]
   df <- sf::st_as_sf(df)
 }
 # Finally check if any commonly used coordinate name exist
 if(any( nx %in% names(df))){
   df <- sf::st_as_sf(df, coords = c(nx[which(nx %in% names(df))],
                                     ny[which(ny %in% names(df))])
                      )
 }
 # If at this point df is still not a sf object, then it is unlikely to be converted
 assertthat::assert_that(inherits(df, 'sf'),
                         msg = "Point object could not be converted to an sf object.")
 if(attr(df, "sf_column") != geom_name){
   names(df)[which(names(df) == attr(df, "sf_column"))] <- geom_name
   sf::st_geometry(df) <- geom_name
 }
 return(df)
}

#' Kernel density estimation of coordinates
#'
#' @description
#' Takes input point coordinates as [`sf`] layer and estimates the Gaussian Kernel density over
#' a specified bandwidth for constructing a bivariate Gaussian kernel (see also [`MASS::kde2d()`]).
#' @details
#' Requires the `MASS` R-package to be installed!
#' @param points A \code{POINTS} [`sf`] object.
#' @param background A template [`SpatRaster`] object describing the background.
#' @param bandwidth A [`numeric`] of the input bandwidth (Default \code{2}).
#' @returns A [`SpatRaster`] with the density of point observations.
#' @keywords utils, internal
#' @noRd
st_kde <- function(points, background, bandwidth = 3){
  assertthat::assert_that(
    inherits(points, "sf"),
    is.numeric(bandwidth)
  )
  check_package("MASS")

  # Get extent and cellsize
  cellsize <- terra::res(background)[1]
  extent_vec <- sf::st_bbox(background)[c(1,3,2,4)]

  n_y <- ceiling((extent_vec[4]-extent_vec[3])/cellsize)
  n_x <- ceiling((extent_vec[2]-extent_vec[1])/cellsize)

  extent_vec[2] <- extent_vec[1]+(n_x*cellsize)-cellsize
  extent_vec[4] <- extent_vec[3]+(n_y*cellsize)-cellsize

  # Make
  coords <- sf::st_coordinates(points)
  matrix <- MASS::kde2d(coords[,1],coords[,2],
                        h = bandwidth, n = c(n_x, n_y), lims = extent_vec)

  out <- expand.grid(x = matrix$x, y = matrix$y, KEEP.OUT.ATTRS = FALSE)
  out$z <- as.vector(matrix$z)*(1e11)
  out <- terra::rast(out, crs = terra::crs(background))

  # Resample output for small point mismatches
  if(!terra::compareGeom(out, background, stopOnError = FALSE)){
    out <- terra::resample(out, background, threads = getOption("ibis.nthread"))
  }
  out <- terra::mask(out, background)
  names(out) <- "kde__coordinates"
  rm(matrix, coords)
  return( out )
}

#' Polygon to points
#'
#' @description
#' Converts a polygon [`sf`] layer to a point layer by rasterizing it
#' over a provided [SpatRaster].
#' @param poly A \code{POLYGON} or \code{MULTIPOLYGON} [`sf`] object.
#' @param template A template [`SpatRaster`] object.
#' @param field_occurrence A [`character`] specifying the occurrence field. Should contain information on the type.
#' @keywords utils, internal
#' @noRd
polygon_to_points <- function(poly, template, field_occurrence ) {
  assertthat::assert_that(
    inherits(poly, 'sf'),
    is.Raster(template),
    is.character(field_occurrence),
    assertthat::has_name(poly, field_occurrence)
  )

  # Rasterize the polygon to
  out <- terra::rasterize(x = poly, y = template, field = field_occurrence)

  # Construct new point data
  co <- terra::xyFromCell(out, cell = which(!is.na(out[])) ) |> as.data.frame()
  co[[field_occurrence]] <- out[!is.na(out[])]
  co <- guess_sf(co) # Convert to sf and add coordinates
  co$x <- sf::st_coordinates(co)[,1]
  co$y <- sf::st_coordinates(co)[,2]
  sf::st_crs(co) <- sf::st_crs(template)

  assertthat::assert_that(inherits(co, 'sf'))
  return(co)
}

#' Calculate the dimensions from a provided extent object
#'
#' @description Calculate the dimensions of an extent
#' (either an extent object or four-element vector in the right order), either in projected or spherical space.
#' @param ex Either a [`vector`], a [`SpatExtent`] or alternatively a [`SpatRaster`], [`Spatial*`] or [`sf`] object.
#' @param lonlat A [`logical`] indication whether the extent is WGS 84 projection (Default: \code{TRUE}).
#' @param output_unit [`character`] determining the units. Allowed is 'm' and 'km' (Default: \code{'km'}).
#' @keywords utils, internal
#' @noRd
extent_dimensions <- function(ex, lonlat = terra::is.lonlat(ex), output_unit = 'km') {
  assertthat::assert_that(inherits(ex, 'SpatExtent') || inherits(ex, 'numeric') || inherits(ex, 'sf') || inherits(ex, 'SpatRaster') || inherits(ex, 'Spatial'),
                          is.logical(lonlat),
                          is.character(output_unit) && output_unit %in% c('m','km'))
  # Coerce to vector if necessary
  if(is.Raster(ex)) ex <- terra::ext(ex)
  if(is.vector(ex)) assertthat::assert_that(length(ex)==4, is.numeric(ex),msg = 'No valid extent object supplied!')

  # Convert to vector
  ex <- switch(class(ex)[1],
               Extent = as.vector(ex),
               Raster = as.vector( terra::ext(ex) ),
               sf = as.vector( terra::ext(ex) ),
               numeric = ex
               )
  # Rename the vector
  names(ex) <- c("xmin", "xmax", "ymin", "ymax")

  # Procedures for longlat raster
  if(lonlat) {
    # Dimensions in degrees
    height <- as.numeric( abs(diff(ex[1:2])) )
    width <-  as.numeric( abs(diff(cos(ex[3:4]))) )
    # Scaling to get spherical surface area in km2
    scaling <- (6371 ^ 2 * pi) / 180
    surface_area <- width * height * scaling

    # Ratio between GCD height and width
    # Get ratio between height and width in great circle distance, given an extent vector in lat/long
    lonLatRatio <- function(extent) {
      # lower left point
      p1 <- matrix(extent[c(1, 3)], nrow = 1)
      # upper left and lower right points
      p2 <- rbind(extent[c(1, 4)], extent[c(2, 3)])
      # get ratio between distances
      dists <- terra::distance(p1, p2, lonlat = TRUE)
      ratio <- dists[1] / dists[2]
      return (ratio)
    }
    ratio <- lonLatRatio( as.vector(ex) )
    # calculate equivalent dimensions in km
    w <- sqrt(surface_area / ratio)
    dim <- c(w, w * ratio)
    if(output_unit == 'm') dim * 1000
  } else {
    # else assume a rectangle in m and convert to km
    dim <- abs(diff(ext)[c(1, 3)])
    if(output_unit=='km'){
      dim <- dim * 0.1 ^ 3
    }
  }
  return(dim)
}

#' Align a [`SpatRaster-class`] object to another by harmonizing geometry and extend.
#'
#' If the data is not in the same projection as the template, the alignment
#' will be computed by reprojection only. If the data has already the same
#' projection, the data set will be cropped and aggregated prior to resampling
#' in order to reduce computation time.
#'
#' @param data [`SpatRaster-class`] object to be resampled.
#' @param template [`SpatRaster-class`] or [`sf`] object from which geometry can be extracted.
#' @param method method for resampling (Options: \code{"near"} or \code{"bilinear"}).
#' @param func function for resampling (Default: [mean]).
#' @param cl [`logical`] value if multicore computation should be used (Default: \code{TRUE}).
#' @keywords utils
#' @details
#' Nearest Neighbour resampling (near) is recommended for discrete and bilinear
#' resampling recommended for continuous data. See also help from [terra::resample] for other options.
#' @return New [`SpatRaster`] object aligned to the supplied template layer.
#' @examples
#' \dontrun{
#'  # Align one raster to another
#'  ras1 <- alignRasters( ras1, ras2, method = "near", cl = FALSE)
#' }
#' @export
alignRasters <- function(data, template, method = "bilinear", func = mean, cl = TRUE){
  # Security checks
  assertthat::assert_that(
    is.Raster(data),
    is.Raster(template) || inherits(template, "sf"),
    is.character(method),
    is.logical(cl)
  )
  method <- match.arg(method, c("bilinear", "ngb"), several.ok = FALSE)

  if(sf::st_crs(data) != sf::st_crs(template)){
    # Project Raster layer
    data <- terra::project(data, terra::crs(template), method = method, threads = getOption("ibis.nthread"))
  }

  # Crop raster to template
  data <- terra::crop(data, template, snap = "out")

  # Aggregate to minimal scale
  if(is.Raster(template)){
    if(terra::ncol(data) / terra::ncol(template) >= 2){
      factor <- floor(data@ncols/template@ncols)
      data <- terra::aggregate(data, fact = factor, fun = func, cores = ifelse(cl, getOption("ibis.nthread"), 1))
    }
  } else {
    # Resample with target method
    ras <- terra::rasterize(template, data)
    data <- terra::resample(data, ras, method = method, threads = getOption("ibis.nthread"))
  }
  return(data)
}

#' @title Create an empty \code{SpatRaster} based on a template
#'
#' @description
#' This function creates an empty copy of a provided \code{SpatRaster} object. It
#' is primarily used in the package to create the outputs for the predictions.
#' @param x A \code{SpatRaster*} object corresponding.
#' @param ... other arguments that can be passed to \code{\link{terra}}
#' @return an empty [`SpatRaster`], i.e. all cells are \code{NA}.
#' @import terra
#' @keywords terra, utils
#' @examples
#' require(terra)
#' r <- rast(matrix(1:100, 5, 20))
#' emptyraster(r)
#' @export
emptyraster <- function(x, ...) { # add name, filename,
  assertthat::assert_that(is.Raster(x) || inherits(x, "stars"))
  if(is.Raster(x)){
    terra::rast(nrows = nrow(x), ncols = ncol(x),
                crs = terra::crs(x),
                ext = terra::ext(x), ...)
  } else {
    emptyraster( stars_to_raster(x)[[1]] )
  }
}

#' Function to extract nearest neighbour predictor values of provided points
#'
#' @description
#' This function performs nearest neighbour matching between biodiversity observations and independent
#' predictors, and operates directly on provided data.frames.
#' **Note that despite being parallized this function can be rather slow for large data volumes of data!**
#' @param coords A [`matrix`], [`data.frame`] or [`sf`] object.
#' @param env A [`data.frame`] object with the predictors.
#' @param longlat A [`logical`] variable indicating whether the projection is long-lat.
#' @param field_space A [`vector`] highlight the columns from which coordinates are to be extracted (Default: \code{c('x','y')}).
#' @param cheap A [`logical`] variable whether the dataset is considered to be large and faster computation could help.
#' @param ... other options.
#' @return A [`data.frame`] with the extracted covariate data from each provided data point.
#' @details Nearest neighbour matching is done via the [geodist] R-package (\code{geodist::geodist}).
#' @note If multiple values are of equal distance during the nearest neighbour check, then the results is by default averaged.
#' @examples
#' \dontrun{
#'  # Create matchup table
#' tab <- get_ngbvalue( coords = coords, # Coordinates
#'                      env = env # Data.frame with covariates and coordinates
#'                   )
#' }
#'
#' @references
#' * Mark Padgham and Michael D. Sumner (2021). geodist: Fast, Dependency-Free Geodesic Distance Calculations. R package version 0.0.7. https://CRAN.R-project.org/package=geodist
#' @keywords utils
#' @export
get_ngbvalue <- function(coords, env, longlat = TRUE, field_space = c('x','y'), cheap = FALSE, ...) {
  # Security checks
  assertthat::assert_that(
    is.data.frame(coords) || inherits(coords,'sf') || inherits(coords,'matrix'),
    assertthat::is.flag(longlat),
    is.data.frame(env), assertthat::has_name(env, field_space),
    length(field_space) == 2, is.vector(field_space)
  )
  # Convert to matrices
  coords <- as.matrix(coords)
  coords_env <- as.matrix(env[,field_space])

  # If either of the matrices are larger than 10000 records, process in parallel
  if(is.null( getOption('ibis.runparallel') ) || getOption('ibis.runparallel') == TRUE ){
    process_in_parallel = ifelse(nrow(coords) > 10000 || nrow(coords_env) > 100000, TRUE, FALSE)
  } else {
    process_in_parallel = FALSE
  }

  # Pairwise distance function
  # FIXME: Potentially evaluate whether sf::st_distance is of similar speed for very large matrices.
  # Thus making this dependency suggested and optional
  # disfun <- geosphere::distHaversine
  if(longlat){
    disfun <- function(x1, x2, m = ifelse(cheap,'cheap', 'haversine')) geodist::geodist(x1,x2, measure = m)
  } else {
    disfun <- function(x1, x2, m = NULL) terra::distance(x1, x2, lonlat = longlat)
  }

  if(process_in_parallel){

    check_package("doParallel")

    # Split coordinates into equal size batches of 10
    coords_split <- ggplot2::cut_width(1:nrow(coords),10,boundary=0)

    cl <- doParallel::registerDoParallel(cores = getOption('ibis.nthread'))
    out <- foreach::foreach(z = iterators::iter(unique(coords_split)),
                            .combine = 'rbind',
                            .inorder = FALSE,
                            .multicombine = TRUE,
                            .errorhandling = 'stop',
                            .export = c('coords','coords_env','coords_split', 'disfun'),
                            .packages = c('geodist')
    ) %dopar% {
      o <-
        apply(coords[which(coords_split==z),], 1, function(xy1, xy2){
          dists <- disfun(xy2, xy1)
          # In a few cases these can be multiple in equal distance
          d <- which(dists==min(dists))
          if(length(d)>=2){
            # Average them both
            o <- as.data.frame(
              t(
                apply(env[d, ,drop = FALSE], 2, function(x) mean(x, na.rm = TRUE) )
              )
            )
            return(o)
          } else return( env[d, ,drop = FALSE] )
        }, xy2 = coords_env)
      return(do.call(rbind, o))
    }
    doParallel::stopImplicitCluster()
    rm(cl)
  } else {
    env_sub <- apply(coords, 1, function(xy1, xy2) {
      dists <- disfun(xy2, xy1)
      # In a few cases these can be multiple in equal distance
      d <- which(dists==min(dists))
      if(length(d)>=2){
        # Average them both
        o <- as.data.frame(
          t(
            apply(env[d, ,drop = FALSE], 2, function(x) mean(x, na.rm = TRUE) )
          )
        )
        return(o)
      } else return( env[d, ,drop = FALSE] )
    }, xy2 = coords_env)
    # Combine
    out <- do.call(rbind, env_sub)
    out[,field_space] <- as.data.frame(coords) # Ensure that coordinates are back in
  }
  return(out)
}

#' Function to extract directly the raster value of provided points
#'
#' @description
#' This function simply extracts the values from a provided [`SpatRaster`], [`SpatRasterDataset`] or
#' [`SpatRasterCollection`] object. For points where or NA values were extracted
#' a small buffer is applied to try and obtain the remaining values.
#' @details
#' It is essentially a wrapper for [`terra::extract`].
#' @param coords A [`Spatial`], [`data.frame`], [`matrix`] or [`sf`] object.
#' @param env A [`SpatRaster`] object with the provided predictors.
#' @param ngb_fill [`logical`] on whether cells should be interpolated from neighbouring values.
#' @param rm.na [`logical`] parameter which - if set - removes all rows with a missing data point (\code{NA}) from the result.
#' @return A [`data.frame`] with the extracted covariate data from each provided data point.
#' @keywords utils
#' @examples
#' \dontrun{
#' # Extract values
#' vals <- get_rastervalue(coords, env)
#' }
#' @export
get_rastervalue <- function(coords, env, ngb_fill = TRUE, rm.na = FALSE){
  assertthat::assert_that(
    inherits(coords,"sf") || inherits(coords, "Spatial") || (is.data.frame(coords) || is.matrix(coords)),
    is.Raster(env),
    is.logical(ngb_fill),
    is.logical(rm.na)
    )

  # Try an extraction
  ex <- try({terra::extract(x = env,
                            y = coords,
                            method = "simple")}, silent = FALSE)
  if(inherits(ex, "try-error")) stop(paste("SpatRaster valueextraction failed: ", ex))
  # Find those that have NA in there
  check_again <- apply(ex, 1, function(x) anyNA(x))
  if(any(check_again)){
    # Re-extract but with a small buffer
    if(inherits(coords, "sf")){
      coords_sub <- coords[check_again,]
    } else {
      coords_sub <- coords[which(check_again),] |> as.data.frame()
    }
    if(any(is.factor(env))) ngb_fill <- FALSE
    ex_sub <- try({terra::extract(x = env,
                                  y = coords_sub,
                                  # FIXME: This could fail if values are factors?
                                  method = ifelse(ngb_fill, "bilinear", "simple"),
                                  touches = TRUE)}, silent = FALSE)
    if(inherits(ex_sub, "try-error")) stop(paste("Raster extraction failed!"))
    ex[which(check_again),] <- ex_sub[, names(ex)]
  }
  # Add coordinate fields to the predictors as these might be needed later
  if(!any(assertthat::has_name(ex, c("x", "y")))){
    if(inherits(coords,"sf")) coords <- sf::st_coordinates(coords)
    ex[["x"]] <- as.numeric(coords[,1]); ex[["y"]] <- as.numeric(coords[,2])
  }
  # Convert to factor if any
  if(any(is.factor(env))){
    ex[,names(env)[which(is.factor(env))]] <- factor(ex[,names(env)[which(is.factor(env))]])
  }

  if(rm.na){
    ex <- subset(ex, stats::complete.cases(ex))
  }
  assertthat::assert_that(is.data.frame(ex),
                          nrow(ex)>0,
                          msg = "Something went wrong with the extraction or all points had missing data.")
  return(ex)
}

#' Create new raster stack from a given data.frame
#'
#' @param post A data.frame
#' @param background A [`SpatRaster-class`] object for the background raster.
#' @keywords internal, utils
#' @return A [`SpatRaster-class`] object with number of columns equal to \code{ncol(post)}.
#' @noRd
fill_rasters <- function(post, background){
  assertthat::assert_that(
    is.data.frame(post),
    is.Raster(background),
    nrow(post) == terra::ncell(background)
  )
  # Make names to be sure
  names(post) <- base::make.names(names(post))

  # If only one raster
  if(ncol(post)==1){
    out <- emptyraster(background)
    out[] <- post[,1]
  } else {
    # Loop through each column
    out <- terra::rast()
    for(co in 1:ncol(post)){
      o <- emptyraster(background)
      o[] <- post[,co] # Assign values
      # Add to stack
      suppressWarnings( out <- c(out, o) )
      rm(o)
    }
  }
  # Assign names
  names(out) <- names(post)

  # Check that derivate attributes if existing are passed
  if(length( grep("deriv", names(attributes(post)) ))>0){
    attr(out, grep("deriv", names(attributes(post)), value = TRUE) ) <- attr(post,
                                                                             grep("deriv", names(attributes(post)),
                                                                                  value = TRUE) )
  }

  # Final check
  assertthat::assert_that(
    is.Raster(out),
    terra::nlyr(out) == ncol(post)
  )
  return(out)
}

#' Create a polynomial transformation from coordinates
#'
#' @description This function transforms the coordinates of a supplied file through a polynomial transform.
#' By default it applies weights and a QR decomposition for numerical stability.
#' @param coords A [`data.frame`], [`matrix`] or [`sf`] object with coordinates (2 columns named x-y).
#' @param degree The number of degrees used for polynominal transformation (Default: \code{2}).
#' @param weights Set by default to the inverse of the number of coordinates.
#' @returns A data.frame with transformed coordinates.
#' @keywords utils
#' @keywords internal
#' @references
#' * Dray S., Plissier R., Couteron P., Fortin M.J., Legendre P., Peres-Neto P.R., Bellier E., Bivand R., Blanchet F.G., De Caceres M., Dufour A.B., Heegaard E., Jombart T., Munoz F., Oksanen J., Thioulouse J., Wagner H.H. (2012). Community ecology in the age of multivariate multiscale spatial analysis. Ecological Monographs 82, 257–275.
#' @noRd
polynominal_transform <- function(coords, degree = 2, weights = rep(1/nrow(coords), nrow(coords)) ){
  assertthat::assert_that(
    inherits(coords, 'data.frame') || inherits(coords, 'matrix') || inherits(coords, 'sf'),
    is.numeric(degree),
    !is.null(weights) && length(weights) == nrow(coords)
  )
  # If spatial get coordinates
  if(inherits(coords, 'sf')){
    coords <- sf::st_coordinates(coords)
  }
  # Polynomial transform
  a0 <- stats::poly(x = as.matrix( coords ), degree = degree, simple = TRUE)
  # Standardize colnames
  poly.names <- colnames(a0) # Column names for later
  poly.names <- paste0("spatialtrend_", gsub("\\.","_",poly.names) )

  # Standardize the weights
  weights <- weights/sum(weights)
  a0 <- cbind(weights, a0) # Add to polynominal transform
  a0 <- base::qr.Q(base::qr(a0)) # QR decomposition for better numerical stability
  a0 <- as.data.frame(a0[, -1])/sqrt(weights) # Weighting

  # Rename
  colnames(a0) <- poly.names
  return(a0)
}

#' Clean up raster layer from disk
#'
#' Completely deletes for instance a temporary created [`SpatRaster`] object.
#' @param A [`SpatRaster`] object.
#' @param verbose Print progress (Default: \code{FALSE}).
#' @keywords utils
#' @noRd
clean_rasterfile <- function(x, verbose = FALSE)
{
  stopifnot(grepl("SpatRaster", class(x)))
  sink(tempfile())
  tdir = terra::terraOptions()[["tmpdir"]]
  sink(NULL)
  if (inherits(x, "SpatRaster"))
    files = basename( x@ptr$filenames()[1] )
  if (inherits(x, "SpatRaster") & terra::nlyr(x)>1)
    files = do.call(c, lapply(methods::slot(x, "layers"),
                              function(x) x@ptr$filenames()[1]))
  files = files[file.exists(files)]
  if (length(files) == 0)
    return(NULL)
  lapply(files, function(f) {
    if (file.exists(f))
      file.remove(f, sub("grd", "gri", f))
    if (verbose) {
      print(paste("Deleted: ", f))
      print(paste("Deleted: ", sub("grd", "gri",
                                   f)))
    }
  })
  parent.var.name <- deparse(substitute(x))
  rm(list = parent.var.name, envir = sys.frame(-1))
}

#' Split raster factor levels to stack
#'
#'  @description Takes a single raster that is a [`factor`] and creates
#'  a new [`SpatRaster`] that contains the individual levels.
#'  @param ras A [`SpatRaster`] object that is a [`factor`]. Alternatively a [`SpatRaster`] object
#'  can be supplied in which only factor variables are 'exploded'.
#'  @param name An optional [`character`] name for the [`SpatRaster`].
#'  @param ... Other parameters (not used).
#'  @returns A [`SpatRaster`] object.
#'  @keywords utils, internal
#'  @noRd
explode_factorized_raster <- function(ras, name = NULL, ...){
  assertthat::assert_that(is.Raster(ras),
                          is.null(name) || is.character(name))

  # Simply return the input if there are no factors
  if(!any(is.factor(ras))) return(ras)

  # If input is a SpatRaster
  if(inherits(ras, 'SpatRaster')){
    # Get name
    # Create output template
    temp <- emptyraster(ras)
    if(is.null(name)) name <- names(ras)

    # Extract data
    o <- data.frame(val = values(ras));names(o) <- name;o[[name]] <- factor(o[[name]])

    # Check if there is an NaN and remove if present.
    if( any(is.nan( terra::levels(o[[name]]) ) | terra::levels(o[[name]]) == "NaN") ){
      lvl <- terra::levels(o[[name]])
      if(any(is.nan(lvl))) lvl <- lvl[-which(is.nan(lvl))]
      if(any(lvl == "NaN")) lvl <- lvl[-which(lvl == "NaN")]
    } else {
      lvl <- terra::levels(o[[name]])
    }

    # Make function that converts all factors to split rasters
    f <- as.data.frame(
      outer(o[[name]], lvl, function(w, f) ifelse(w == f, 1, 0))
    )

    # Fill template rasters
    out <- fill_rasters(f, temp)
    names(out) <- paste(name, lvl, sep = ".")

  } else if(terra::nlyr(ras)>1){
    # Alternatively if input is stack
    fcts <- is.factor(ras)

    # Get non-factor variables
    if(length( which(!fcts) ) >0){
      out <- ras[[which(!fcts)]]
    } else {
      out <- terra::rast()
    }
    for(k in which(fcts)){

      sub <- ras[[k]]

      temp <- emptyraster(sub)
      if(is.null(name)) new_name <- names(sub)

      # Extract data
      o <- data.frame(val = values(sub));names(o) <- new_name;o[[new_name]] <- factor(o[[new_name]])

      # Check if there is an NaN and remove if present.
      if( any(is.nan(terra::levels(o[[name]])) | terra::levels(o[[name]]) == "NaN") ){
        lvl <- terra::levels(o[[name]])
        if(any(is.nan(lvl))) lvl <- lvl[-which(is.nan(lvl))]
        if(any(lvl == "Nan")) lvl <- lvl[-which(lvl == "NaN")]
      } else {
        lvl <- terra::levels(o[[name]])
      }

      # Make function that converts all factors to split rasters
      f <- as.data.frame(
        outer(o[[new_name]], lvl, function(w, f) ifelse(w == f, 1, 0))
      )

      # Fill template rasters
      new <- fill_rasters(f, temp)
      names(new) <- paste(new_name, lvl, sep = ".")
      suppressWarnings( out <- c(out, new) )
    }
  }
  return(out) # Return the result
}

#' Functionality for geographic and environmental thinning
#'
#' @description
#' For most species distribution modelling approaches it is assumed that occurrence records are unbiased, which
#' is rarely the case. While model-based control can alleviate some of the effects of sampling bias, it can often be
#' desirable to account for some sampling biases through spatial thinning (Aiello‐Lammens et al. 2015). This
#' is an approach based on the assumption that oversampled grid cells contribute little more than bias, rather than
#' strengthing any environmental responses.
#' This function provides some methods to apply spatial thinning approaches. Note that this effectively removes
#' data prior to any estimation and its use should be considered with care (see also Steen et al. 2021).
#'
#' @details
#' Currently implemented thinning methods:
#'
#'  * \code{"random"}: Samples at random up to number of \code{"minpoints"} across all occupied grid cells.
#'  Does not account for any spatial or environmental distance between observations.
#'  * \code{"bias"}: This option removed explicitly points that are considered biased (parameter \code{"env"}) only.
#'  Points are preferentially thinned from grid cells which are in the 25% most biased (larger values assumed greater bias)
#'  and have high point density. Thins the observations up to \code{"minpoints"}.
#'  * \code{"zones"}: Assesses for each observation that it falls with a maximum of \code{"minpoints"} into
#'  each occupied zone. Careful: If the zones are relatively wide this can remove quite a few observations.
#'  * \code{"environmental"}: This approach creates an observation-wide clustering (k-means) under the assumption
#'  that the full environmental niche has been comprehensively sampled and is covered by the provided covariates \code{env}.
#'  We then obtain an number equal to (\code{"minpoints"}) of observations for each cluster.
#'  * \code{"spatial"}: Calculates the spatial distance between all observations. Then points are removed
#'  iteratively until the minimum distance between points is crossed.  The \code{"mindistance"} parameter has to
#'  be set for this function to work.
#'
#' @param df A [`sf`] or [`data.frame`] object with observed occurrence points. All methods threat presence-only
#' and presence-absence occurrence points equally.
#' @param background A [`SpatRaster`] object with the background of the study region. Use for assessing point density.
#' @param env A [`SpatRaster`] object with environmental covariates. Needed when method is set to \code{"environmental"}
#' or \code{"bias"} (Default: \code{NULL}).
#' @param method A [`character`] of the method to be applied (Default: \code{"random"}).
#' @param minpoints A [`numeric`] giving the number of data points at minimum to take (Default: \code{10}).
#' @param mindistance A [`numeric`] for the minimum distance of neighbouring observations (Default: \code{NULL}).
#' @param zones A [`SpatRaster`] to be supplied when option \code{"method"} is chosen (Default: \code{NULL}).
#' @param verbose [`logical`] of whether to print some statistics about the thinning outcome (Default: \code{TRUE}).
#' @examples
#' \dontrun{
#'  # Thin a certain number of observations
#'  # At random
#'  thin_points <- thin_observations(points, background, method = "random")
#'  # using a bias layer
#'  thin_points <- thin_observations(points, background, method = "bias", env = bias)
#' }
#' @references
#' * Aiello‐Lammens, M. E., Boria, R. A., Radosavljevic, A., Vilela, B., & Anderson, R. P. (2015). spThin: an R package for spatial thinning of species occurrence records for use in ecological niche models. Ecography, 38(5), 541-545.
#' * Steen, V. A., Tingley, M. W., Paton, P. W., & Elphick, C. S. (2021). Spatial thinning and class balancing: Key choices lead to variation in the performance of species distribution models with citizen science data. Methods in Ecology and Evolution, 12(2), 216-226.
#' @keywords utils
#' @export
thin_observations <- function(df, background, env = NULL, method = "random", minpoints = 10, mindistance = NULL,
                              zones = NULL, verbose = TRUE){
  assertthat::assert_that(
    inherits(df, "sf") || inherits(df, "data.frame"),
    nrow(df) > 0,
    is.Raster(background),
    is.Raster(env) || is.null(env),
    is.character(method),
    is.numeric(minpoints) && minpoints > 0,
    is.null(mindistance) || is.numeric(mindistance),
    is.Raster(zones) || is.null(zones)
  )
  check_package("dplyr")
  # Match method
  method <- match.arg(method, choices = c("random", "spatial", "bias", "environmental", "zones"), several.ok = FALSE)

  # Label background with id
  bg <- background
  bg[] <- 1:terra::ncell(bg)
  bg <- terra::mask(bg, background)

  # Check that environment has the same projection
  if(is.Raster(env) && method == "environmental"){
    assertthat::assert_that( terra::compareGeom(bg, env, stopOnError = FALSE) )
  }
  # Check that CRS is the same as background
  if(sf::st_crs(df) != sf::st_crs(bg)){
    message("Projection is different from input data. Reprojecting!")
    df <- df |> sf::st_transform(crs = sf::st_crs(bg))
  }

  # Take coordinates of supplied data and rasterize
  coords <- sf::st_coordinates( df )
  ras <- terra::rasterize(coords, bg) # Get the number of observations per grid cell

  # Bounds for thining
  totake <- c(lower = minpoints, upper = max( terra::global(ras, "min", na.rm = TRUE)[,1], minpoints))

  # -- #
  if(method == "random"){
    # For each unique grid cell id, get the minimum value up to a maximum of the points
    # by sampling at random from the occupied grid cells

    # Output vector
    sel <- vector()

    ex <- data.frame(id = 1:nrow(coords),
                     cid = terra::extract(bg, coords)[,1]
    )
    ex <- subset(ex, stats::complete.cases(ex)) # Don't need missing points

    ex <- dplyr::left_join(ex,
                           ex |> dplyr::group_by(cid) |> dplyr::summarise(N = dplyr::n()),
                           by = "cid"
    )
    # Points to take
    sel <- append(sel, ex$id[which(ex$N <= min(totake))] )

    # For those where we have more than the minimum, take at random the upper limits of observations
    ex$oversampled <- ifelse(ex$N >= totake["upper"], 1, 0)
    if(dplyr::n_distinct(ex$oversampled) > 1){
      # If there any oversampled
      # Now sample at random up to the maximum amount. Got tired of doing this outside tidyverse
      o <- ex  |> dplyr::filter(oversampled == 1) |>
        dplyr::group_by(cid) |>
        dplyr::slice_sample(n = min(totake))
      if(nrow(o)>0) sel <- append(sel, o$id)
      rm(o)
    }
    if(anyDuplicated(sel)) sel <- unique(sel)

    try({rm(ex)},silent = TRUE)
  } else if(method == "bias"){
    assertthat::assert_that(is.Raster(env),
                            terra::nlyr(env)==1,
                            msg = "Bias requires a single SpatRaster layer given to env.")
    sel <- vector()

    # Convert bias layer into percentile (largest being)
    bias_perc <- terra::global(env, fun = quantile, na.rm = TRUE)[["X75."]]

    # Now extract
    ex <- data.frame(id = 1:nrow(coords),
                     cid = terra::extract(bg, coords)[,1],
                     pres = terra::extract(ras, coords)[,1],
                     bias = terra::extract(env, coords)[,1]
    )
    ex <- subset(ex, stats::complete.cases(ex)) # Don't need missing points
    # Now identify those to be thinned
    ex$tothin <- ifelse((ex$bias >= bias_perc) & (ex$pres < totake[1]), 1, 0)
    # Now thin those points that are to be thinned
    if(length(unique(ex$tothin))>1){
      ss <- ex |> dplyr::filter(tothin == 1) |>
        dplyr::group_by(cid) |>
        dplyr::slice_sample(n = totake[1], weight_by = bias, replace = T) |>
        dplyr::distinct()

      # Points to take
      sel <- append(sel, ex$id[ex$tothin==0] )
      sel <- append(sel, ss$id )
      try({rm(ss, ex)},silent = TRUE)
    }

  } else if(method == "zones"){
    # Thinning by zones
    assertthat::assert_that(is.Raster(zones),
                            is.factor(zones))

    if(!terra::compareGeom(bg, zones, stopOnError = FALSE)){
      zones <- alignRasters(zones, bg, method = "near", func = terra::modal, cl = FALSE)
    }

    # Output vector
    sel <- vector()

    ex <- data.frame(id = 1:nrow(coords),
                     cid = terra::extract(bg, coords)[,1],
                     zones = terra::extract(zones, coords)[,1]
    )
    # Now for each zone, take the minimum amount at random
    ss <- ex |>
      dplyr::group_by(zones) |>
      dplyr::slice_sample(n = max(totake[1]), replace = TRUE) |>
      dplyr::distinct()

    # Take the zone data points
    sel <- append(sel, ss$id )
    try({rm(ss, ex)},silent = TRUE)

  } else if(method == "environmental"){
    # Environmental clustering

    if(!terra::compareGeom(bg, env, stopOnError = FALSE)){
      env <- alignRasters(env, bg, method = "near", func = terra::modal, cl = FALSE)
    }
    # If there are any factors, explode
    if(any(is.factor(env))){
      env <- explode_factorized_raster(env)
    }

    # Output vector
    sel <- vector()

    # Get a matrix of all environmental data, also with coordinates
    # However first normalize all data
    stk <- terra::as.data.frame(
        predictor_transform(env, option = "norm"),
      xy = TRUE)

    stk$cid <- 1:nrow(stk)
    stk <- subset(stk, stats::complete.cases(stk))

    # Cluster
    E <- stats::kmeans(x = subset(stk, select = -cid),
                       centers = ncol(stk)-1, iter.max = 10)

    stk$cluster <- E$cluster

    # Now fill an empty raster and re-xtract
    new <- emptyraster(env)
    new[stk$cid] <- stk$cluster

    # Now re-extract and sampling points
    ex <- data.frame(id = 1:nrow(coords),
                     cid = terra::extract(bg, coords)[,1],
                     zones = terra::extract(new, coords)[,1]
    )

    # Now for each zone, take the minimum amount at random
    ss <- ex |>
      dplyr::group_by(zones) |>
      dplyr::slice_sample(n = max(totake[1]), replace = TRUE) |>
      dplyr::distinct()

    # Take the zone data points
    sel <- append(sel, ss$id )

    try({rm(new, stk, ss, ex, E)},silent = TRUE)
  } else if(method == "spatial"){
    # Spatial thinning
    stop("Not yet implemented!")
  }

  # Return subsampled coordinates
  out <- df[sel,]
  if(nrow(out)==0) {
    message("Thinning failed for some reason")
    return(df)
  } else {
    if(verbose){
      message(paste0(
        "(", method, ")", " thinning completed! \n",
        "Original number of records: ", nrow(df), "\n",
        "Number of retained records: ", nrow(out))
      )
    }
    return(out)
  }
}
