#' Are rasters comparable?
#'
#' This function checks if two `Raster-class` objects
#' are comparable.
#'
#' @param x [`Raster-class`] object.
#' @param y [`Raster-class`] object.
#' @keywords internal, utils
#' @return [`logical`] indicating if the two [`Raster-class`] objects have the same
#'   resolution, extent, dimensionality, and coordinate system.
is_comparable_raster <- function(x, y) {
  assertthat::assert_that(inherits(x, "Raster"), inherits(y, "Raster")) &&
    sf::st_crs(x@crs) == sf::st_crs(y@crs) &&
    raster::compareRaster(x, y, crs = FALSE, res = TRUE, tolerance = 1e-5,
                          stopiffalse = FALSE)
}

#' Do extents intersect?
#'
#' Verify if the extents of two spatial objects intersect or not.
#'
#' @param x [`Raster-class`], [`Spatial-class`] or [`sf::sf()`] object.
#' @param y [`Raster-class`], [`Spatial-class`] or [`sf::sf()`] object.
#' @keywords internal, utils
#' @return [`logical`].
intersecting_extents <- function(x, y) {
  assertthat::assert_that(
    inherits(x, c("Raster", "Spatial", "sf")),
    inherits(y, c("Raster", "Spatial", "sf")))
  isTRUE(sf::st_intersects(
    sf::st_as_sf(methods::as(raster::extent(x), "SpatialPolygons")),
    sf::st_as_sf(methods::as(raster::extent(y), "SpatialPolygons")),
    sparse = FALSE)[[1]])
}


#' Extract polygon data from intersecting point data
#' @param poly A [sf] object
#' @param points A [`Spatial`] or [sf] object
#' @param coords A [vector] pointing to the coordinate columns. (Default: \code{c("x", "y")})
#' @keywords utils
#' @return An object with the spatial intersection
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
      split_results <- split(sf_df, split_vector) %>%
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
#' This function takes available point data and intersects it with a zonal layer
#' of the
#' @param df A [`sf`] object with point information.
#' @param zones A [`sf`] or [`RasterLayer`] object with polygons of the zones to be used for occurrence masking.
#' @param column A [`character`] giving the column in which zonal ids are found. Only used when zones is of
#' type [`sf`].  (Default: \code{"limits"}).
#' @param template An optional [`RasterLayer`] object on which which the zones should be rasterized (Default: \code{NULL}).
#' @returns A [`sf`] or [`RasterLayer`] object.
#' @keywords utils
#' @noRd
create_zonaloccurrence_mask <- function(df, zones, column = "limits", template = NULL){
  assertthat::assert_that(
    inherits(df, "sf"),
    unique(sf::st_geometry_type(df)) %in% "POINT",
    is.character(column),
    inherits(zones, "sf") || is.Raster(zones),
    is.null(template) || is.Raster(template)
  )
  # If zones is sf, check that it is of type polygon
  if(inherits(zones, "sf")) assertthat::assert_that( all( unique(sf::st_geometry_type(zones)) %in% c("POLYGON", "MULTIPOLYGON") ) )

  if(inherits(zones, "sf")){
    # Get zones from the limiting area, e.g. those intersecting with input
    suppressMessages(
      suppressWarnings(
        zones <- sf::st_intersection(df, zones)
      )
    )
    # Limit zones
    zones <- subset(zones, limit %in% unique(zones[[column]]) )

    # Finally rasterize if template is set
    if(!is.null(template)) zones <- raster::rasterize(zones, template, field = column)
  } else {
    # Extract values from zonal raster layer
    ex <- raster::extract(zones, df) |> unique()
    # Remove NA if found
    if(anyNA(ex)) ex <- ex[-which(is.na(ex))]

    # Now create copy of zonal raster and set all values other than ex to NA
    new <- emptyraster(zones)
    new[zones %in% ex] <- 1
    zones <- new
    # Align with template if set
    if(!is.null(template)){
      zones <- alignRasters(zones, template, method = "ngb", func = raster::modal, cl = FALSE)
    }
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
bbox2wkt <- function(minx=NA, miny=NA, maxx=NA, maxy=NA, bbox=NULL){
  if(is.null(bbox)) bbox <- c(minx, miny, maxx, maxy)
  assertthat::assert_that(length(bbox)==4) #check for 4 digits
  assertthat::assert_that(noNA(bbox)) #check for NAs
  assertthat::assert_that(is.numeric(as.numeric(bbox))) #check for numeric-ness

  paste('POLYGON((',
        sprintf('%s %s',bbox[1],bbox[2]), ',', sprintf('%s %s',bbox[3],bbox[2]), ',',
        sprintf('%s %s',bbox[3],bbox[4]), ',', sprintf('%s %s',bbox[1],bbox[4]), ',',
        sprintf('%s %s',bbox[1],bbox[2]),
        '))', sep="")
}

#' Expand an extent by a certain number
#' @param e an [`extent`] object
#' @param f [`numeric`] value to increase the extent (Default: \code{0.1})
#' @keywords utils
#' @return Returns the unified total [`extent`] object
extent_expand <- function(e,f=0.1){
  assertthat::assert_that(inherits(e,'Extent'))
  xi <- (e@xmax-e@xmin)*(f/2)
  yi <- (e@ymax-e@ymin)*(f/2)

  xmin <- e@xmin-xi
  xmax <- e@xmax+xi
  ymin <- e@ymin-yi
  ymax <- e@ymax+yi

  return(extent(c(xmin,xmax,ymin,ymax)))
}

#' Helper function rename the geometry of a provided
#'
#' @param g A [`sf`] object containing some data.
#' @param name A [`character`] with the new name for the geometry.
#' @source https://gis.stackexchange.com/questions/386584/sf-geometry-column-naming-differences-r
#' @keywords internal
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
#' to a simpel feature
#' @param df A [`data.frame`], [`tibble`] or [`sf`] object
#' @param geom_name A [`character`] indicating the name of the geometry column. Default: 'geometry'
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

#' Polygon to points
#'
#' @description
#' Converts a polygon [`sf`] layer to a point layer by rasterizing it
#' over a provided Raster layer.
#' @param poly A \code{POLYGON} or \code{MULTIPOLYGON} [`sf`] object.
#' @param template A template [`Raster`] object.
#' @param field_occurrence A [`character`] specifying the occurrence field. Should contain information on the type.
#' @keywords utils
#' @noRd
polygon_to_points <- function(poly, template, field_occurrence ) {
  assertthat::assert_that(
    inherits(poly, 'sf'),
    is.Raster(template),
    is.character(field_occurrence),
    assertthat::has_name(poly, field_occurrence)
  )

  # Rasterize the polygon to
  out <- raster::rasterize(poly, template, field = field_occurrence)

  # Construct new point data
  co <- raster::xyFromCell(out, cell = which(!is.na(out[])) ) |> as.data.frame()
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
#' (either an extent object orfour-element vector in the right order), either in projected or spherical space
#' @param ex Either a [`vector`], a [`extent`] or alternatively a [`Raster`],[`Spatial*`] or [`sf`] object
#' @param lonlat A [`logical`] indication whether the extent is WGS 84 projection (Default: TRUE)
#' @param output_unit [`character`] determining the units. Allowed is 'm' and 'km' (Default: 'km')
#' @keywords utils
#' @noRd
extent_dimensions <- function(ex, lonlat = TRUE, output_unit = 'km') {
  assertthat::assert_that(inherits(ex, 'Extent') || inherits(ex, 'numeric') || inherits(ex, 'sf') || inherits(ex, 'Raster') || inherits(ex, 'Spatial'),
                          is.logical(lonlat),
                          is.character(output_unit) && output_unit %in% c('m','km'))
  # Coerce to vector if necessary
  if(is.Raster(ex)) ex <- raster::extent(ex)
  if(is.vector(ex)) assertthat::assert_that(length(ex)==4, is.numeric(ex),msg = 'No valid extent object supplied!')

  # Convert to vector
  ex <- switch(class(ex)[1],
               Extent = as.vector(ex),
               Raster = as.vector( raster::extent(ex) ),
               sf = as.vector( raster::extent(ex) ),
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
      dists <- raster::pointDistance(p1,p2,lonlat = TRUE)
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
    dim <- abs(diff(extent)[c(1, 3)])
    if(output_unit=='km'){
      dim <- dim * 0.1 ^ 3
    }
  }
  return(dim)
}

#' Align a [`Raster-class`] object to another by harmonizing geometry and extend.
#'
#' If the data is not in the same projection as the template, the alignment
#' will be computed by reprojection only. If the data has already the same
#' projection, the data set will be cropped and aggregated prior to resampling
#' in order to reduce computation time.
#'
#' @param data [`Raster-class`] object to be resampled
#' @param template [`Raster-class`] or [`Spatial-class`] object from which geometry can be extracted
#' @param method method for resampling (Options: \code{"ngb"} or \code{"bilinear"})
#' @param func function for resampling (Default: [mean])
#' @param cl [`logical`] value if multicore computation should be used (Default: \code{TRUE})
#' @keywords utils
#' @details
#' Nearest Neighbour resampling (ngb) is recommended for discrete and Bilinear
#' resampling for continuous data.
#' @return New [`Raster`] object aligned to the supplied template layer
alignRasters <- function(data, template, method = "bilinear",func = mean,cl = TRUE){
  # Security checks
  assertthat::assert_that(
    inherits(data,'Raster'), inherits(template, c("Raster", "Spatial", "sf")),
    is.character(method),
    is.logical(cl)
  )
  # Start cluster if necessary
  if(cl) raster::beginCluster(parallel::detectCores()-1)
  if(raster::projection(data) == raster::projection(template)){
    # Crop raster to template
    data <- raster::crop(data, template, snap = "out")
    if(class(template) == "RasterLayer"){
      # Aggregate to minimal scale
      if(data@ncols / template@ncols >= 2){
        factor <- floor(data@ncols/template@ncols)
        data <- aggregate(data, fact = factor, fun = func,
                          expand=TRUE)
      }
      # Resample with target method
      data <- raster::resample(data, template, method = method)
    }
  } else {
    # Project Raster layer
    data <- projectRaster(data, template, method = method)
  }
  # Stop cluster
  if(cl) endCluster()
  return(data)
}

#' @title Create an empty \code{rasterLayer} based on a template
#'
#' @param x a \code{raster*} object corresponding.
#' @param ... other arguments that can be passed to \code{\link{raster}}
#' @return an empty raster, i.e. all cells are \code{NA}.
#' @import raster
#' @keywords raster, utils
#' @examples
#' require(raster)
#' r <- raster(matrix(1:100, 5, 20))
#' emptyraster(r)
#' @export
emptyraster <- function(x, ...) { # add name, filename,
  assertthat::assert_that(is.Raster(x))
  raster::raster(nrows = nrow(x), ncols = ncol(x),
                        crs = x@crs,
                        ext = raster::extent(x), ...)
}

#' Function to extract nearest neighbour predictor values of provided points
#'
#' @description
#' This function performs nearest neighbour matching between biodiversity observations and independent
#' predictors, and operates directly on provided data.frames.
#' **Note that despite being parallized this function can be rather slow for large data volumes of data!**
#' @param coords A [`matrix`], [`data.frame`] or [`sf`] object.
#' @param env A [`data.frame`] object with the predictors
#' @param longlat A [`logical`] variable indicating whether the projection is long-lat
#' @param field_space A [`vector`] highlight the columns from which coordinates are to be extracted (default: \code{c('x','y')})
#' @param cheap A [`logical`] variable whether the dataset is considered to be large and faster computation could help
#' @return A [`data.frame`] with the extracted covariate data from each provided data point.
#' @details Nearest neighbour matching is done via the [geodist] R-package (\code{geodist::geodist})
#' @note If multiple values are of equal distance during the nearest neighbour check, then the results is by default averaged.
#' @references
#' * Mark Padgham and Michael D. Sumner (2021). geodist: Fast, Dependency-Free Geodesic Distance Calculations. R package version 0.0.7. https://CRAN.R-project.org/package=geodist
#' @keywords utils
#' @export
get_ngbvalue <- function(coords, env, longlat = TRUE, field_space = c('x','y'), cheap = FALSE, ...) {
  # Security checks
  assertthat::assert_that(
    is.data.frame(coords) || inherits(coords,'sf') || inherits(coords,'matrix'),
    assertthat::is.flag(longlat),
    is.data.frame(env),assertthat::has_name(env, field_space),
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
    disfun <- function(x1,x2, m = ifelse(cheap,'cheap','haversine')) geodist::geodist(x1,x2, measure = m)
  } else {
    disfun <- function(x1, x2) raster::pointDistance(x1, x2, lonlat = longlat)
  }

  if(process_in_parallel){
    check_package("doParallel")
    suppressPackageStartupMessages(require(doParallel))

    # Split coordinates into equal size batches of 10
    coords_split <- ggplot2::cut_width(1:nrow(coords),10,boundary=0)

    cl <- doParallel::registerDoParallel(cores = getOption('ibis.nthread'))
    out <- foreach(z = iterators::iter(unique(coords_split)),
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
#' This function simply extracts the values from a provided [`RasterLayer`],
#' [`RasterStack`] or [`RasterBrick`] object. For points where or NA values were extracted
#' a small buffer is applied to try and obtain the remaining values.
#' @param coords A [`Spatial`] or [`sf`] object.
#' @param env A [`Raster`] object with the provided predictors.
#' @param rm.na [`logical`] parameter which - if set - removes all rows with a missing data point (\code{NA}) from the result.
#' @return A [`data.frame`] with the extracted covariate data from each provided data point.
#' @keywords utils
#' @export
get_rastervalue <- function(coords, env, rm.na = FALSE){
  assertthat::assert_that(
    inherits(coords,"sf") || inherits(coords, "Spatial") || (is.data.frame(coords) || is.matrix(coords)),
    is.Raster(env),
    is.logical(rm.na)
    )

  # Try an extraction
  try({ex <- raster::extract(x = env,
                             y = coords,
                             method = "simple",
                             df = TRUE)},silent = FALSE)
  if(class(ex) == "try-error") stop(paste("Raster extraction failed: ", ex))
  # Find those that have NA in there
  check_again <- apply(ex, 1, function(x) anyNA(x))
  if(any(check_again)){
    # Re-extract but with a small buffer
    coords_sub <- coords[which(check_again),]
    try({ex_sub <- raster::extract(x = env,
                               y = coords_sub,
                               method = "simple",
                               small = TRUE,
                               df = TRUE)},silent = FALSE)
    if(class(ex_sub) == "try-error") stop(paste("Raster extraction failed: ", ex_sub))
    ex[which(check_again),] <- ex_sub
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
    ex <- subset(ex, complete.cases(ex))
  }
  assertthat::assert_that(is.data.frame(ex),
                          nrow(ex)>0,
                          msg = "Something went wrong with the extraction or all points had missing data.")
  return(ex)
}

#' Spatial adjustment of environmental predictors and raster stacks
#'
#' @description
#' This function allows the transformation of provided environmental predictors (in [`Raster`] format).
#' A common use case is for instance the standardization (or scaling) of all predictors prior to model fitting.
#' This function works both with [`Raster`] as well as with [`stars`] objects.
#' @details
#' Available options are:
#' * \code{'none'} The original layer(s) are returned.
#' * \code{'scale'} This run the [`scale()`] function with default settings (1 Standard deviation) across all predictors.
#' A sensible default to for most model fitting.
#' * \code{'norm'} This normalizes all predictors to a range from \code{0-1}.
#' * \code{'pca'} This option runs a principal component decomposition of all predictors (via [`prcomp()`]).
#' It returns new predictors resembling all components in order of the most important ones. Can be useful to
#' reduce collinearity, however note that this changes all predictor names to 'PCX', where X is the number of the component.
#' The parameter \code{'pca.var'} can be modified to specify the minimum variance to be covered by the axes.
#' * \code{'revjack'} Removes outliers from the supplied stack via a reverse jackknife procedure.
#' Identified outliers are by default set to \code{NA}.
#' @param env A [`Raster`] object.
#' @param option A [`vector`] stating whether predictors should be preprocessed in any way (Options: \code{'none'},
#' \code{'scale'}, \code{'norm'}, \code{'pca'}, \code{'revjack'}). See Details.
#' @param windsor_props A [`numeric`] vector specifying the proportions to be clipped for windsorization (Default: \code{c(.05,.95)}).
#' @param pca.var A [`numeric`] value between \code{>0} and \code{1} stating the minimum amount of variance to be covered (Default: \code{0.8}).
#' @returns Returns a adjusted [`Raster`] object of identical resolution.
#' @seealso predictor_derivate
#' @examples
#' \dontrun{
#' # Where x is a rasterstack
#' new_x <- predictor_transform(x, option = 'scale')
#' }
#' @keywords utils
#' @export
predictor_transform <- function(env, option, windsor_props = c(.05,.95), pca.var = 0.8, ...){
   assertthat::assert_that(
     inherits(env,'Raster') || inherits(env, 'stars'),
     is.character(option),
     base::length(option) == 1,   # TODO: incorporate possibility of doing multiple options at once and in the order they are supplied?
     is.numeric(windsor_props) & length(windsor_props)==2,
     is.numeric(pca.var)
   )
  # Match option
  option <- match.arg(option, c('none','pca', 'scale', 'norm','windsor', 'revjack'), several.ok = FALSE)

  # Nothing to be done
  if(option == 'none') return(env)

  # If stars see if we can convert it to a stack
  if(inherits(env, 'stars')){
    lyrs <- names(env) # Names of predictors
    times <- stars::st_get_dimension_values(env, which = 3) # Time attribute
    # Convert to list
    env_list <- list()
    for(name in lyrs) env_list[[name]] <- as(env[name], 'Raster')
  }

  # Normalization
  if(option == 'norm'){
    if(is.Raster(env)){
      out <- (env - raster::cellStats(env, stat="min")) /
        (raster::cellStats(env, stat="max") -
           raster::cellStats(env, stat="min"))
    } else {
      out <- lapply(env_list, function(x) {
        (x - raster::cellStats(x, stat="min")) /
          (raster::cellStats(x, stat="max") -
             raster::cellStats(x, stat="min"))
      })
    }
  }
  # Scaling
  if(option == 'scale'){
    if(is.Raster(env)){
      out <- raster::scale(env, center = TRUE, scale = TRUE)
    } else {
      out <- lapply(env_list, function(x) raster::scale(x, center = TRUE, scale = TRUE))
    }
  }

  # Windsorization
  if(option == 'windsor'){
    win <- function(x, windsor_props){
      xq <- stats::quantile(x = x[], probs = windsor_props, na.rm = TRUE)
      min.value <- xq[1]
      max.value <- xq[2]
      if(is.vector(env)) out <- units::drop_units(env) else out <- env
      out[out < min.value] <- min.value
      out[out > max.value] <- max.value
      out
    }
    if(is.Raster(env)){
      out <- win(env, windsor_props )
    } else {
      out <- lapply(env_list, function(x) win(x, windsor_props))
    }
  }

  # Reverse jackknife removal of outliers
  if(option == 'revjack'){
    rj <- function(x){
      o <- emptyraster(x)
      o[] <- rm_outlier_revjack(x[], procedure = "missing")
      return(o)
    }
    if(is.Raster(env)){
      out <- env
      for(n in 1:nlayers(out)){
        out <- raster::addLayer(out, rj(env[[n]]) )
      }
    } else {
      out <- lapply(env_list, function(x) rj(x))
    }
  }

  # Principle component separation of variables
  # Inspiration taken from RSToolbox package
  if(option == 'pca'){
    if(is.Raster(env)){
      assertthat::assert_that(raster::nlayers(env)>=2,msg = 'Need at least two predictors to calculate PCA.')

      # FIXME: Allow a reduction to few components than nr of layers?
      nComp <- nlayers(env)
      # Construct mask of all cells
      envMask <- !sum(raster::calc(env, is.na))
      assertthat::assert_that(cellStats(envMask, sum)>0,msg = 'A predictor is either NA only or no valid values across all layers')
      env <- raster::mask(env, envMask, maskvalue = 0)

      # Sample covariance from stack and fit PCA
      covMat <- raster::layerStats(env, stat = "cov", na.rm = TRUE)
      pca <- stats::princomp(covmat = covMat[[1]], cor = FALSE)
      # Add means and grid cells
      pca$center <- covMat$mean
      pca$n.obs <- raster::ncell(env)

      # Check how many components are requested:
      if(pca.var<1){
        sums <- loadings( summary(pca) )[]
        props <- cumsum(colSums(sums^2) / nrow(sums)) # Cumulative explained variance
        nComp <- length( which(props <= pca.var) )
      }
      # Predict principle components
      out <- raster::predict(env, pca,na.rm = TRUE, index = 1:nComp)
      names(out) <- paste0("PC", 1:nComp)

      return(out)
    } else {
      # TODO:
      stop("Principal component transformation for stars objects is not yet implemented. Pre-process externally!")
    }
  }

  # If stars convert back to stars object
  if(inherits(env, 'stars')){
    # Convert list back to stars
    out <- do.call(
      stars:::c.stars,
      lapply(out, function(x) stars::st_as_stars(x))
    )
    # Reset names of attributes
    names(out) <- lyrs
    out <- stars::st_set_dimensions(out, which = 3, values = times, names = "time")
  } else {
    # Final security checks
    assertthat::assert_that(
      raster::nlayers(env) == raster::nlayers(out),
      is_comparable_raster(out, env)
    )
    return(out)
  }
}

#' Create spatial derivative of raster stacks
#'
#' @description
#' This function creates derivatives of existing covariates and returns them in Raster format.
#' Derivative variables can in the machine learning literature commonly be understood as one aspect of feature
#' engineering. They can be particularly powerful in introducing non-linearities in otherwise linear models,
#' for example is often done in the popular Maxent framework.
#' @details
#' Available options are:
#' * \code{'none'} The original layer(s) are returned.
#' * \code{'quadratic'} A quadratic transformation (\eqn{x^{2}}) is created of the provided layers.
#' * \code{'hinge'} Creates hinge transformation of covariates, which set all values lower than a set threshold to \code{0}
#' and all others to a range of \eqn{[0,1]}. The number of thresholds and thus new derivates is specified
#' via the parameter \code{'nknots'} (Default: \code{4}).
#' * \code{'thresh'} A threshold transformation of covariates, which sets all values lower than a set threshold ot
#' \code{0} and those larger to \code{1}.
#' The number of thresholds and thus new derivates is specified via the parameter \code{'nknots'} (Default: \code{4}).
#' * \code{'bin'} Creates a factor representation of a covariates by cutting the range of covariates by their percentiles.
#' The number of percentile cuts and thus new derivates is specified via the parameter \code{'nknots'} (Default: \code{4}).
#' @param env A [`Raster`] object.
#' @param option A [`vector`] stating whether predictors should be preprocessed in any way
#' (Options: 'none','quadratic', 'hinge', 'thresh', 'bin').
#' @param nknots The number of knots to be used for the transformation (Default: \code{4}).
#' @param deriv A [`vector`] with [`characters`] of specific derivates to create (Default: \code{NULL}).
#' @return Returns the derived adjusted [`Raster`] objects of identical resolution.
#' @seealso predictor_derivate
#' @examples
#' \dontrun{
#' # Create a hinge transformation of one or multiple RasterLayers.
#' predictor_derivate(covs, option = "hinge", knots = 4)
#' }
#' @keywords utils
#' @export
predictor_derivate <- function(env, option, nknots = 4, deriv = NULL, ...){
  assertthat::assert_that(
    inherits(env,'Raster') || inherits(env, "stars"),
    !missing(env),
    is.numeric(nknots) && nknots > 0,
    is.character(option),
    is.null(deriv) || is.character(deriv)
  )
  # Match argument.
  option <- match.arg(option, c('none','quadratic', 'hinge', 'thresh', 'bin'), several.ok = FALSE)

  # None, return as is
  if(option == 'none') return(env)

  # If stars see if we can convert it to a stack
  if(inherits(env, 'stars')){
    assertthat::assert_that(!is.null(deriv),msg = "Derivate names could not be found!")
    # Decompose derivate variable names if set
    deriv <- grep(paste0(option, "__"), deriv, value = TRUE)
    if(length(deriv)==0){
      if(getOption('ibis.setupmessages')) myLog('[Setup]','red','Predictors with derivates not found!')
      return(NULL)
    }
    cutoffs <- do.call(rbind,strsplit(deriv, "__")) |> as.data.frame()
    cutoffs$deriv <- deriv

    lyrs <- names(env) # Names of predictors
    times <- stars::st_get_dimension_values(env, which = 3) # Time attribute
    # Create a list to house the results
    env_list <- list()
    for(name in cutoffs$deriv){
      env_list[[name]] <- as(env[cutoffs[which(cutoffs$deriv==name),2]], 'Raster') # Specify original raster
    }
    assertthat::assert_that(length(env_list) > 0)
  } else {cutoffs <- NULL}

  # Simple quadratic transformation
  if(option == 'quadratic'){
    if(is.Raster(env)){
      new_env <- raster::calc(env, function(x) I(x^2))
      names(new_env) <- paste0('quad__', names(env))
    } else {
      # Stars processing
      new_env <- lapply(env_list, function(x) {
        raster::calc(x, function(z) I(z^2))
      })
    }
  }

  # Hinge transformation
  # From`maxnet` package
  if(option == 'hinge'){
    if(is.Raster(env)){
      # Build new stacks
      new_env <- raster::stack()
      for(val in names(env)){
        o <- makeHinge(env[[val]], n = val, nknots = nknots, cutoffs = cutoffs)
        if(is.null(o)) next()
        new_env <- raster::addLayer(new_env,
                                    fill_rasters(o, emptyraster(env) )
        )
        rm(o)
      }
    } else {
      # Stars object
      for(val in names(env_list)){
        # Format cutoffs
        cu <- cutoffs[which(cutoffs$deriv == val), 3]
        cu <- strsplit(cu, "_") |> unlist()
        # Remove any leading points
        if(any(substr(cu,1, 1)==".")){
          cu[which(substr(cu,1, 1)==".")] <- gsub("^.","",cu[which(substr(cu,1, 1)==".")])
        }
        cu <- as.numeric(cu)
        assertthat::assert_that(!anyNA(cu), is.numeric(cu))
        for(k in 1:nlayers(env_list[[val]])){
          o <- emptyraster(env_list[[val]][[k]])
          o[] <- hingeval(env_list[[val]][[k]][], cu[1], cu[2])
          env_list[[val]][[k]] <- o
          rm(o)
        }
      }
      invisible(gc())
    }
  }

  # For thresholds
  # Take functionality in maxnet package
  if(option == 'thresh'){
    if(is.Raster(env)){
      new_env <- raster::stack()
      for(val in names(env)){
        o <- makeThresh(env[[val]],n = val,nknots = nknots, cutoffs = cutoffs)
        if(is.null(o)) next()
        new_env <- raster::addLayer(new_env,
                                    fill_rasters(o, emptyraster(env))
        )
        rm(o)
      }
    } else {
      # For stats layers
      for(val in names(env_list)){
        # Format cutoffs
        cu <- cutoffs[which(cutoffs$deriv == val), 3]
        cu <- strsplit(cu, "_") |> unlist()
        # Remove any leading points
        if(any(substr(cu,1, 1)==".")){
          cu[which(substr(cu,1, 1)==".")] <- gsub("^.","",cu[which(substr(cu,1, 1)==".")])
        }
        cu <- as.numeric(cu)
        assertthat::assert_that(!anyNA(cu), is.numeric(cu))
        for(k in 1:nlayers(env_list[[val]])){
          o <- emptyraster(env_list[[val]][[k]])
          o[] <- thresholdval(env_list[[val]][[k]][], cu)
          env_list[[val]][[k]] <- o
          rm(o)
        }
      }
      invisible(gc())
    }
  }

  # For binning, calculate cuts of thresholds
  if(option == 'bin'){
    if(is.Raster(env)){
      new_env <- raster::stack()
      for(val in names(env)){
        o <- makeBin(env[[val]],n = val,nknots = nknots, cutoffs = cutoffs)
        if(is.null(o)) next()
        new_env <- raster::addLayer(new_env, o)
        rm(o)
      }
    } else {
      # For stats layers
      for(val in names(env_list)){
        # Format cutoffs
        cu <- cutoffs[which(cutoffs$deriv == val), 3]
        cu <- strsplit(cu, "_") |> unlist()
        # Remove any leading points
        if(any(substr(cu,1, 1)==".")){
          cu[which(substr(cu,1, 1)==".")] <- gsub("^.","",cu[which(substr(cu,1, 1)==".")])
        }
        cu <- as.numeric(cu)
        assertthat::assert_that(!anyNA(cu), is.numeric(cu))
        for(k in 1:nlayers(env_list[[val]])){
          o <- emptyraster(env_list[[val]][[k]])
          o <- raster::cut(env_list[[val]][[k]], cu)
          o[is.na(o)] <- 0
          o <- raster::mask(o, env_list[[val]][[k]] )
          env_list[[val]][[k]] <- o
          rm(o)
        }
      }
      invisible(gc())
    }
  }

  # If stars convert back to stars object
  if(inherits(env, 'stars')){
    # Add the original layers back
    for(name in names(env)){
      env_list[[name]] <- as(env[name], 'Raster') # Specify original raster
    }

    # Convert list back to stars
    new_env <- do.call(
      stars:::c.stars,
      lapply(env_list, function(x) stars::st_as_stars(x))
    )
    # Reset names of attributes
    names(new_env) <- c( cutoffs$deriv, names(env))
    new_env <- stars::st_set_dimensions(new_env, which = 3, values = times, names = "time")
  }
  return(new_env)
}


#' Hinge transformation of a given predictor
#'
#' @description
#' This function transforms a provided predictor variable with a hinge transformation,
#' e.g. a new range of values where any values lower than a certain knot are set to 0,
#' while the remainder is left at the original values.
#' @param v A [`Raster`] object.
#' @param n A [`character`] describing the name of the variable. Used as basis for new names.
#' @param nknots The number of knots to be used for the transformation (Default: \code{4}).
#' @param cutoffs A [`numeric`] vector of optionally used cutoffs to be used instead (Default: \code{NULL}).
#' @keywords utils, internal
#' @concept Concept taken from the [maxnet] package.
#' @returns A hinge transformed [`data.frame`].
#' @noRd
makeHinge <- function(v, n, nknots = 4, cutoffs = NULL){
  assertthat::assert_that(is.Raster(v),
                          is.character(n),
                          is.numeric(nknots),
                          is.numeric(cutoffs) || is.null(cutoffs))
  # Get stats
  v.min <- raster::cellStats(v, min)
  v.max <- raster::cellStats(v, max)
  if(is.null(cutoffs)){
    k <- seq(v.min, v.max, length = nknots)
  } else {
    k <- cutoffs
  }
  if(length(k)<=1) return(NULL)

  # Hinge up to max
  lh <- outer(v[], utils::head(k, -1), function(w, h) hingeval(w,h, v.max))
  # Hinge starting from min
  rh <- outer(v[], k[-1], function(w, h) hingeval(w, v.min, h))
  colnames(lh) <- paste0("hinge__",n,'__', round( utils::head(k, -1), 2),'_', round(v.max, 2))
  colnames(rh) <- paste0("hinge__",n,'__', round( v.min, 2),'_', round(k[-1], 2))
  o <- as.data.frame(
    cbind(lh, rh)
  )
  # Kick out first (min) and last (max) col as those are perfectly correlated
  o <- o[,-c(1,ncol(o))]
  attr(o, "deriv.hinge") <- k
  return(o)
}

#' Threshold transformation of a given predictor
#'
#' @description
#' This function transforms a provided predictor variable with a threshold transformation,
#' e.g. a new range of values where any values lower than a certain knot are set to 0,
#' while the remainder is set to 1.
#' @param v A [`Raster`] object.
#' @param n A [`character`] describing the name of the variable. Used as basis for new names.
#' @param nknots The number of knots to be used for the transformation (Default: \code{4}).
#' @param cutoffs A [`numeric`] vector of optionally used cutoffs to be used instead (Default: \code{NULL}).
#' @keywords utils, internal
#' @concept Concept taken from the [maxnet] package.
#' @returns A threshold transformed [`data.frame`].
#' @noRd
makeThresh <- function(v, n, nknots = 4, cutoffs = NULL){
  assertthat::assert_that(is.Raster(v),
                          is.character(n),
                          is.numeric(nknots),
                          is.numeric(cutoffs) || is.null(cutoffs))
  if(is.null(cutoffs)){
    # Get min max
    v.min <- raster::cellStats(v,min)
    v.max <- raster::cellStats(v,max)
    k <- seq(v.min, v.max, length = nknots + 2)[2:nknots + 1]
  } else {
    k <- cutoffs
  }
  if(length(k)<=1) return(NULL)
  f <- outer(v[], k, function(w, t) ifelse(w >= t, 1, 0))
  colnames(f) <- paste0("thresh__", n, "__",  round(k, 2))
  f <- as.data.frame(f)
  attr(f, "deriv.thresh") <- k
  return(f)
}

#' Binned transformation of a given predictor
#'
#' @description
#' This function takes predictor values and 'bins' them into categories based on a
#' percentile split.
#' @param v A [`Raster`] object.
#' @param n A [`character`] describing the name of the variable. Used as basis for new names.
#' @param nknots The number of knots to be used for the transformation (Default: \code{4}).
#' @param cutoffs A [`numeric`] vector of optionally used cutoffs to be used instead (Default: \code{NULL}).
#' @keywords utils, internal
#' @returns A binned transformed [`data.frame`] with columns representing each bin.
#' @noRd
makeBin <- function(v, n, nknots, cutoffs = NULL){
  assertthat::assert_that(is.Raster(v),
                          is.character(n),
                          is.numeric(nknots),
                          is.numeric(cutoffs) || is.null(cutoffs))
  if(is.null(cutoffs)){
    # Calculate cuts
    cu <- raster::quantile(v, probs = seq(0, 1, by = 1/nknots) )
  } else { cu <- cutoffs}

  if(anyDuplicated(cu)){
    # If duplicated quantiles (e.g. 0, 0, 0.2..), sample from a larger number
    cu <- raster::quantile(v, probs = seq(0, 1, by = 1/(nknots*2)) )
    cu <- cu[-which(duplicated(cu))] # Remove duplicated cuts
    if(length(cu)<=2) return( NULL )
    if(length(cu) > nknots){
      cu <- cu[(length(cu)-(nknots)):length(cu)]
    }
  }
  # Make cuts and explode
  out <- explode_factorized_raster(
    raster::ratify(
      raster::cut(v, cu)
    )
  )
  # Format threshold names
  cu.brk <- as.character(cut(cu[-1], cu))
  cu.brk <- gsub(",","_",cu.brk)
  cu.brk <- gsub("\\(|\\]", "", cu.brk)
  # names(out) <- paste0("bin__",n, "__", gsub(x = names(cu)[-1], pattern = "\\D", replacement = ""),"__", cu.brk )
  names(out) <- paste0("bin__",n, "__", cu.brk )
  for(i in 1:nlayers(out)){
    attr(out[[i]], "deriv.bin") <- cu[i:(i+1)]
  }
  return(out)
}

#' Homogenize NA values across a set of predictors.
#'
#' @description This method allows the homogenization of missing data across a set of environmental predictors.
#' It is by default called when predictors are added to [´BiodiversityDistribution´] object. Only grid cells with NAs that contain
#' values at some raster layers are homogenized.
#' Additional parameters allow instead of homogenization to fill the missing data with neighbouring values
#' @param env A [`Raster`] object with the predictors
#' @param fill A [`logical`] value indicating whether missing data are to be filled (Default: FALSE).
#' @param fill_method A [`character`] of the method for filling gaps to be used (Default: 'ngb')
#' @param return_na_cells A [`logical`] value of whether the ids of grid cells with NA values is to be returned instead (Default: FALSE)
#' @returns A [`Raster`] object with the same number of layers as the input.
#' @keywords utils
#' @noRd
predictor_homogenize_na <- function(env, fill = FALSE, fill_method = 'ngb', return_na_cells = FALSE){
  assertthat::assert_that(
    is.Raster(env) || inherits(env, 'stars'),
    is.logical(fill),
    is.character(fill_method), fill_method %in% c('ngb'),
    is.logical(return_na_cells)
  )
  # Workflow for raster layers
  if(is.Raster(env)){
    nl <- raster::nlayers(env)
    # If the number of layers is 1, no need for homogenization
    if(nl > 1){
      # Calculate number of NA grid cells per stack
      mask_na <- sum( is.na(env) )
      # Remove grid cells that are equal to the number of layers (all values NA)
      none_area <- mask_na == nl
      none_area[none_area == 0 ] <- NA
      mask_na <- raster::mask(mask_na,
                              mask = none_area,inverse = TRUE)

      # Should any fill be conducted?
      if(fill){
        stop('Not yet implemented!')
      } else {
        # Otherwise just homogenize NA values across predictors
        if(cellStats(mask_na,'max')>0){
          mask_all <- mask_na == 0; mask_all[mask_all == 0] <- NA
          env <- raster::mask(env, mask = mask_all)
        }
      }
      # Should NA coordinates of cells where 1 or more predictor is NA be returned?
      # FIXME: One could directly return a data.frame with the predictor names to allow easier lookup.
      if(return_na_cells){
        vals <- which((mask_na>0)[])
        env <- list(cells_na = vals, env = env)
      }
      rm(mask_na, none_area) # Cleanup
    }
  } else if(inherits(env, 'stars')){
    stop('Not implemented yet.')
  }
  # Security checks
  assertthat::assert_that(
    is.Raster(env) || is.list(env) || inherits(env, 'stars')
  )
  # Return the result
  return(env)
}

#' Create new raster stack from a given data.frame
#'
#' @param post A data.frame
#' @param background A [`Raster-class`] object for the background raster
#' @keywords internal, utils
#' @return A [`Raster-class`] object with number of columns equal to ncol(post)
#' @noRd
fill_rasters <- function(post, background){
  assertthat::assert_that(
    is.data.frame(post),ncol(post)>1,
    inherits(background,'Raster'),
    nrow(post) == ncell(background)
  )
  # Make names to be sure
  names(post) <- base::make.names(names(post))

  # If only one raster
  if(ncol(post)==1){
    out <- emptyraster(background)
    out[] <- post[,1]
  } else {
    # Loop through each column
    out <- raster::stack()
    for(co in 1:ncol(post)){
      o <- emptyraster(background)
      o[] <- post[,co] # Assign values
      # Add to stack
      out <- raster::addLayer(out, o)
      rm(o)
    }
  }
  # Assign names
  names(out) <- names(post)

  # Check that derivate attributes if existing are passed
  if(length( grep("deriv", names(attributes(post)) ))>0){
    attr(out, grep("deriv", names(attributes(post)),value = TRUE) ) <- attr(post, grep("deriv", names(attributes(post)),value = TRUE) )
  }

  # Final check
  assertthat::assert_that(
    inherits(out,'Raster'),
    nlayers(out) == ncol(post)
  )
  return(out)
}

#' Create a polynomial transformation from coordinates
#'
#' @description This function transforms the coordinates of a supplied file through a polynomial transform.
#' By default it applies weights and a QR decomposition for numerical stability.
#' @param coords A [`data.frame`], [`matrix`] or [`sf`] object with coordinates (2 columns named x-y)
#' @param degree The number of degrees used for polynominal transformation (Default: 2)
#' @param weights Set by default to the inverse of the number of coordinates.
#' @returns A data.frame with transformed coordinates.
#' @keywords utils
#' @keywords internal
#' @references Dray S., Plissier R., Couteron P., Fortin M.J., Legendre P., Peres-Neto P.R., Bellier E., Bivand R., Blanchet F.G., De Caceres M., Dufour A.B., Heegaard E., Jombart T., Munoz F., Oksanen J., Thioulouse J., Wagner H.H. (2012). Community ecology in the age of multivariate multiscale spatial analysis. Ecological Monographs 82, 257–275.
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
  a0 <- poly(x = as.matrix( coords ), degree = degree, simple = TRUE)
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
#' Completely deletes for instance a temporary created raster file.
#' @param A [`raster`] object.
#' @param verbose Print progress (Default: \code{FALSE})
#' @keywords utils
clean_rasterfile <- function(x, verbose = FALSE)
{
  stopifnot(grepl("Raster", class(x)))
  if (!fromDisk(x))
    return(NULL)
  sink(tempfile())
  tdir = rasterOptions()[["tmpdir"]]
  sink(NULL)
  if (class(x) == "RasterLayer")
    files = basename(x@file@name)
  if (class(x) == "RasterStack")
    files = do.call(c, lapply(methods::slot(x, "layers"),
                              function(x) x@file@name))
  files = files[file.exists(files)]
  if (length(files) == 0)
    return(NULL)
  lapply(files, function(f) {
    if (fromDisk(x) & file.exists(f))
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
#'  a new [`RasterStack`] that contains the individual levels.
#'  @param ras A [`RasterLayer`] object that is a [`factor`]. Alternatively a [`RasterStack`] object
#'  can be supplied in which only factor variables are 'exploded'
#'  @param name An optional [`character`] name for the [raster].
#'  @param ... Other parameters (not used).
#'  @returns A [`RasterStack`] object
#'  @keywords utils
#'  @noRd
explode_factorized_raster <- function(ras, name = NULL, ...){
  assertthat::assert_that(is.Raster(ras),
                          is.null(name) || is.character(name))

  # Simply return the input if there are no factors
  if(!any(is.factor(ras))) return(ras)

  # If input is a RasterLayer
  if(inherits(ras, 'RasterLayer')){
    # Get name
    # Create output template
    temp <- emptyraster(ras)
    if(is.null(name)) name <- names(ras)

    # Extract data
    o <- data.frame(val = values(ras));names(o) <- name;o[[name]] <- factor(o[[name]])

    # Make function that converts all factors to split rasters
    f <- as.data.frame(
      outer(o[[name]], levels(o[[name]]), function(w, f) ifelse(w == f, 1, 0))
    )

    # Fill template rasters
    out <- fill_rasters(f,temp)
    names(out) <- paste(name, levels(o[[name]]), sep = ".")

  } else if(inherits(ras, 'RasterStack') || inherits(ras, 'RasterBrick')){
    # Alternatively if input is stack
    fcts <- is.factor(ras)

    # Get non-factor variables
    out <- ras[[which(!fcts)]]
    for(k in which(fcts)){

      sub <- ras[[k]]

      temp <- emptyraster(sub)
      if(is.null(name)) new_name <- names(sub)

      # Extract data
      o <- data.frame(val = values(sub));names(o) <- new_name;o[[new_name]] <- factor(o[[new_name]])

      # Make function that converts all factors to split rasters
      f <- as.data.frame(
        outer(o[[new_name]], levels(o[[new_name]]), function(w, f) ifelse(w == f, 1, 0))
      )

      # Fill template rasters
      new <- fill_rasters(f, temp)
      names(new) <- paste(new_name, levels(o[[new_name]]), sep = ".")
      out <- raster::addLayer(out, new)
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
#'  [*] \code{"random"}: Samples at random up to number of \code{"minpoints"} across all occupied grid cells.
#'  Does not account for any spatial or environmental distance between observations.
#'  [*] \code{"bias"}: As for random, but here occurrence points are incrementally removed from the dataset, starting
#'  with the ones that are most biased (parameter \code{"bias"}) first. Thins the observations up to \code{"minpoints"}.
#'  [*] \code{"zones"}: Assesses for each observation that it falls with a maximum of \code{"minpoints"} into
#'  each occupied zone. Careful: If the zones are relatively wide this can remove quite a few observations.
#'  [*] \code{"spatial"}: Calculates the spatial distance between all observations. Then points are removed
#'  iteratively until the minimum distance between points is crossed.  The \code{"mindistance"} parameter has to
#'  be set for this function to work.
#'  [*] \code{"environmental"}: This approach creates an observation-wide clustering under the assumption that the
#'  environmental niche has been comprehensively sampled. We then obtain an equal number (\code{"minpoints"}) of
#'  observations  for clustered zone. This ensures that environmental conditions are approximately covered.
#'
#' @param df A [`sf`] or [`data.frame`] object with observed occurrence points. All methods threat presence-only
#' and presence-absence occurrence points equally.
#' @param background A [`RasterLayer`] object with the background of the study region. Use for assessing point density.
#' @param env A [`Raster`] object with environmental covaraites. Needed when method is set to \code{"environmental"} (Default: \code{NULL}).
#' @param method A [`character`] of the method to be applied (Default: \code{"random"}).
#' @param minpoints A [`numeric`] giving the number of data points at minimum to take (Default: \code{10}).
#' @param mindistance A [`numeric`] for the minimum distance of neighbouring observations (Default: \code{NULL}).
#' @param bias A supplied [`RasterLayer`] bias layer for method \code{"bias"} (Default: \code{NULL}).
#' @param zones A [`RasterLayer`] to be supplied when option \code{"method"} is chosen (Default: \code{NULL}).
#' @param verbose [`logical`] of whether to print some statistics about the thinning outcome (Default: \code{TRUE}).
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
  method <- match.arg(method, choices = c("random", "spatial", "environmental", "zones"), several.ok = FALSE)

  # Label background with id
  bg <- background
  bg[] <- 1:raster::ncell(bg)
  bg <- raster::mask(bg, background)

  # Take coordinates of supplied data and rasterize
  coords <- sf::st_coordinates( df )
  ras <- raster::rasterize(coords, bg) # Get the number of observations per grid cell

  # Bounds for thining
  totake <- c(lower = minpoints, upper = max(raster::cellStats(ras,"min"), minpoints))

  # -- #
  if(method == "random"){
    # For each unique grid cell id, get the minimum value up to a maximum of the points
    # by sampling at random from the occupied grid cells

    # Output vector
    sel <- vector()

    ex <- data.frame(id = 1:nrow(coords),
                     cid = raster::extract(bg, coords)
    )
    ex <- subset(ex, complete.cases(ex)) # Don't need missing points

    ex <- dplyr::left_join(ex,
                           ex %>% dplyr::group_by(cid) %>% dplyr::summarise(N = dplyr::n()),
                           by = "cid"
    )
    # Points to take
    sel <- append(sel, ex$id[which(ex$N <= min(totake))] )

    # For those where we have more than the minimum, take at random the upper limits of observations
    ex$oversampled <- ifelse(ex$N >= totake["upper"], 1, 0)
    # Now sample at random up to the maximum amount. Got tired of doing this outside tidyverse
    o <- ex %>% dplyr::filter(oversampled == 1) %>%
      dplyr::group_by(cid) %>%
      dplyr::slice_sample(n = min(totake))
    if(nrow(o)>0) sel <- append(sel, o$id)
    if(anyDuplicated(sel)) sel <- unique(sel)
    rm(o,ex)
  } else if(method == "environmental"){
    # Environmental clustering

  } else if(method == "spatial"){
    # Spatial thinning

  } else if(method == "zones"){
    # Thinning by zones
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
