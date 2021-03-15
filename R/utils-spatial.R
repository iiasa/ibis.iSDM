#' Are rasters comparable?
#'
#' This function checks if two `Raster-class` objects
#' are comparable.
#'
#' @param x [`Raster-class`] object.
#' @param y [`Raster-class`] object.
#' @keywords internal
#' @return [`logical`] indicating if the two [`Raster-class`] objects have the same
#'   resolution, extent, dimensionality, and coordinate system.
#'
#' @noRd
is_comparable_raster <- function(x, y) {
  assertthat::assert_that(inherits(x, "Raster"), inherits(y, "Raster")) &&
    sf::st_crs(x@crs) == sf::st_crs(y@crs) &&
    raster::compareRaster(x, y, crs = FALSE, res = TRUE, tolerance = 1e-5,
                          stopiffalse = FALSE)
}

assertthat::on_failure(is_comparable_raster) <- function(call, env) {
  paste0(deparse(call$x), " and ", deparse(call$y),  " are not comparable: ",
         "they have different spatial resolutions, extents, ",
         "coordinate reference systems, or dimensionality (rows / columns)")
}

#' Do extents intersect?
#'
#' Verify if the extents of two spatial objects intersect or not.
#'
#' @param x [`Raster-class`], [`Spatial-class`] or [`sf::sf()`] object.
#' @param y [`Raster-class`], [`Spatial-class`] or [`sf::sf()`] object.
#' @keywords internal
#' @return [`logical`].
#' @noRd
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
#' @param poly A [`sf`] object
#' @param points A [`Spatial`] or [`sf`] object
point_in_polygon <- function(poly, points){
  assertthat::assert_that(
    inherits(poly,'sf'),
    inherits(points,'sf') || inherits(points,'Spatial') || inherits(points,'data.frame')
  )
  # Convert to sf
  points <- sf::st_as_sf(points) %>%
    # Convert to be sure
    st_transform(crs = st_crs(poly))
  assertthat::assert_that(
    st_crs(poly) == st_crs(points)
  )

  # Within test
  ov <- sf::st_join(points, poly, join = st_within)
  return(ov)
}

#' Converts a bounding box to a Well Known Text polygon
#'
#' @param minx Minimum x value, or the most western longitude
#' @param miny Minimum y value, or the most southern latitude
#' @param maxx Maximum x value, or the most eastern longitude
#' @param maxy Maximum y value, or the most northern latitude
#' @param all A [`vector`] of length 4, with the elements: minx, miny, maxx, maxy
#' @return An object of class [`character`], a Well Known Text (WKT) string of the form
#' 'POLYGON((minx miny, maxx miny, maxx maxy, minx maxy, minx miny))'
#' @keywords internal
#' @noRd

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
#' @param f value to increase the extent (Default = 0.1)
#' @return Returns the unified total extent object
#' @noRd

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

#' Align a [`Raster-class`] object to another by harmonizing geometry and extend.
#'
#' If the data is not in the same projection as the template, the alignment
#' will be computed by reprojection only. If the data has already the same
#' projection, the data set will be cropped and aggregated prior to resampling
#' in order to reduce computation time.
#'
#' Nearest Neighbour resampling (NGB) is recommended for discrete and Bilinear
#' resampling for continious data.
#'
#' @param data [`Raster-class`] object to be resampled
#' @param template [`Raster-class`] or [`Spatial-class`] object from which geometry can be extracted
#' @param method method for resampling ("ngb" or "bilinear")
#' @param func function for resampling (Default: mean)
#' @param cl Boolean value if multicore computation should be used (Default: TRUE)
#' @return New [`Raster`] object aligned to the supplied template layer
#' @noRd

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
#' @description .
#'
#' @param x a \code{raster*} object corresponding.
#' @param ... other arguments that can be passed to \code{\link{raster}}
#' @return an empty raster, i.e. all cells are \code{NA}.
#' @import raster
#' @keywords raster
#' @examples
#' require(raster)
#' r <- raster(matrix(1:100, 5, 20))
#' emptyraster(r)
#' @noRd

emptyraster <- function(x, ...) { # add name, filename,
  assertthat::assert_that(inherits(x,'Raster'))
  raster::raster(nrows=nrow(x), ncols=ncol(x),
                        crs=x@crs,
                        ext=extent(x), ...)
}


#' Function to extract nearest neighbour predictor values of provided points
#'
#'
#' @param coords A [`matrix`], [`data.frame`] or [`sf`] object.
#' @param env A [`data.frame`] object with the predictors
#' @param longlat A boolean variable indicating whether the projection is long-lat
#' @param field_space A [`vector`] highlight the columns from which coordinates
#'  are to be extracted (default: c('X','Y'))
#' @return Extracted data from each point
#' @note If multiple values are of equal distance, average them
#' @noRd

get_ngbvalue <- function(coords, env, longlat = TRUE, field_space = c('X','Y'),...) {

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

  # Pairwise distance function
  # Use geosphere distance calculation if available as faster
  if('geosphere' %in% installed.packages()[,1]){
    disfun <- geosphere::distGeo
  } else { disfun <- sp::spDistsN1 }

  env_sub <- apply(coords, 1, function(xy1, xy2) {
                    dists <- disfun(xy2, xy1, TRUE)
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

  out <- do.call(rbind,env_sub)
  out[,field_space] <- as.data.frame(coords) # Ensure that coordinates are back in
  return(out)
}

#' Spatial adjustment of raster stacks
#'
#' @param env A [`Raster`] object
#' @param option A [`vector`] stating whether predictors should be preprocessed in any way (Options: 'none','scale','norm','pca')
#' @return Returns a adjusted [`Raster`] object of identical resolution
#' @noRd
predictor_transform <- function(env, option, ...){
   assertthat::assert_that(
     inherits(env,'Raster'),
     is.character(option),
     option %in% c('none','pca', 'scale', 'norm')
   )
  # TODO: Another option would be a windsoriation, e.g. cut of extremes
  # TODO: incorporate possipibiltiy of doing multiple options at once

  # Nothing to be done
  if(option == 'none') return(env)

  # Normalization
  if(option == 'norm'){
    out <- (env - raster::cellStats(env, stat="min")) /
      (raster::cellStats(env, stat="max") -
         raster::cellStats(env, stat="min"))
  }
  # Scaling
  if(option == 'scale'){
    out <- raster::scale(env, center = TRUE, scale = TRUE)
  }

  # Principle component separation of variables
  # Inspiration taken from RSToolbox package
  if(option == 'pca'){
    assertthat::assert_that(raster::nlayers(env)>=2,msg = 'Need at least two predictors to calculate PCA.')

    # FIXME: Allow a reduction to few components than nr of layers?
    nComp <- nlayers(env)
    # Construct mask of all cells
    envMask <- !sum(raster::calc(env, is.na))
    assertthat::assert_that(cellStats(envMask,sum)>0,msg = 'A predictor is either NA only or no valid values across all layers')
    env <- raster::mask(env, envMask, maskvalue = 0)

    # Sample covariance from stack and fit PCA
    covMat <- raster::layerStats(env, stat = "cov", na.rm = TRUE)
    pca <- stats::princomp(covmat = covMat[[1]], cor = FALSE)
    # Add means and grid cells
    pca$center <- covMat$mean
    pca$n.obs <- raster::ncell(env)

    # FIXME: Allow parallel processing for rather large files and check how many nodes are available
    # Predict principle components
    out <- raster::predict(env, pca,na.rm = TRUE, index = 1:nComp)

    names(out) <- paste0("PC", 1:nComp)
    return(out)
  }

  # Final security checks
  assertthat::assert_that(
    raster::nlayers(env) == raster::nlayers(out),
    is_comparable_raster(out,env)
  )
  return(out)
}

#' Create spatial derivate of raster stacks
#'
#' @param env A [`Raster`] object
#' @param option A [`vector`] stating whether predictors should be preprocessed in any way (Options: 'none','quadratic', 'hinge', 'thresh')
#' @return Returns the derived adjusted [`Raster`] objects of identical resolution
#' @noRd
predictor_derivate <- function(env, option, ...){
  assertthat::assert_that(
    inherits(env,'Raster'),
    is.character(option),
    option %in% c('none','quadratic', 'hinge', 'thresh')
  )

  # Simple quadratic transformation
  if(option == 'quadratic'){
    new_env <- calc(env, function(x) I(x^2))
    names(new_env) <- paste0('quad.', names(new_env))
  }

  # Hinge transformation
  # From`maxnet` package
  if(option == 'hinge'){
    # Function to create hingevalue
    makeHinge <- function(v, n, nknots = 4){
      # Function to create hingeval
      hingeval <- function (x, min, max) ifelse(is.na(x),NA, pmin(1, pmax(0, (x - min)/(max - min),na.rm = TRUE),na.rm = TRUE))
      # Get stats
      v.min <- raster::cellStats(v, min)
      v.max <- raster::cellStats(v, max)
      k <- seq(v.min, v.max, length = nknots)
      # Hinge up to max
      lh <- outer(v[], utils::head(k, -1), function(w, h) hingeval(w,h, v.max))
      # Hinge starting from min
      rh <- outer(v[], k[-1], function(w, h) hingeval(w, v.min, h))
      colnames(lh) <- paste("hinge",n, round( utils::head(k, -1), 2),'__', round(v.max, 2),sep = ":")
      colnames(rh) <- paste("hinge",n, round( v.min, 2),'__', round(k[-1], 2), sep = ":")
      o <- as.data.frame(
        cbind(lh, rh)
      )
      # Kick out first (min) and last (max) col as those are perfectly correlated
      o <- o[,-c(1,ncol(o))]
      return(o)
    }

    # Build new stacks
    new_env <- raster::stack()
    for(val in names(env)){
      new_env <- raster::addLayer(new_env,
                              fill_rasters(
                                makeHinge(env[[val]],n = val,nknots = 4),
                                           emptyraster(env)
                                )
                             )
    }
  }

  # For thresholds
  # Take functionality in maxnet package
  if(option == 'thresh'){
    # Function to create derivative thresholds
    makeThresh <- function(v,n, nknots = 4){
      # Get min max
      v.min <- cellStats(v,min)
      v.max <- cellStats(v,max)
      k <- seq(v.min, v.max, length = nknots + 2)[2:nknots + 1]
      f <- outer(v[], k, function(w, t) ifelse(w >= t, 1, 0))
      colnames(f) <- paste("thresh", n, round(k, 2), sep = ":")
      return(as.data.frame(f))
    }

    new_env <- raster::stack()
    for(val in names(env)){
      new_env <- raster::addLayer(new_env,
                                  fill_rasters(
                                    makeThresh(env[[val]],n = val,nknots = 4),
                                    emptyraster(env)
                                  )
      )
    }


  }

  return(new_env)
}
#' Create new raster stack from a given data.frame
#'
#' @param post A data.frame
#' @param background A [`Raster-class`] object for the background raster
#' @keywords internal
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

  # Final check
  assertthat::assert_that(
    inherits(out,'Raster'),
    nlayers(out) == ncol(post)
  )
  return(out)
}

#' Clean up raster layer from disk
#'
#' @param A [`raster`] object
#' @noRd
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
