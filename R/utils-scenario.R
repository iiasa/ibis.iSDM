#' Approximate missing time steps between dates
#'
#' @description
#' This function linearly approximates shares between time steps, so that gaps for instance
#' between 2010 and 2020 are filled with data for 2010, 2011, 2012, etc.
#' @param env A [`stars`] object.
#' @param date_interpolation [`character`] on how missing dates between events should be interpolated. See [`project()`].
#' @return [`logical`] indicating if the two [`Raster-class`] objects have the same
#' @keywords scenario
#' @noRd
approximate_gaps <- function(env, date_interpolation = "annual"){
  assertthat::assert_that(
    inherits(env, "stars"),
    is.character(date_interpolation)
  )
  check_package("dplyr")
  date_interpolation <- match.arg(date_interpolation, c("none", "yearly", "annual", "monthly", "daily"), several.ok = FALSE)
  if(date_interpolation=="none") return(env)

  stop("Still in progress")
  # --- #
  # Get individual time steps at interval
  times <- stars::st_get_dimension_values(env, which = names(dim(env))[3], center = TRUE)
  times <- to_POSIXct(times)
  tzone <- attr(as.POSIXlt(times), "tzone")[2] # Get timezone
  assertthat::assert_that(tzone != "", length(times)>=2)
  # Interpolate time steps
  inc <- switch (date_interpolation,
                 "yearly" = "year",
                 "annual" = "year",
                 "monthly" = "month",
                 "daily" = "day"
  )
  new_times <- seq.Date(from = as.Date(times[1],tz = tzone), to = as.Date(times[length(times)],tz = tzone), by = inc)
  new_times <- to_POSIXct(new_times)

  # Linearly approximate all attributes for new object
  # FIXME: Probably terribly memory inefficient but works
  new <- as.data.frame(env)
  assertthat::assert_that(assertthat::has_name(new,c("x","y","time")))
  new <- dplyr::right_join(new, expand.grid(x = unique(new$x), y = unique(new$y), time = new_times),
                          by = c("x", "y","time"))
  # Sort by time
  new <- new[order(new$time),]
  # Now linearly interpolate the missing values per grid cell
  new2 <- apply(new[,4:ncol(new)], 2, function(z){
    # if(inherits(z, "POSIXct")) return(z)
    if(all(is.na(z))) return(z)
    approx(y = z, x = as.numeric(new$time), method = "linear")
  })

  # Steps:
  # empty_stars
  # Join with existing one
  # approxNA

  tt <- as.numeric(new_times)
  # Calc pixel-wise linear slope
  out <- stars::st_apply(
    env,
    1:2,
    function(x) {
      if (anyNA(x))
        NA_real_
      else
        lm.fit(cbind(1, tt), x)$coefficients[2]
    }
  )

  # new <- stars::st_redimension(out, along = list(time = new_times))
}

#' Aggregate stars variables across dimensions
#'
#' @description
#' Small helper function that acts a wrapper to combine
#' 2 or more variables in a `stars` object together.
#' @note Currently only works via matrix manipulation
#' @param obj A [`stars`] object or a [`list`] that can be coerced to one.
#' @param vars A [`vector`] describing the variables to be combined. Has to be of
#' length two or greater.
#' @param newname A [`character`] with the new name for the variable.
#' @param fun A function how the respective layers should be combined.
#' @examples
#' \dontrun{
#'  st_reduce(obj, vars = c('forestShare', 'forestShare.2'),
#'   newname = "forest",fun = "sum")
#' }
#' @keywords scenario, internal
st_reduce <- function(obj, vars, newname, fun = 'sum'){
  assertthat::assert_that(
    is.list(obj) || inherits(obj, 'stars'),
    is.vector(vars) && is.character(vars),
    length(vars) >=2,
    is.character(fun),
    is.character(newname)
  )
  fun <- match.arg(fun, c("sum", "multiply", "divide", "subtract", "mean"), several.ok = FALSE)
  check_package('stars')
  # Convert to stars if not already
  if(!inherits(obj, 'stars')) obj <- stars::st_as_stars(obj)
  # Check that variables are present
  assertthat::assert_that(
    all(vars %in% names(obj))
  )
  # Future?
  if(foreach:::getDoParRegistered()){
    ibis_future(cores = getOption("ibis.nthread"), strategy = getOption("ibis.futurestrategy"))
    fut <- TRUE
  } else { fut <- FALSE }
  # --- #
  # First get all target variables and non-target variables
  target <- stars:::select.stars(obj, vars)
  non_target <- stars:::select.stars(obj, -vars)

  # Now apply the function on the target
  # FIXME: Not working as intended and horribly slow
  # target <- stars::st_apply(target, dims, fun,
  #                           PROGRESS = TRUE,
  #                           FUTURE = fut,
  #                           ...
  #                           )
  what <- switch (fun,
    "sum" = "+",
    "multiply" = "*",
    "divide" = "/",
    "subtract" = "-",
    "mean" = "+"
  )
  new <- stars::st_as_stars( Reduce(what, target) )
  if(what == "mean") new <- new / length(target)
  stars::st_dimensions(new) <- stars::st_dimensions(target)
  # Rename to newname
  names(new) <- newname

  # Combine the result again with non-target
  out <- c(non_target, new)
  rm(new, non_target, target)

  return(
    out
  )
}

#' Converts a stars object to list of rasters
#'
#' @description
#' This is a small helper function to convert a [`stars`] object
#' to a [`Raster`] object. It is possible to select the time frame as well.
#' If multiple \code{"which"} entries are specified, then a [`list`] will be returned.
#' @param obj A [`stars`] object with a \code{"time"} dimension at least.
#' @param which The time entry to use for subsetting. Can be single [`numeric`] or a [`vector`]
#' of numeric time entries corresponding to the time dimension (Default: \code{NULL}).
#' @param template An optional [`Raster`] template to which the output should be aligned too.
#' @returns A [`list`] containing [`Raster`] objects.
#' @keywords scenario, internal
stars_to_raster <- function(obj, which = NULL, template = NULL){
  assertthat::assert_that(
    inherits(obj, 'stars'),
    is.null(which) || is.numeric(which),
    is.null(template) || is.Raster(template)
  )
  assertthat::assert_that(
    length(which) <= dim(obj)['time']
  )
  # Get time dimension and correct if specific entries are requested
  times <- stars::st_get_dimension_values(obj, "time")
  if(is.null(which)) {
    which <- 1:length(times) # Use default length
  }

  # Output type raster
  out <- list()
  for(tt in which){
    # Slice to a specific time frame for each
    o <- obj %>% stars:::slice.stars("time" , tt) |>
      as("Raster")

    # Reset times to the correct ones
    o <- raster::setZ(o, rep(times[tt], raster::nlayers(o)))

    # Now transform the out put if template is set
    if(!is.null(template)){
      if(is.Raster(template)){
        # Check again if necessary to rotate
        if(!raster::compareCRS(o, template)){
          o <- raster::projectRaster(from = o, crs = template, method = "bilinear")
          names(o) <- class_units
        }
        # Now crop and resample to target extent if necessary
        if(!compareRaster(o, template, stopiffalse = FALSE)){
          o <- raster::crop(o, template)
          o2 <- try({alignRasters(data = o,
                                  template = template,
                                  method = "bilinear",
                                  func = "mean", cl = FALSE)
          },silent = TRUE)
          if(inherits(o2,"try-error")){
            o <- raster::resample(o, template,
                                  method = "bilinear")
          } else { o <- o2; rm(o2)}
        }
      }
    } # End of template adjustments
    out[[paste0("time",times[tt])]] <- o
  }
  return( out )
}

#' Quick handy function to calculate the centre of a range
#'
#' @param ras A [`RasterLayer`] object for which the centre of the range is to be calculated.
#' If the distribution is continious, then the centre is calculated as the value centre to all non-NA values.
#' @keywords scenario, internal
#' @noRd
calculate_range_centre <- function(ras) {
  assertthat::assert_that(
    is.Raster(ras)
  )
  r_wt <- area(ras)
  values(r_wt)[is.na(values(ras))] <- NA

  spdf <- rasterToPoints(stack(ras,r_wt), spatial=T)


  # spatialEco::wt.centroid(spdf, sp = F)
  wt.centroid <- function(x, p, sp = TRUE) {
    # if(class(x)[1] == "sf") { x <- as(x, "Spatial") }
    if (!inherits(x, "SpatialPointsDataFrame"))
      stop(deparse(substitute(x)), " MUST BE A SpatialPointsDataFrame OBJECT")
    p <- x@data[, p]
    Xw <- sum(sp::coordinates(x)[, 1] * p)
    Yw <- sum(sp::coordinates(x)[, 2] * p)
    wX <- Xw/sum(p)
    wY <- Yw/sum(p)
    if (sp == FALSE) {
      return(c(wX, wY))
    } else {
      xy <- sp::SpatialPoints(matrix(c(wX, wY), nrow = 1, ncol = 2))
      sp::proj4string(xy) <- sp::proj4string(x)
      return(xy)
    }
  }
}
