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
  # Take name of third band, assuming it to be time
  time_band <- names(dim(obj))[3]

  assertthat::assert_that(
    length(which) <= dim(obj)[time_band]
  )
  # Get time dimension and correct if specific entries are requested
  times <- stars::st_get_dimension_values(obj, time_band, center = TRUE)
  if(is.null(which)) {
    which <- 1:length(times) # Use default length
  }

  # Output type raster
  out <- list()
  for(tt in which){
    # Slice to a specific time frame for each
    o <- obj %>% stars:::slice.stars({{time_band}}, tt) |>
      as("Raster")

    # Reset times to the correct ones
    o <- raster::setZ(o, rep(times[tt], raster::nlayers(o)))

    # Now transform the out put if template is set
    if(!is.null(template)){
      if(is.Raster(template)){
        # Check again if necessary to rotate
        if(!raster::compareCRS(o, template)){
          o <- raster::projectRaster(from = o, crs = template, method = "ngb")
          names(o) <- names(obj)
        }
        # Now crop and resample to target extent if necessary
        if(!compareRaster(o, template, stopiffalse = FALSE)){
          o <- raster::crop(o, template)
          o2 <- try({alignRasters(data = o,
                                  template = template,
                                  method = "ngb",
                                  func = "mean", cl = FALSE)
          },silent = TRUE)
          if(inherits(o2,"try-error")){
            o <- raster::resample(o, template,
                                  method = "ngb")
          } else { o <- o2; rm(o2)}
        }
      }
    } # End of template adjustments
    out[[paste0("time",times[tt])]] <- o
  }
  return( out )
}

#' This function add layers from a RasterStack to a stars object
#'
#' @description
#' Often it is necessary to add static variables to existing stars objects.
#' These will be replicated across the time dimension. This function is a small helper function
#' that allows the addition of said raster stacks to a stars object.
#' @param obj A [`stars`] object with a time dimension (\code{"time"}).
#' @param new A [`RasterStack`] object with additional covariates to be added.
#' @returns A [`stars`] object with the names of the [`Raster`] object added.
#' @keywords scenario, internal
st_add_raster <- function(obj, new){
  assertthat::assert_that(
    inherits(obj, "stars"),
    is.Raster(new),
    raster::nlayers(new) >= 1
  )

  # Check whether there are any variables in the stars object already, if so drop
  if(any(names(new) %in% names(obj))){
    myLog("[Starting]", "yellow", "Duplicate variables in stars and new objects.")
    new <- raster::dropLayer(new, which( names(new) %in% names(obj) ) )
  }

  # Get times objects
  times <- rep(stars::st_get_dimension_values(obj, "time"))
  # Overwrite time dimension
  full_dims <- stars::st_dimensions(obj)

  # Now loop through each layer and add it to the target file
  for(lyr in names(new)){
    s <- raster::stack(replicate(length(times), new[[lyr]])) |>
      stars::st_as_stars()
    names(s) <- lyr

    dims <- stars::st_dimensions(s)
    # Replace the band variable with the original one
    names(dims)[3] <- "time"
    dims$time <- full_dims$time
    stars::st_dimensions(s) <- dims
    obj <- c(obj, s)
  }
  assertthat::assert_that(
    all(all(names(new) %in% names(obj)))
  )
  return(obj)
}

#' Summarize results from scenario projection object
#'
#' @description
#' This is a wrapper function to summarize the output of a scenario projection. The
#' output will contain the average change in the layer per time step.
#' A parameter called \code{"relative"} can be set to calculate relative change instead.
#' @param scenario A [`stars`] object with a time dimension.
#' @param relative A [`logical`] check whether to calculate relative changes instead.
#' @keywords internal, scenario
#' @noRd
summarise_projection <- function(scenario, fun = "mean", relative = TRUE){
  assertthat::assert_that(
    is.list(scenario) || inherits(scenario, "stars"),
    length(dim(scenario))==3,
    is.logical(relative)
  )
  fun <- match.arg(fun, c("mean", "sum"),several.ok = FALSE)

  # Convert to scenarios to data.frame
  df <- stars:::as.data.frame.stars(stars:::st_as_stars(scenario)) %>% subset(., complete.cases(.))
  names(df) <- c("x", "y", "band", "suitability")
  # Add grid cell grouping
  df <- df %>% dplyr::group_by(x,y) %>% dplyr::mutate(id = dplyr::cur_group_id()) %>%
    dplyr::ungroup() %>% dplyr::select(-x,-y) %>%
    dplyr::arrange(id, band)

  # Summarize the overall moments
  if(fun == "mean"){
    out <- df %>%
      dplyr::filter(suitability > 0) %>%
      dplyr::group_by(band) %>%
      dplyr::summarise(suitability_mean = mean(suitability, na.rm = TRUE),
                       suitability_q25 = quantile(suitability, .25),
                       suitability_q50 = quantile(suitability, .5),
                       suitability_q75 = quantile(suitability, .75))
    # Total amount of area lost / gained / stable since previous time step
    totchange_occ <- df %>%
      dplyr::group_by(id) %>%
      dplyr::mutate(change = (suitability - dplyr::lag(suitability)) ) %>% dplyr::ungroup()
    o <- totchange_occ %>% dplyr::group_by(band) %>%
      dplyr::summarise(suitability_avggain = mean(change[change > 0]),
                       suitability_avgloss = mean(change[change < 0]))

    out <- out %>% dplyr::left_join(o, by = "band")
    if(relative){
      # Finally calculate relative change to baseline (first entry) for all entries where this is possible
      relChange <- function(v, fac = 100) (((v- v[1]) / v[1]) *fac)
      out[,c("suitability_mean","suitability_q25", "suitability_q50", "suitability_q75")] <- apply(
        out[,c("suitability_mean","suitability_q25", "suitability_q50", "suitability_q75")], 2, relChange)
      out$suitability_avggain[is.na(out$suitability_avggain)] <- 0
      out$suitability_avggain <- out$suitability_avggain+1
      out$suitability_avgloss[is.na(out$suitability_avgloss)] <- 0
      out$suitability_avgloss <- out$suitability_avgloss+1
    }
  } else {
    stop("TBD")
  }

  # Return output
  return(
    out
    )
}

#' Crop and project a stars raster `HACK`
#'
#' @description
#' The reprojection of WGS84 currently fails due to some unforeseen bug.
#' This function is meant to reproject back the lasyer
#' @param obj A ['stars'] object to be clipped and cropped.
#' @param template A ['Raster'] or ['sf'] object to which the object should be projected.
#' @keywords internal, scenario
#' @noRd
hack_project_stars <- function(obj, template){
  assertthat::assert_that(
    inherits(obj, "stars"),
    is.Raster(template) || inherits(template, "sf")
  )
  # Get tempdir
  td <- raster::tmpDir()

  # Get resolution
  bg <- stars::st_as_stars(template)

  # Get full dis
  full_dis <- stars::st_dimensions(obj)
  assertthat::assert_that(length(full_dis)<=3,msg = "Stars object can only have x,y,z dimension.")

  # Output
  out <- c()
  for(v in names(obj)){
    sub <- obj[v]
    stars::write_stars(sub, file.path(td, "ReprojectedStars.tif"))

    suppressWarnings(
      gdalUtils::gdalwarp(srcfile = file.path(td, "ReprojectedStars.tif"),
                          dstfile = file.path(td, "ReprojectedStars_temp.tif"),
                          s_srs = "EPSG:4296",
                          tr = raster::res(template),
                          te = raster::bbox(template),
                          t_srs = sp::proj4string(template))
    )
    oo <- stars::read_stars(file.path(td, "ReprojectedStars_temp.tif"),proxy = F)
    names(oo) <- v # Rename

    # provide to output
    out <- c(out, oo)
    rm(oo)
    try({file.remove(file.path(td, "ReprojectedStars.tif"),
                file.path(td, "ReprojectedStars_temp.tif"))},silent = TRUE)
  }
  # Reformat again
  out <- stars::st_as_stars(out)
  assertthat::assert_that(
    length(stars::st_get_dimension_values(bg, "x")) == length(stars::st_get_dimension_values(out, "x"))
  )
  # Now reset the dimensions and add to output
  dims <- stars::st_dimensions(out)
  # Replace the band variable with the original one
  names(dims)[3] <- "time"
  dims$time <- full_dis$time
  # And the x-y dimensions by the template values
  bg_dim <- stars::st_dimensions(bg)
  dims$x <- bg_dim$x; dims$y <- bg_dim$y
  stars::st_dimensions(out) <- dims
  out <-  stars::st_set_dimensions(out, xy = c("x","y"))
  assertthat::assert_that(
    length(out) == length(obj),
    stars:::is_regular_grid(out)
  )
  return(out)
}

#' Quick handy function to calculate the centre of a range
#'
#' @param ras A [`RasterLayer`] object for which the centre of the range is to be calculated.
#' If the distribution is continuous, then the centre is calculated as the value centre to all non-NA values.
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
