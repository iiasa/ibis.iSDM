#' Approximate missing time steps between dates
#'
#' @description This function linearly approximates shares between time steps,
#' so that gaps for instance between 2010 and 2020 are filled with data for
#' 2010, 2011, 2012, etc.
#'
#' @param env A [`stars`] object.
#' @param date_interpolation [`character`] on how missing dates between events
#' should be interpolated. See [`project()`].
#'
#' @return [`logical`] indicating if the two [`SpatRaster-class`] objects have the same
#'
#' @keywords scenario
#'
#' @examples
#' \dontrun{
#'   # Interpolate stars stack
#'   sc <- interpolate_gaps( stack, "annual")
#' }
#'
#' @export
interpolate_gaps <- function(env, date_interpolation = "annual"){
  assertthat::assert_that(
    inherits(env, "stars"),
    is.character(date_interpolation)
  )
  check_package("dplyr")

  date_interpolation <- match.arg(date_interpolation,
                                  c("none", "yearly", "annual", "monthly", "daily"),
                                  several.ok = FALSE)
  if(date_interpolation=="none") return(env)

  # --- #
  # Get individual time steps at interval
  times <- stars::st_get_dimension_values(env, which = names(dim(env))[3], center = TRUE)
  times <- to_POSIXct(times)
  tzone <- attr(as.POSIXlt(times), "tzone") # Get timezone
  assertthat::assert_that(tzone != "", length(times)>=2)
  # Interpolate time steps
  inc <- switch (date_interpolation,
                 "yearly" = "year",
                 "annual" = "year",
                 "monthly" = "month",
                 "daily" = "day"
  )
  # Create new time layer
  new_times <- seq.Date(from = as.Date(times[1],tz = tzone),
                        to = as.Date(times[length(times)],tz = tzone),
                        by = inc)
  new_times <- to_POSIXct(new_times)
  ori.dims <- names(stars::st_dimensions(env))

  # Now for each variable, interpolate
  out <- list()
  if(getOption('ibis.setupmessages')) pb <- progress::progress_bar$new(total = length(env))
  for(v in names(env)){
    if(getOption('ibis.setupmessages')) pb$tick()

    # Get the variable
    o <- Reduce(c, stars_to_raster(env[v]) )

    # Create empty copies per times
    nt <- new_times[new_times %notin% terra::time(o)]
    oo <- rep(emptyraster(o, vals = NA), length(nt))
    terra::time(oo) <- nt
    oo <- c(o,oo)
    names(oo) <- paste0(v,'_',as.numeric(terra::time(oo)))

    # Linearly approximate gaps and save in list
    out[[v]] <- terra::approximate(x = oo,
                            method = "linear",
                            yleft = NA,
                            yright = NA) |> stars::st_as_stars()
  }

  # Combine
  out <- Reduce(c, out)
  # Rename again
  names(out) <- names(env)

  # Checks
  assertthat::assert_that(
    is.stars(out),
    all( names(out) %in% names(env) )
  )
  return(out)
}

#' Aggregate stars variables across dimensions
#'
#' @description Small helper function that acts a wrapper to combine 2 or more
#' variables in a `stars` object together.
#'
#' @param obj A [`stars`] object or a [`list`] that can be coerced to one.
#' @param vars A [`vector`] describing the variables to be combined. Has to be
#'   of length two or greater.
#' @param newname A [`character`] with the new name for the variable.
#' @param weights An optional variable layer to use for weighting.
#' @param fun A function how the respective layers should be combined.
#'
#' @note Currently only works via matrix manipulation
#'
#' @keywords scenario, internal
#'
#' @examples
#' \dontrun{
#'  st_reduce(obj, vars = c('forestShare', 'forestShare.2'),
#'   newname = "forest",fun = "sum")
#' }
#'
#' @noRd
st_reduce <- function(obj, vars, newname, weights = NULL, fun = 'sum'){
  assertthat::assert_that(
    is.list(obj) || inherits(obj, 'stars'),
    is.vector(vars) && is.character(vars),
    length(vars) >=2,
    is.null(weights) || is.character(weights),
    is.character(fun),
    is.character(newname)
  )
  fun <- match.arg(fun, c("sum", "multiply", "divide", "subtract", "mean", "weighted.mean"), several.ok = FALSE)
  check_package('stars')
  # Convert to stars if not already
  if(!inherits(obj, 'stars')) obj <- stars::st_as_stars(obj)
  # Check that variables are present
  assertthat::assert_that(
    all(vars %in% names(obj))
  )
  # If weighted mean is chosen, ensure that weights are set and of equal length to vars.
  if(fun == "weighted.mean") assertthat::assert_that(is.character(weights),
                                                     length(vars) == length(weights),
                                                     all(weights %in% names(obj)))

  # Future?
  if(foreach::getDoParRegistered()){
    ibis_future(cores = getOption("ibis.nthread"), strategy = getOption("ibis.futurestrategy"))
    fut <- TRUE
  } else { fut <- FALSE }
  # --- #
  # First get all target variables and non-target variables
  target <- obj |> dplyr::select( dplyr::all_of(vars))
  non_target <- obj |> dplyr::select(-dplyr::all_of(vars))
  if(fun == "weighted.mean"){
    target_weights <- obj |> dplyr::select( dplyr::all_of(weights))
    # Normalize the weights
    no <- stars::st_as_stars( Reduce("+", target_weights) )
    stars::st_dimensions(no) <- stars::st_dimensions(target_weights)
    target_weights <- target_weights / no
  }

  # Now apply the function on the target
  # FIXME: Not working as intended and horribly slow
  # target <- stars::st_apply(target, dims, fun,
  #                           PROGRESS = TRUE,
  #                           FUTURE = fut,
  #                           ...
  #                           )
  if(fun == "weighted.mean") target <- target * target_weights

  what <- switch (fun,
    "sum" = "+",
    "multiply" = "*",
    "divide" = "/",
    "subtract" = "-",
    "mean" = "+",
    "weighted.mean" = "+"
  )
  new <- stars::st_as_stars( Reduce(what, target) )
  if(what == "mean") new <- new / length(target)
  if(what == "weighted.mean") new <- new / ( stars::st_as_stars( Reduce("+", target_weights) ) )
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
#' @description This is a small helper function to convert a [`stars`] object to
#' a [`SpatRaster`] object. It is possible to select the time frame as well. If
#' multiple \code{"which"} entries are specified, then a [`list`] will be
#' returned.
#'
#' @param obj A [`stars`] object with a \code{"time"} dimension at least.
#' @param which The time entry to use for subsetting. Can be single [`numeric`]
#' or a [`vector`] of numeric time entries corresponding to the time dimension
#' (Default: \code{NULL}).
#' @param template An optional [`SpatRaster`] template to which the output should
#' be aligned too.
#'
#' @returns A [`list`] containing [`SpatRaster`] objects.
#'
#' @keywords scenario, internal
stars_to_raster <- function(obj, which = NULL, template = NULL){
  assertthat::assert_that(
    inherits(obj, 'stars'),
    is.null(which) || is.numeric(which),
    is.null(template) || is.Raster(template)
  )
  # Take name of third band, assuming it to be time
  time_band <- names(dim(obj))[3]

  # if timeband is na and dimensions equal to 2, return SpatRaster
  if(is.na(time_band) && length(stars::st_dimensions(obj))<=2){
    out <- list()
    out[["SingleTimeStep"]] <- terra::rast(obj)
    return(out)
  }

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
    o <- obj |> stars:::slice.stars({{time_band}}, tt) |>
      terra::rast() # Or alternatively rast
    names(o) <- names(obj)

    # Reset times to the correct ones
    terra::time(o) <- as.Date( rep(times[tt], terra::nlyr(o)) )

    # Now transform the out put if template is set
    if(!is.null(template)){
      if(is.Raster(template)){
        # Check again if necessary to rotate
        if(sf::st_crs(o) != sf::st_crs(template)){
          o <- terra::project(x = o, y = template, method = "near")
          names(o) <- names(obj)
        }
        # Now crop and resample to target extent if necessary
        if(!terra::compareGeom(o, template, stopOnError = FALSE)){
          o <- terra::crop(o, template)
          o2 <- try({alignRasters(data = o,
                                  template = template,
                                  method = "near",
                                  func = "mean",
                                  cl = FALSE)
          },silent = TRUE)
          if(inherits(o2,"try-error")){
            o <- terra::resample(o, template,
                                 method = "near",
                                 threads = getOption("ibis.nthread"))
          } else { o <- o2; rm(o2)}
        }
      }
    } # End of template adjustments
    out[[paste0("time", times[tt])]] <- o
  }
  return( out )
}

#' Converts a raster object to stars
#'
#' @description This is a small helper function to convert a to a [`SpatRaster`]
#' object.
#'
#' @param obj A [`SpatRaster`] object with a \code{"time"} dimension at least
#' (checked via [`time`]).
#'
#' @returns A [`stars`] object with the formatted data.
#'
#' @seealso `stars_to_raster`
#' @keywords scenario, internal
#'
#' @examples
#' \dontrun{
#'  # Convert stars to SpatRaster
#'  stars_to_raster(obj)
#' }
#'
raster_to_stars <- function(obj){
  assertthat::assert_that(
    is.Raster(obj)
  )
  # Check that time dimension exist
  assertthat::assert_that( !is.null( terra::time(obj) ),
                           msg = "The supplied object requires a z dimension! Preferably provide a stars object.")
  assertthat::assert_that(!is.na(terra::crs(obj)),
                           msg = "Uniform projection for input raster is missing!")

  # Get time dimension
  times <- terra::time(obj)
  if(all(is.na( as.character(times) ))) stop("Predictor covariates are missing a time dimension! See terra::time() ")
  if(!all(inherits(times, "Date"))) times <- as.Date(times)
  prj <- sf::st_crs( terra::crs(obj) )

  # Convert to RasterStack and reset time dimension
  terra::time(obj) <- times
  # stars::make_intervals(times[1], times[2]) # For making intervals from start to end

  # Convert to stars step by step
  new_env <- list()
  for(i in 1:terra::nlyr(obj)){
    oo <- subset(obj, i)
    # Check if times are unique
    if(length(unique(times))==1) terra::time(oo) <- NULL
    suppressWarnings(  o <- stars::st_as_stars(oo) )
    # If CRS is NA
    if(is.na(sf::st_crs(o))) sf::st_crs(o) <- prj

    # Some hacky stuff since stars is not behaving as intended
    dims <- stars::st_dimensions(o)
    dims$time <- stars:::create_dimension(values = times[i])
    o <- stars::st_redimension(o,new_dims = dims)

    new_env[[names(obj)[i]]] <- o
  }

  new_env <- do.call(c, new_env)
  assertthat::assert_that(inherits(new_env, "stars"),
                          stars::st_dimensions(new_env) |> length() == 3)

  return(new_env)
}

#' This function add layers from a SpatRaster to a stars object
#'
#' @description Often it is necessary to add static variables to existing stars
#' objects. These will be replicated across the time dimension. This function is
#' a small helper function that allows the addition of said raster stacks to a
#' stars object.
#'
#' @param obj A [`stars`] object with a time dimension (\code{"time"}).
#' @param new A [`SpatRaster`] object with additional covariates to be added.
#'
#' @returns A [`stars`] object with the names of the [`SpatRaster`] object added.
#'
#' @keywords scenario, internal
st_add_raster <- function(obj, new){
  assertthat::assert_that(
    inherits(obj, "stars"),
    is.Raster(new),
    terra::nlyr(new) >= 1
  )

  # Check whether there are any variables in the stars object already, if so drop
  if(any(names(new) %in% names(obj))){
    myLog("[Starting]", "yellow", "Duplicate variables in stars and new objects.")
    new <- new[[ which( names(new) %in% names(obj)) ]]
  }

  full_dims <- stars::st_dimensions(obj)
  # Get times objects
  time_name <- names(full_dims)[3]
  times <- rep(stars::st_get_dimension_values(obj, time_name))

  # Now loop through each layer and add it to the target file
  for(lyr in names(new)){
    s <- terra::rast(replicate(length(times), new[[lyr]])) |>
      stars::st_as_stars()
    names(s) <- lyr

    stars::st_dimensions(s) <- full_dims
    obj <- c(obj, s)
  }
  assertthat::assert_that(
    all(all(names(new) %in% names(obj)))
  )
  return(obj)
}

#' Summarize results from scenario projection object
#'
#' @description This is a wrapper function to summarize the output of a scenario
#' projection. The output will contain the average change in the layer per time
#' step. A parameter called \code{"relative"} can be set to calculate relative
#' change instead.
#'
#' @param scenario A [`stars`] object with a time dimension.
#' @param relative A [`logical`] check whether to calculate relative changes instead.
#'
#' @keywords internal, scenario
#'
#' @noRd
summarise_projection <- function(scenario, fun = "mean", relative = TRUE){
  assertthat::assert_that(
    is.list(scenario) || inherits(scenario, "stars"),
    length(dim(scenario))==3,
    is.logical(relative)
  )
  fun <- match.arg(fun, c("mean", "sum"),several.ok = FALSE)

  # Convert to scenarios to data.frame
  df <- stars:::as.data.frame.stars(stars::st_as_stars(scenario)) |> (\(.) subset(., stats::complete.cases(.)))()
  names(df) <- c("x", "y", "band", "suitability")
  # Add grid cell grouping
  df <- df |> dplyr::group_by(x,y) |> dplyr::mutate(id = dplyr::cur_group_id()) |>
    dplyr::ungroup() |> dplyr::select(-x,-y) |>
    dplyr::arrange(id, band)

  # Summarize the overall moments
  if(fun == "mean"){
    # Check if has unit, if so deparse
    if(inherits(df$suitability, 'units')) df$suitability <- as.numeric(df$suitability)
    out <- df  |>
      dplyr::filter(suitability > 0) |>
      dplyr::group_by(band) |>
      dplyr::summarise(suitability_mean = mean(suitability, na.rm = TRUE),
                       suitability_q25 = quantile(suitability, .25),
                       suitability_q50 = quantile(suitability, .5),
                       suitability_q75 = quantile(suitability, .75))
    # Total amount of area lost / gained / stable since previous time step
    totchange_occ <- df |>
      dplyr::group_by(id)  |>
      dplyr::mutate(change = (suitability - dplyr::lag(suitability)) )  |>  dplyr::ungroup()
    o <- totchange_occ |> dplyr::group_by(band) |>
      dplyr::summarise(suitability_avggain = mean(change[change > 0]),
                       suitability_avgloss = mean(change[change < 0]))

    out <- out |> dplyr::left_join(o, by = "band")
    if(relative){
      # Finally calculate relative change to baseline (first entry) for all entries where this is possible
      relChange <- function(v, fac = 100) (((v- v[1]) / v[1]) *fac)
      out[,c("suitability_mean","suitability_q25", "suitability_q50", "suitability_q75")] <- apply(
        out[,c("suitability_mean","suitability_q25", "suitability_q50", "suitability_q75")], 2, relChange)
    }
  } else if(fun == "sum") {
    # Check if has unit, if so deparse
    if(inherits(df$suitability, 'units')) df$suitability <- as.numeric(df$suitability)
    out <- df |>
      dplyr::filter(suitability > 0) |>
      dplyr::group_by(band) |>
      dplyr::summarise(suitability_sum = sum(suitability, na.rm = TRUE),
                       suitability_q25 = quantile(suitability, .25),
                       suitability_q50 = quantile(suitability, .5),
                       suitability_q75 = quantile(suitability, .75))
    # Total amount of area lost / gained / stable since previous time step
    totchange_occ <- df |>
      dplyr::group_by(id) |>
      dplyr::mutate(change = (suitability - dplyr::lag(suitability)) ) |> dplyr::ungroup()
    o <- totchange_occ |> dplyr::group_by(band) |>
      dplyr::summarise(suitability_avggain = sum(change[change > 0]),
                       suitability_avgloss = sum(change[change < 0]))

    out <- out |> dplyr::left_join(o, by = "band")
    if(relative){
      # Finally calculate relative change to baseline (first entry) for all entries where this is possible
      relChange <- function(v, fac = 100) (((v- v[1]) / v[1]) *fac)
      out[,c("suitability_sum","suitability_q25", "suitability_q50", "suitability_q75")] <- apply(
        out[,c("suitability_sum","suitability_q25", "suitability_q50", "suitability_q75")], 2, relChange)
    }
  }

  # Return output
  return(out)
}

#' Duplicate a provided stars raster
#'
#' @description This function duplicates a provded [`stars`] object along a
#' given dimension (usually time). This can be useful for simply multiplication
#' tasks, e.g. multiplying every attribute by another stars object that needs to
#' have the same dimensions.
#'
#' @param obj A [`stars`] object that is to be duplicated.
#' @param dim A dimensions file return from \code{st_dimensions(...)}.
#' @param dimname A [`character`] of the dimension name to be used.
#'
#' @returns A [`stars`] object.
#'
#' @examples
#' \dontrun{
#' o <- st_rep(obj, dim)
#' }
#'
#' @keywords internal, scenario
#'
#' @noRd
st_rep <- function(obj, dim, dimname = "time"){
  assertthat::assert_that(
    inherits(obj, "stars"),
    inherits(dim, "dimensions"),
    dimname %in% names(dim)
  )
  # Check that coordinate systems and values are identical
  assertthat::assert_that(
    dim[[1]]$to == stars::st_dimensions(obj)[[1]]$to,
    dim[[1]]$refsys == sf::st_crs(obj),
    msg = "Coordinate system between provided dimension and obj are not identical!"
  )

  # Make a sequence of the given dimension
  it <- seq(dim[[dimname]]$from, dim[[dimname]]$to-1)

  # Make a dummy
  new <- obj
  for(i in it) new <- c(new, obj,along = 3)

  stars::st_dimensions(new) <- dim # Redimension

  assertthat::assert_that(inherits(new, "stars"))
  return(new)
}

#' Summarize change before to after
#'
#' @description This is a wrapper function to summarize the output of a scenario
#' projection, but specifically calculates statistics of change for two time
#' steps, a before and after step.
#'
#' @param scenario A [`stars`] object with a time dimension.
#'
#' @note This function currently requires the \code{"geosphere"} package
#' installed.
#'
#' @references
#' * Godsoe, W. (2014). Inferring the similarity of species distributions using
#' Species’ Distribution Models. Ecography, 37(2), 130-136.
#'
#' @keywords internal, scenario
#'
#' @noRd
summarise_change <- function(scenario){
  assertthat::assert_that(
    inherits(scenario, "stars")
  )
  # Check that geosphere is installed and loaded
  check_package("geosphere")
  if(!("geosphere" %in% loadedNamespaces()) || ('geosphere' %notin% utils::sessionInfo()$otherPkgs) ) {
    try({requireNamespace('geosphere');attachNamespace("geosphere")},silent = TRUE)
  } else {
    if(getOption("ibis.setupmessages")) myLog("[Summary]","red","This summary function requires the geosphere package.")
  }

  # Get the current and future
  ss <- stars_to_raster(scenario)
  # Time period
  times <- stars::st_get_dimension_values(scenario, 3,center = TRUE)
  current <- ss[[1]]
  future <- ss[[length(ss)]]
  times_length <- round(as.numeric(difftime(times[length(times)], times[1], units = "weeks"))/52.25,0)
  rm(ss)

  # Calculate the area and  units
  ar <- st_area(scenario)
  ar_unit <- units::deparse_unit(ar$area)
  if(ar_unit == "m2"){
    ar_unit <- "ha"
    mult <- 0.0001
  } else { mult <- 1}
  ar <- terra::rast(ar) # or sds if that one fails

  # --- #
  val <- c("Current range", "Future range", "Unsuitable",
           "Loss", "Gain", "Stable", "Percent loss",
           "Percent gain", "Range change", "Percent change",
           "Sorensen index", "Centroid distance", "Centroid change direction")
  change <- data.frame(category = val,
                       period = c(times[1] |> as.character(),
                                  times[length(times)] |> as.character(), rep(paste0(times_length, " years"), 11 ) ),
                       value = NA,
                       unit = c(rep(ar_unit,6), "%", "%", ar_unit, "%", "similarity", NA, "deg"))
  change$value[1] <- terra::global((current) * terra::cellSize(current, unit = "km"), "sum", na.rm = TRUE)[,1] * mult
  change$value[2] <- terra::global((future) * terra::cellSize(future, unit = "km"), "sum", na.rm = TRUE)[,1] * mult
  assertthat::assert_that(is.numeric(change$value))

  # Check that is binary thresholded
  rr <- terra::lapp(c(current, future), fun = function(x, y){x + y * 2})
  change$value[3] <- terra::global((rr == 0) * terra::cellSize(current, unit = "km"), "sum", na.rm = TRUE)[,1] * mult
  change$value[4] <- terra::global((rr == 1) * terra::cellSize(current, unit = "km"), "sum", na.rm = TRUE)[,1] * mult
  change$value[5] <- terra::global((rr == 2) * terra::cellSize(current, unit = "km"), "sum", na.rm = TRUE)[,1] * mult
  change$value[6] <- terra::global((rr == 3) * terra::cellSize(current, unit = "km"), "sum", na.rm = TRUE)[,1] * mult
  change$value[7] <- change$value[4] / change$value[1] * 100
  change$value[8] <- change$value[5] / change$value[1] * 100
  change$value[9] <- change$value[2] - change$value[1]
  change$value[10] <- change$value[9] / sum(c(change$value[3], change$value[4])) * 100

  # Sorensen similarity index
  change$value[11] <- 2 * terra::global(rr == 3, "sum", na.rm = TRUE)[,1] / (terra::global(current, "sum", na.rm = TRUE)[,1] + terra::global(future, "sum", na.rm = TRUE)[,1])

  # Calculate distance between centroids
  sf1 <- calculate_range_centre(current, spatial = TRUE)
  sf2 <- calculate_range_centre(future, spatial = TRUE)
  dis <- sf::st_distance(sf1, sf2, by_element = FALSE)
  dis_unit <- units::deparse_unit(dis)
  # Convert units if meter
  if( dis_unit == "m") {mult <- 0.001; dis_unit = "km" } else { mult <- 1}
  change$value[12] <- as.vector(dis) * mult
  change$unit[12] <- dis_unit

  # Calculate direction between centroids
  change$value[13] <- geosphere::finalBearing(as_Spatial(sf1 |> sf::st_transform(crs = sf::st_crs(4326))),
                                              as_Spatial(sf2 |> sf::st_transform(crs = sf::st_crs(4326))))

  change <- change |> tibble::as_tibble()
  return(change)
}

#' Label patches and apply a minimum size constraint
#'
#' @description This is a
#'
#' @param obj A ['stars'] object to be clipped and cropped.
#'
#' @keywords internal, scenario
#'
#' @noRd
st_minsize <- function(obj, value, unit = "km2",  establishment_step = FALSE){
  assertthat::assert_that(
    inherits(obj, "stars") || is.Raster(obj),
    is.numeric(value),
    is.character(unit),
    is.logical(establishment_step)
  )

  # Match unit
  unit <- match.arg(arg = unit,
                    choices = c("km2", "ha", "pixel"), several.ok = FALSE)

  # If ha convert to km2
  if(unit == "ha"){
    value <- value / 100
    unit <- "km2"
  }

  # Original
  ori.obj <- obj

  # Procedure for raster object
  if(is.Raster(obj)){
    assertthat::assert_that(utils::hasName(obj, "threshold"))
    obj <- obj["threshold"]

    # Check unique values
    assertthat::assert_that(length(unique(obj)[,1]) <=2)

    if(unit == "pixel"){
      # Sieve by pixel size
      assertthat::assert_that(value>1,
                              msg = "Pixel units need to be larger than 1")
      obj[obj == 0] <- NA
      new <- terra::sieve(obj,
                          value,
                          directions = 8)
      new <- terra::mask(new, ori.obj)
    } else {
      # Now first label
      labs <- terra::patches(obj)

      # Then calculate area in km2
      ar <- terra::mask(terra::cellSize(labs, unit = "km"), labs)

      # Now summarize per label the area
      zt <- terra::zonal(ar, labs)

      # Identify all grid cells smaller than value and
      not_big_enough <- zt$patches[which(zt$area < value)]
      if(length(not_big_enough)>0){
        new <- labs
        new[new %in% not_big_enough] <- 0
        new[new>0] <- 1 # Set remaining cells back to threshold
        names(new) <- names(obj)
      } else {
        new <- ori.obj
      }
    }

    if(terra::global(new, "max", na.rm = TRUE)[,1] == 0){
      if(getOption('ibis.setupmessages')) myLog('[Done]','yellow',
                                                paste0('Min size constraint removed all areas!'))
    }
  } else {
    # Check that threshold is present
    assertthat::assert_that(utils::hasName(obj, "threshold"))
    obj <- obj["threshold"]

    # Also check that time dimension is present
    assertthat::assert_that(
      names(stars::st_dimensions(obj)[3]) %in% c("time", "band")
    )

    # Save dimensions
    ori.dims <- stars::st_dimensions(obj)

    # Now for each time element, convert to raster and label
    times <- stars::st_get_dimension_values(obj, 3)
    new <- terra::rast()
    for(step in 1:length(times)){
      ras <- stars_to_raster(obj, which = step)[[1]]

      # For stars objects, convert to raster and self-call
      suppressWarnings(
        new <- c(new, st_minsize(obj = ras,
                 value = value,
                 unit = unit)
        )
      )
    }
    assertthat::assert_that(is.Raster(new))
    # Convert back to stars
    new <- stars::st_as_stars(new,crs = sf::st_crs(obj))
    names(new) <- 'threshold'
    # Correct band if different
    if(all(!stars::st_get_dimension_values(obj, 3) != stars::st_get_dimension_values(new, 3 ))){
      new <- stars::st_set_dimensions(new, 3, values = stars::st_get_dimension_values(obj, 3))
    }
  }
  # Return the result
  return(new)
}

#' Crop and project a stars raster `HACK`
#'
#' @description The reprojection of WGS84 currently fails due to some unforeseen
#' bug. This function is meant to reproject back the layer.
#'
#' @param obj A ['stars'] object to be clipped and cropped.
#' @param template A ['SpatRaster'] or ['sf'] object to which the object should
#' be projected.
#' @param use_gdalutils (Deprecated) [`logical`] on to use gdalutils hack around.
#'
#' @keywords internal, scenario
#'
#' @noRd
hack_project_stars <- function(obj, template, use_gdalutils = TRUE){
  assertthat::assert_that(
    inherits(obj, "stars"),
    is.Raster(template) || inherits(template, "sf"),
    is.logical(use_gdalutils)
  )
  # Get tempdir
  td <- terra::terraOptions(print = FALSE)[['tempdir']]

  # Get resolution
  bg <- stars::st_as_stars(template)

  # Get full dis
  full_dis <- stars::st_dimensions(obj)
  assertthat::assert_that(length(full_dis)<=3, msg = "Stars object can only have x,y,z dimension.")

  # Output
  out <- c()
  for(v in names(obj)){
    sub <- obj[v]

    if(use_gdalutils){
      check_package("gdalUtilities")
      # Write output
      stars::write_stars(sub, file.path(td, "ReprojectedStars.tif"))
      suppressWarnings(
        gdalUtilities::gdalwarp(srcfile = file.path(td, "ReprojectedStars.tif"),
                            dstfile = file.path(td, "ReprojectedStars_temp.tif"),
                            s_srs = "EPSG:4296",
                            tr = terra::res(template),
                            te = terra::ext(template) |> sf::st_bbox(),
                            t_srs = sf::st_crs(template)$proj4string)
      )
    } else {
      # Try and use terra
      stars::write_stars(sub, file.path(td, "ReprojectedStars.tif"))
      # Re project with terra
      temp <- terra::rast(x = file.path(td, "ReprojectedStars.tif"))
      temp <- terra::project(x = temp,
                             y = template,
                             align = TRUE,
                             gdal = TRUE,
                             threads = FALSE)
      terra::writeRaster(temp, file.path(td, "ReprojectedStars_temp.tif"),overwrite = TRUE)
    }
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

#' Quick handy function to calculate an area-weighted centre of a range
#'
#' @param layer A [`SpatRaster`] or [`sf`] object for which the centre of the
#'   range is to be calculated. If the distribution is continuous, then the
#'   centre is calculated as the value centre to all non-NA values.
#' @param spatial A [`logical`] of whether outputs should be returned as
#'   spatial.
#' @aliases calculate_range_centre
#' @keywords scenario, internal
#' @noRd
calculate_range_centre <- function(layer, spatial = TRUE) {
  assertthat::assert_that(
    is.Raster(layer) || inherits(layer, "sf")
  )

  # If layer is a raster
  if(is.Raster(layer)){
    assertthat::assert_that(
      length( unique(layer)[,1] ) == 2,
      terra::global(layer, 'max', na.rm = TRUE) == 1
    )
    # Calculate area-weighted centre
    r_wt <- terra::cellSize(layer, unit = "km")
    values(r_wt)[is.na(values(layer))] <- NA

    # Make a spatial point layer
    spdf <- terra::as.points( c(layer, r_wt)) |> sf::st_as_sf()
    spdf <- spdf[which(spdf[[1]]>0), ] # Get only non-zero values

    if(is.na(sf::st_crs(spdf))) stop("Unprojected layer found. Check projections throughout!")
    # If long-latitude, convert to google mercator for calculating the centroids
    if(sf::st_is_longlat(spdf) ){
      ori.proj <- sf::st_crs(spdf)
      spdf <- sf::st_transform( spdf, crs = sf::st_crs(3857))
    } else { ori.proj <- sf::st_crs(spdf) }

    p <- sf::st_drop_geometry(spdf[, names(spdf)[2] ])[,1]
    # Calculate weighted centroid
    Xw <- sum(sf::st_coordinates(spdf)[,1] * p)
    Yw <- sum(sf::st_coordinates(spdf)[,2] * p)
    wX <- Xw/sum(p)
    wY <- Yw/sum(p)
    xy <- data.frame(ID = 1, name = names(layer), X=wX, Y=wY)
    cent <- sf::st_as_sf(xy, coords = c("X", "Y"),
                       crs = sf::st_crs(spdf), agr = "constant")
    # Convert back to original projection
    cent <- sf::st_transform(cent, ori.proj)

  } else {
    if(is.na(sf::st_crs(layer))) stop("Unprojected layer found. Check projections throughout!")
    # If long-latitude, convert to google mercator for calculating the centroids
    if(sf::st_is_longlat(layer) ){
      ori.proj <- sf::st_crs(layer)
      layer <- sf::st_transform( layer, crs = sf::st_crs(3857))
    } else { ori.proj <- sf::st_crs(layer) }

    if(unique(sf::st_geometry_type(layer)) %in% c("POLYGON", "MULTIPOLYGON")){
      # Cast them into a multi-polygon
      cent <- sf::st_combine(layer) |> sf::st_centroid() |> sf::st_as_sf()
    } else if(unique(sf::st_geometry_type(layer)) %in% c("POINT", "MULTIPOINT")){
      cent <- sf::st_combine(layer) |> sf::st_centroid() |> sf::st_as_sf()
    } else {
      stop("Centroid calculations not implemented!")
    }
    # Convert back to original projection
    cent <- sf::st_transform(cent, ori.proj)
    cent$ID = 1
  }

  if(!spatial){
    cent$X <- sf::st_coordinates(cent)[,1]
    cent$Y <- sf::st_coordinates(cent)[,2]
    cent <- sf::st_drop_geometry(cent)
  }
  return(cent)
}
