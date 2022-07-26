#' Approximate missing time steps between dates
#'
#' @description
#' This function linearly approximates shares between time steps, so that gaps for instance
#' between 2010 and 2020 are filled with data for 2010, 2011, 2012, etc.
#' @param env A [`Raster-class`] object.
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

  stop("TBD. Not yet implemented")

  # --- #
  # Get individual time steps at interval
  times <- stars::st_get_dimension_values(env, which = names(dim(env))[3], center = TRUE)
  times <- to_POSIXct(times)
  tzone <- attr(as.POSIXlt(times), "tzone")[2] # Get timezone
  assertthat::assert_that(tzone != "")
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
  # new <- stars::st_redimension(out, along = list(time = new_times))
}
