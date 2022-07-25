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
approximate_gaps <- function(env, date_interpolation = "year"){
  assertthat::assert_that(
    inherits(env, "stars"),
    is.character(date_interpolation)
  )
  # Get individual time steps at interval
  # Check whether any are missing
  # Linearly approximate all attributes for new object

}
