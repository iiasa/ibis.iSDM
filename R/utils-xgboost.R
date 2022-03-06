#' Split up a factor variable into binary components
#'
#' @param df A [`vector`] object containing the factor variables
#' @param name Name for the new object
#' @keywords utils, internal
#' @noRd
explode_factor <- function(df, name = "facvar"){
  assertthat::assert_that(
    is.data.frame(df) || is.factor(df),
    all(is.factor(df)),
    is.character(name)
  )
  z <- as.data.frame(
    outer(df, levels(df), function(w, f) ifelse(w == f, 1, 0))
    )
  names(z) <- paste(name, levels(df), sep = ".")
  return(z)
}
