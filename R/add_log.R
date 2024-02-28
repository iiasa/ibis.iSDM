#' @include class-biodiversitydistribution.R class-log.R
NULL

#' Adds a log file to distribution object
#'
#' @description This function allows to specify a file as [Log-class] file,
#' which is used to save all console outputs, prints and messages.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param filename A [`character`] object. The destination must be writeable and
#' filename ends with \code{'txt'}.
#'
#' @returns Adds a log file to a [`distribution`] object.
#'
#' @examples
#' \dontrun{
#'  x <- distribution(background) |>
#'     add_log()
#'  x
#' }
#'
#' @name add_log
NULL

#' @rdname add_log
#' @export
methods::setGeneric(
  "add_log",
  signature = methods::signature("x", "filename"),
  function(x, filename) standardGeneric("add_log"))

#' @rdname add_log
methods::setMethod(
  "add_log",
  methods::signature(x = "BiodiversityDistribution", filename = "character"),
  function(x, filename ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.character(filename),
                            assertthat::has_extension(filename,'txt') )

    # Messenger
    if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Setup]','green','Adding log file...')

    # Check whether a log is already present in the distribution file
    if(!is.Waiver(x$log)) myLog('[Setup]','yellow','Overwriting previous set log file.')

    # Finally set the data to the BiodiversityDistribution object
    l <- Log$new(filename = filename,
                 output = new_waiver())

    # Make a clone copy of the object
    y <- x$clone(deep = TRUE)

    x$set_log(l)
  }
)
