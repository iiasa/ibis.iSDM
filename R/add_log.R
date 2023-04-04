#' @include utils.R bdproto.R bdproto-biodiversitydistribution.R bdproto-log.R
NULL

#' Adds a log file to distribution object
#'
#' @description
#' This function allows to specify a file as [Log-class] file, which is used to save all
#' console outputs, prints and messages.
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param filename A [`character`] object.
#' The destination must be writeable and filename ends with \code{'txt'}.
#' @aliases add_log
#' @returns Adds a log file to a [`distribution`] object.
#' @examples
#' \dontrun{
#'  x <- distribution(background) |>
#'     add_log()
#'  x
#' }
#' @name add_log
NULL

#' @name add_log
#' @rdname add_log
#' @exportMethod add_log
#' @export
methods::setGeneric(
  "add_log",
  signature = methods::signature("x", "filename"),
  function(x, filename) standardGeneric("add_log"))

#' @name add_log
#' @rdname add_log
#' @usage \S4method{add_log}{BiodiversityDistribution, character}(x, filename)
methods::setMethod(
  "add_log",
  methods::signature(x = "BiodiversityDistribution", filename = "character"),
  function(x, filename ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.character(filename),
                            assertthat::has_extension(filename,'txt') )

    # Messenger
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding log file...')

    # Check whether a log is already present in the distribution file
    if(!is.Waiver(x$log)) myLog('[Setup]','yellow','Overwriting previous set log file.')

    # Finally set the data to the BiodiversityDistribution object
    x$set_log(
      bdproto(NULL, Log,
              filename = filename,
              output = new_waiver()
      )
    )
  }
)
