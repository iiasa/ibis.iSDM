#' @include utils.R bdproto.R bdproto-biodiversitydistribution.R bdproto-log.R
NULL

#' Adds a log file to distribution object to write messages and output in there
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param filename A [`character`] object. Destination must be writeable and filename ends with 'txt'
#' @details TBD
#' @section Notes:
#' @aliases add_log
#' @references
#'
#' @examples
#' \dontrun{
#'  TBD
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

    # Check whether a log is already present in the distribution file
    if(!is.Waiver(x$log)) message('Overwriting previous set log file.')

    # Finally set the data to the BiodiversityDistribution object
    x$set_log(
      bdproto(NULL, Log,
              filename = filename,
              output = new_waiver()
      )
    )
  }
)
