#' @include utils.R bdproto.R
NULL

#' @export
if (!methods::isClass("Engine")) methods::setOldClass("Engine")
NULL

#' Engine prototype.
#'
#' Basic object for engine, all other engines inherit from here
#'
#' @name Engine-class
#'
#' @aliases Engine
NULL

#' @export
Engine <- bdproto(
  "Engine",
  name = character(0),
  data = list(),
  # Functions to be overwritten when mesh is specified
  setup = function(...) stop("Engine is missing functions to setup settings."),
  train = function(...) stop("Engine is missing a fitting method."),
  # Print message with name of engine
  print = function(self) {
    message( paste0('Engine: ', self$name) )
  },
  # Show the name of the Engine
  show = function(self) {
    self$name
  },
  # Get specific data from this engine
  get_data = function(self, x) {
    if (!x %in% names(self$data))
      return(new_waiver())
    return(self$data[[x]])
  },
  # Set data for this engine
  set_data = function(self, x, value) {
    self$data[[x]] <- value
    invisible()
  }
)
