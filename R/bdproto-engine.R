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
  # TODO: Set up parameter objects or cut here
#  parameters = parameters(),
  data = list(),
  calculate = function(...) stop("Engine is missing a fitting method"),
  # Print message with name of engine
  print = function(self) {
    message(self$name)
  },
  # Show the name of the Engine
  show = function(self) {
    self$name
  },
  # Get run parameters
  get_parameter = function(self, x) {
    self$parameters$get(x)
  },
  # Set run parameters
  set_parameter = function(self, x, value) {
    self$parameters$set(x, value)
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
