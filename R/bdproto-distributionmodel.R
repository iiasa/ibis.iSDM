#' @include utils.R bdproto.R
NULL

#' @export
if (!methods::isClass("DistributionModel")) methods::setOldClass("DistributionModel")
NULL

#' Prototype for the trained Model object
#'
#' All trained Models should inherit the options here
#'
#' @name DistributionModel-class
#'
#' @aliases DistributionModel
NULL

#' @export
DistributionModel <- bdproto(
  "DistributionModel",
  name = character(0),
  id = new_waiver(),
  fits = list(), # List of fits with data
  # Print message with name of Model
  print = function(self) {
    message(self$name)
  },
  # Show the name of the Model
  show = function(self) {
    self$name
  },
  # Get specific fit from this Model
  get_data = function(self, x) {
    if (!x %in% names(self$fits))
      return(new_waiver())
    return(self$fits[[x]])
  },
  # Set fit for this Model
  set_data = function(self, x, value) {
    self$fits[[x]] <- value
    invisible()
  }
)
