#' @include utils.R bdproto.R
NULL

#' @export
if (!methods::isClass("Settings")) methods::setOldClass("Settings")
NULL

#' Prototype for model settings object
#'
#' Basic proto object for Settings object, a List that stores settings used related to model training
#'
#' @name Settings-class
#'
#' @aliases Settings
NULL

#' @export
Settings <- bdproto(
  "Settings",
  modelid = new_waiver(), # The id to which the model belong
  data = list(),
  # Print
  print = function(self) {
    message(paste0( self$show(),': ',as.character( self$modelid ) ))
  },
  # Show
  show = function(self){
    'Model settings'
  },
  # Number of options
  length = function(self){
    base::length( self$data )
  },
  # Computation duration convenience function
  duration = function(self){
    assertthat::assert_that('start.time' %in% names(self$data),
                            'end.time' %in% names(self$data))
    return(
      self$data$end.time - self$data$start.time
      )
  },
  # Get setting
  get = function(self, what){
    assertthat::assert_that(is.character(what),
                            what %in% names(self$data))
    return(self$data[[what]])
  },
  # Set setting
  set = function(self, what, x){
    assertthat::assert_that(is.character(what))
    self$data[[what]] <- x
    invisible()
  }
)
