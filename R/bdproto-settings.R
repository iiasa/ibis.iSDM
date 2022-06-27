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
#' @aliases Settings
#' @family bdproto
#' @keywords bdproto
NULL

#' @export
Settings <- bdproto(
  "Settings",
  name = character(0),
  modelid = new_waiver(), # The id to which the model belong
  data = list(),
  # Print
  print = function(self) {
    # Get number of settings set
    ss <- sum(sapply(self$data, function(x) !is.Waiver(x)))
    ss <- text_green(ss)
    message(paste0( self$show(),': ', ss, ' parameters' ))
  },
  # Show
  show = function(self){
    paste(self$name, 'Settings')
  },
  # Number of options
  length = function(self){
    base::length( self$data )
  },
  # Computation duration convenience function
  duration = function(self){
    assertthat::assert_that('start.time' %in% names(self$data),
                            'end.time' %in% names(self$data),
                            msg = "Model duration not computed or model not fitted!")
    return(
      self$data$end.time - self$data$start.time
      )
  },
  # Get setting
  get = function(self, what){
    assertthat::assert_that(is.character(what))
    if(what %in% names(self$data)){
      return(self$data[[what]])
    } else { return(new_waiver()) }
  },
  # Set setting
  set = function(self, what, x, copy = FALSE){
    assertthat::assert_that(is.character(what))
    # Get biodiversity dataset collection
    ff <- self$data
    # Set the object
    ff[[what]] <- x
    if(copy){
      # If the object is an actual copy of the original one
      bdproto(NULL, self, data = ff )
    } else {
      self$data <- ff
    }
  }
)
