#' @include utils.R bdproto.R
NULL

#' @export
if (!methods::isClass("Settings")) methods::setOldClass("Settings")
NULL

#' Prototype for model settings object
#'
#' @description
#' Basic proto object for Settings object, a List that stores settings used
#' related to model training
#'
#' @family bdproto
#' @keywords bdproto
#'
#' @name Settings-class
NULL

#' @rdname Settings-class
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
  # Summary call
  summary = function(self){
    if(self$length()==0) return(NULL)
    # Loop through each entry and try to summarize
    # Do so by estimating name, type or class and value (if possible)
    d <- self$data
    o <- data.frame(name = character(), type = character(), value = character() )
    for(entry in names(d)){
      if(inherits(d[[entry]], "sf")){
        val <- sf::st_geometry_type( d[[entry]] )[1]
      } else {
        val <- try({ as.character( d[[entry]]) },silent = TRUE)
      }
      if(inherits(val, 'try-error')) val <- NA
      e <- data.frame(name = entry, type = class(d[[entry]])[1], value = val)
      o <- rbind(o, e)
    }
    # Return the data.frame
    return(
      o
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
