if (!methods::isClass("Settings")) methods::setOldClass("Settings")

#' Prototype for model settings object
#'
#' @description
#' Basic [`R6`] object for Settings object, a List that stores settings used
#' related to model training.
#'
#' @keywords classes
#'
#' @name Settings-class
NULL

#' @rdname Settings-class
#' @export
Settings <- R6::R6Class(
  "Settings",
  public = list(
    #' @field name The default name of this settings as [`character`].
    #' @field modelid A [`character`] of the model id this belongs to.
    #' @field data A [`list`] of contained settings.
    name = character(0),
    modelid = new_waiver(), # The id to which the model belong
    data = list(),

    #' @description
    #' Initializes the object and creates an empty list
    #' @return NULL
    initialize = function(){
    },

    #' @description
    #' Print the names and properties of all Biodiversity datasets contained within
    #' @return A message on screen
    print = function(){
      # Get number of settings set
      ss <- sum(sapply(self$data, function(x) !is.Waiver(x)))
      ss <- text_green(ss)
      message(paste0( self$show(),': ', ss, ' parameters' ))
    },

    #' @description
    #' Shows the name and the settings
    #' @return A [`character`] of the name and settings.
    show = function(){
      paste(self$name, 'Settings')
    },

    #' @description
    #' Number of options
    #' @return A [`numeric`] with the number of options.
    length = function(){
      base::length( self$data )
    },

    #' @description
    #' Computation duration convenience function
    #' @return The amount of time passed for model fitting if found.
    duration = function(){
      assertthat::assert_that('start.time' %in% names(self$data),
                              'end.time' %in% names(self$data),
                              msg = "Model duration not computed or model not fitted!")
      return(
        self$data$end.time - self$data$start.time
      )
    },

    #' @description
    #' Summary call of the contained parameters
    #' @return A [`list`] with the parameters in this object.
    summary = function(){
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

    #' @description
    #' Get a specific setting
    #' @param what A [`character`] with the respective setting.
    #' @return The setting if found in the object.
    get = function(what){
      assertthat::assert_that(is.character(what))
      if(what %in% names(self$data)){
        return(self$data[[what]])
      } else { return(new_waiver()) }
    },

    #' @description
    #' Set new settings
    #' @param what A [`character`] with the name for the new settings.
    #' @param x The new setting to be stored. Can be any object.
    #' @param copy [`logical`] on whether a new settings object is to be created.
    #' @return The setting if found in the object.
    set = function(what, x, copy = FALSE){
      assertthat::assert_that(is.character(what))
      # Get biodiversity dataset collection
      ff <- self$data
      # Set the object
      ff[[what]] <- x
      if(copy){
        # If the object is an actual copy of the original one
        Settings$new(data = ff)
      } else {
        self$data <- ff
      }
      invisible()
    }
  ),

  # Any private entries
  private = list(
    finalize = function() {
    }
  )
)
