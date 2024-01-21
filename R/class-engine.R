if (!methods::isClass("Engine")) methods::setOldClass("Engine")

#' Engine class description
#'
#' @description
#' Basic object for engine, all other engines inherit from here.
#'
#' @keywords classes
#'
#' @name Engine-class
NULL

#' @rdname Engine-class
#' @export
Engine <- R6::R6Class(
  "Engine",
  public = list(
    #' @field name The name of the engine
    #' @field data Any data or parameters necessary to make this engine work.
    name = character(0),
    data = list(),

    #' @description
    #' Initializes the object and creates an empty list
    #' @return NULL
    initialize = function(){
      self$name <- text_red("<Unspecified>")
    },

    #' @description
    #' Print the Engine name
    #' @return A message on screen
    print = function() {
      message( paste0('Engine: ', self$name) )
    },

    #' @description
    #' Functions to be overwritten when engine is specified
    #' @param ... Parameter defined by the engine.
    #' @return Invisible
    setup = function(...) stop("Engine is missing functions to setup settings."),

    #' @description
    #' Functions to be overwritten when engine is specified
    #' @param ... Parameter defined by the engine.
    #' @return Invisible
    train = function(...) stop("Engine is missing a fitting method."),

    #' @description
    #' Aliases that calls print.
    #' @return A message on screen
    show = function() {
      self$name
    },

    #' @description
    #' Get specific data from this engine
    #' @param x A respecified data to be added to the engine.
    #' @return A [`list`] with the data.
    get_data = function(x) {
      if (!x %in% names(self$data))
        return(new_waiver())
      return(self$data[[x]])
    },

    #' @description
    #' List all data
    #' @return A [`character`] vector of the data entries.
    list_data = function(){
      return( names( self$data ) )
    },

    #' @description
    #' Set data for this engine
    #' @param x A [`character`] with the name or id of this dataset.
    #' @param value A new [`list`] of parameters.
    #' @return Invisible
    set_data = function(x, value) {
      self$data[[x]] <- value
      invisible()
    }
  ),

  # Any private entries
  private = list(
    finalize = function() {
    }
  )
)
