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
    #' @field engine The class name of the engine.
    #' @field name The name of the engine
    #' @field data Any data or parameters necessary to make this engine work.
    engine = character(0),
    name = character(0),
    data = list(),

    #' @description
    #' Initializes the object and creates an empty list
    #' @param engine The class name of the engine.
    #' @param name The name of the engine
    #' @return NULL
    initialize = function(engine, name){
      assertthat::assert_that(
        is.character(engine),
        is.character(name)
      )
      if(missing(name)){
        self$name <- text_red("<Unspecified>")
      } else {
        self$name <- name
      }
      self$engine <- engine
    },

    #' @description
    #' Print the Engine name
    #' @return A message on screen
    print = function() {
      message( paste0('Engine: ', self$name) )
    },

    #' @description
    #' Aliases that calls print.
    #' @return A message on screen
    show = function() {
      self$name
    },

    #' @description
    #' Get class description
    #' @return A [`character`] with the class as saved in engine
    get_class = function(){
      return( self$engine )
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
      invisible(self)
    },

    #' @description
    #' Dummy function to get self object
    #' @return This object
    get_self = function(){
      return( self )
    }
  ),

  # Any private entries
  private = list(
    finalize = function() {
    }
  ),
  # Don't lock class this will be altered in engine specific calls
  lock_objects = FALSE,
  lock_class = FALSE
)
