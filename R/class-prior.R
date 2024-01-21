if (!methods::isClass("Prior")) methods::setOldClass("Prior")

#' Base Prior class
#'
#' @description
#' This class sets up the base class for priors which will be inherited by all priors.
#'
#' @name Prior-class
#' @return Defines a Prior object.
#' @aliases Prior
#' @keywords classes
NULL

#' @rdname Prior-class
#' @export
Prior <- R6::R6Class(
  "Prior",
  public = list(
    #' @field id A [`character`] with the id of the prior.
    #' @field name A [`character`] with the name of the prior.
    #' @field type A [`character`] with the type of the prior.
    #' @field variable A [`character`] with the variable name for the prior.
    #' @field distribution A [`character`] with the distribution of the prior if relevant.
    #' @field value A [`numeric`] or [`character`] with the prior value, e.g. the hyper-parameters.
    id = new_waiver(),
    name = character(0),
    type = new_waiver(),
    variable = character(0),
    distribution = new_waiver(),
    value = vector(),

    #' @description
    #' Initializes the object
    #' @return NULL
    initialize = function(){
    },

    #' @description
    #' Print out the prior type and variable.
    #' @return A message on screen
    print = function() {
      message(paste0(
        class(self)[1], ': ', self$type, ' prior for \'', self$variable,'\''
      ))
    },

    #' @description
    #' Generic validation function for a provided value.
    #' @param x A new prior value.
    #' @note This functionality likely is deprecated or checks have been superseeded.
    #' @return Invisible TRUE
    validate = function(x){
      assertthat::assert_that(
        all( (is.numeric(x) || is.character(x)) ),
        all( !is.nan(x) ),
        all( !is.infinite(x) )
      )
      invisible()
    },

    #' @description
    #' Get prior values
    #' @param what A [`character`] with the entry to be returned (Default: \code{value}).
    #' @return Invisible TRUE
    get = function(what = "value"){
      assertthat::assert_that(what %in% names(self))
      return(self[[what]])
    },

    #' @description
    #' Set prior
    #' @param x A new prior value as [`numeric`] or [`character`].
    #' @return Invisible TRUE
    set = function(x){
      self$value <- x
      invisible()
    },

    #' @description
    #' Get a specific ID from a prior.
    #' @return A [`character`] id.
    get_id = function(){
      return( self$id )
    }
  ),

  # Any private entries
  private = list(
    finalize = function() {
    }
  )
)
