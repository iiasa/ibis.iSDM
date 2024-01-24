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
    #' @field prob Another [`numeric`] entry on the prior field. The inclusion probability.
    #' @field lims A limitation on the lower and upper bounds of a numeric value.
    id = new_waiver(),
    name = character(0),
    type = character(0),
    variable = character(0),
    distribution = new_waiver(),
    value = vector(),
    prob = numeric(0L),
    lims = vector(),

    #' @description
    #' Initializes the object and prepared the various prior variables
    #' @param id A [`character`] with the id of the prior.
    #' @param name A [`character`] with the name of the prior.
    #' @param type A [`character`] with the type of the prior.
    #' @param variable A [`character`] with the variable name for the prior.
    #' @param distribution A [`character`] with the distribution of the prior if relevant.
    #' @param value A [`numeric`] or [`character`] with the prior value, e.g. the hyper-parameters.
    #' @param prob Another [`numeric`] entry on the prior field. The inclusion probability.
    #' @param lims A limitation on the lower and upper bounds of a numeric value.
    #' @return NULL
    initialize = function(id, name, variable, value, type = NULL,
                          distribution = NULL, prob = NULL, lims = NULL){
      assertthat::assert_that(
        missing(id) || (is.Id(id) || is.character(id)),
        is.character(name),
        is.character(type) || is.null(type),
        is.character(variable),
        is.numeric(prob) || is.null(prob),
        is.null(lims) || is.numeric(lims)
      )
      # Assign a new id for the object
      if(missing(id)) id <- new_id()

      self$id <- id
      self$name <- name
      self$type <- type
      self$variable <- variable
      self$value <- value
      self$distribution <- distribution
      self$prob <- prob
      self$lims <- lims
    },

    #' @description
    #' Print out the prior type and variable.
    #' @return A message on screen
    print = function() {
      message(paste0(
        self$name, ': ', self$type, ' prior for \'', self$variable,'\''
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
      invisible(self)
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
      invisible(self)
    },

    #' @description
    #' Get a specific ID from a prior.
    #' @return A [`character`] id.
    get_id = function(){
      return( self$id )
    },

    #' @description
    #' Get Name of object
    #' @return Returns a [`character`] with the class name.
    get_name = function(){
      return( self$name )
    }

  ),

  # Any private entries
  private = list(
    finalize = function() {
    }
  ),
  # Allow custom bindings
  lock_objects = FALSE
)
