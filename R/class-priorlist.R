#' @include class-prior.R
NULL
if (!methods::isClass("PriorList")) methods::setOldClass("PriorList")

#' List of Priors supplied to an class
#'
#' @description
#' This class represents a collection of [`Prior-class`] objects. It provides
#' methods for accessing, adding and removing priors from the list
#'
#' @name PriorList-class
#' @examples
#' \dontrun{
#' priors(
#'     INLAPrior('var1','normal',c(0,0.1)),
#'     INLAPrior('var2','normal',c(0,0.1))
#'    )
#' }
#' @return A PriorList object.
#' @aliases PriorList
#' @keywords classes
NULL

#' @rdname PriorList-class
#' @export
PriorList <- R6::R6Class(
  "PriorList",
  public = list(
    #' @field priors A list of [`Prior-class`] object.
    priors = new_waiver(),

    #' @description
    #' Initializes the object
    #' @return NULL
    initialize = function(){
    },

    #' @description
    #' Print out summary statistics
    #' @return A message on screen
    print = function() {
      if(is.Waiver(self$priors)) { message('No priors found.') }
        else{
          # Get priors for variables
          message('Set priors: ', text_green( length(self)) )
        }
    },

    #' @description
    #' Aliases that calls print.
    #' @return A message on screen
    show = function() {
      self$print()
    },

    #' @description
    #' Number of priors in object
    #' @return A [`numeric`] with the number of priors set
    length = function(){
      if(is.Waiver(self$priors)) return(0) else base::length(self$priors)
    },

    #' @description
    #' Ids of prior objects
    #' @return A list with ids of the priors objects for query
    ids = function() {
      if(is.Waiver(self$priors)) return( new_waiver() )
      vapply(self$priors, function(x) x$id, character(1) )
    },

    #' @description
    #' Variable names of priors in object
    #' @return A [`character`] list with the variable names of the priors.
    varnames = function() {
      if(is.Waiver(self$priors)) return( character(0) )
      vapply(self$priors, function(x) x$variable, character(1) )
    },

    #' @description
    #' Function to return the classes of all contained priors
    #' @return A [`character`] list with the class names of the priors.
    classes = function(){
      if(is.Waiver(self$priors)) return( character(0) )
      vapply(self$priors, function(x) class(x)[1], character(1) )
    },

    #' @description
    #' Get types of all contained priors
    #' @return A [`character`] list with the type names of the priors.
    types = function(){
      if(is.Waiver(self$priors)) return( character(0) )
      if(!any(sapply(self$priors, function(x) base::length(x$type)>0))) return( character(0) ) # No types set
      vapply(self$priors, function(x) x$type, character(1) )
    },

    #' @description
    #' Does a certain variable or type combination exist as prior ?
    #' @param variable A [`character`] with the variable name.
    #' @param type A [`character`] with the type.
    #' @return A [`character`] id.
    exists = function(variable, type = NULL){
      assertthat::assert_that(!missing(variable),
                              is.character(variable),
                              is.character(type) || is.null(type) || is.Waiver(type),
                              msg = 'Specify a single variable to query a prior.')
      # overwrite if not set
      if(is.Waiver(type) || is.null(type)) {type <- NULL}
      vn <- self$varnames()
      vt <- self$types()
      # If type is specified
      if(!is.null(type)){
        # Check whether variable is in vn
        if(any(variable %in% vn)){
          # Check type and return the id of the combination
          if(type %in% vt){
            id <- names(vn)[which(variable == as.character(vn) & type == as.character(vt))]
          } else {
            # Get type of variable instead
            id <- names(vn)[which(variable %in% as.character(vn))]
          }
        } else { id <- NULL }
      } else {
        # Simply match against variable names and return id
        if(!all(variable %in% vn)){
          id <- NULL
        } else {
          id <- names(vn)[base::match(variable, vn, nomatch = 0)]
        }
      }
      return(id)
    },

    #' @description
    #' Add a new prior to the object.
    #' @param p A [`Prior-class`] object.
    #' @return Invisible TRUE
    add = function(p) {
      assertthat::assert_that(inherits(p, "Prior"))

      # Check if variable and type or variable exists already?
      ex <- self$exists(p$variable, p$type)
      # In case there is not more one variable (SPDE), also check without type
      if(sum(p$variable == self$varnames())<2){
        ex2 <- self$exists(p$variable, NULL)
      } else { ex2 <- "" }

      # Catch and return NULL in case if not set
      if(is.null(ex)) {
        # Set prior
        self$priors[[as.character(p$id)]] <- p
      } else {
        # Otherwise first delete previous prior and add new one
        if(nchar(ex2)>1){ self$rm(ex2)} else {self$rm(ex)}
        self$priors[[as.character(p$id)]] <- p
      }
      invisible()
    },

    #' @description
    #' Get specific prior values from the list if set
    #' @param variable A [`character`] with the variable name.
    #' @param type A [`character`] with the type name
    #' @param what A [`character`] on the specific entry to return (Default: \code{prior value}).
    #' @return The prior object.
    get = function(variable, type = NULL, what = "value"){
      assertthat::assert_that(!missing(variable),
                              is.character(variable),
                              is.character(type) || is.null(type),
                              is.character(what),
                              msg = 'Specify a variable to query a prior.')
      # Check whether there is an id
      ex <- self$exists(variable, type)
      # Catch and return NULL in case if not set
      if(is.null(ex)) ex else {
        # Get all values
        vals <- lapply(self$priors, function(x) x$get(what) )
        vals[[ex]]
      }
    },

    #' @description
    #' Collect priors for a given id or multiple.
    #' @param id A [`character`] with the prior id.
    #' @return A [`PriorList-class`] object.
    collect = function(id){
      assertthat::assert_that(!missing(id) )
      id = as.character(id)
      assertthat::assert_that( all(id %in% as.character(self$ids())) )
      priors( self$priors[id] )
    },

    #' @description
    #' Remove a set prior by id
    #' @param id A [`character`] with the prior id.
    #' @return Invisible TRUE
    rm = function(id){
      assertthat::assert_that(is.Id(id) || is.character(id),
                              id %in% self$ids()
                              )
      self$priors[[id]] <- NULL
      invisible()
    },

    #' @description
    #' Summary function that lists all priors
    #' @return A [`data.frame`] with the summarized priors.
    summary = function(){
      if(is.Waiver(self$priors)) return(NULL)
      # Loop through each id and summarize the hyper parameters
      ids <- self$ids()
      varn <- self$varnames()
      out <- data.frame(id = character(), class = character(), variable = character(),
                        type = character(), hyper = character(), other = character())
      for(id in as.character(ids)){
        pp <- self$collect(ids)
        co <- data.frame(id = as.character( id ),
                         class = self$classes()[[id]],
                         variable = varn[[id]],
                         type = self$types()[id],
                         hyper = paste0( self$get(varn[id]) ,collapse = " | ") )
        # Add other things if found
        for(v in c("lims", "prop")){
          lim <- try({self$get(varn[id], what = v)},silent = TRUE)
          if(!inherits(lim, "try-error")){
            co$other <- paste0(lim, collapse = " | ")
          }
        }
        out <- rbind(out, co)
      }
      assertthat::assert_that(nrow(out)>0)
      return(out)
    },

    #' @description
    #' Combining function to combine this PriorList with another new one
    #' @param x A new [`PriorList-class`] object.
    #' @return Invisible TRUE
    combine = function(x){
      assertthat::assert_that(inherits(x, 'PriorList'))
      if(is.Waiver(self$priors)){
        self$priors <- x$priors
      } else {
        # Check whether priors on variable and type already exist. If yes, replace
        for(p in x$priors) self$add(p)
      }
      invisible()
    }

    # TODO: Plotting function. Plots the distribution of all priors
    # plot = function(self){
    #   # Queries included prior objects and iteratively calls their plot function
    # }
  ),

  # Any private entries
  private = list(
    finalize = function() {
    }
  )
)
