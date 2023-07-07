#' @include utils.R bdproto.R
NULL

#' @export
if (!methods::isClass("PriorList")) methods::setOldClass("PriorList")
NULL

#' List of Priors supplied to an class
#'
#' This class represents a collection of [`Prior-class`] objects.
#' It provides methods for accessing, adding and removing priors from the list
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
#' @family bdproto
#' @keywords bdproto
NULL

#' @export
PriorList <- bdproto(
  "PriorList",
  priors = new_waiver(),
  # Print out summary statistics
  print = function(self) {
    if(is.Waiver(self$priors)) { message('No priors found.') }
      else{
        # Get priors for variables
        message('Set priors: ', text_green( length(self)) )
      }
  },
  show = function(self) {
    self$print()
  },
  # Number of priors
  length = function(self){
    if(is.Waiver(self$priors)) return(0) else base::length(self$priors)
  },
  # Return ids
  ids = function(self) {
    if(is.Waiver(self$priors)) return( new_waiver() )
    vapply(self$priors, function(x) x$id, character(1) )
  },
  # Variable names
  varnames = function(self) {
    if(is.Waiver(self$priors)) return( character(0) )
    vapply(self$priors, function(x) x$variable, character(1) )
  },
  # Function to return the classes of all contained priors
  classes = function(self){
    if(is.Waiver(self$priors)) return( character(0) )
    vapply(self$priors, function(x) class(x)[1], character(1) )
  },
  # Get types of all contained priors
  types = function(self){
    if(is.Waiver(self$priors)) return( character(0) )
    if(!any(sapply(self$priors, function(x) base::length(x$type)>0))) return( character(0) ) # No types set
    vapply(self$priors, function(x) x$type, character(1) )
  },
  # Exists a specific prior for a variable and type combination?
  exists = function(self, variable, type = NULL){
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
  # Add a new prior
  add = function(self, p) {
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
  # Get specific prior values from the list if set
  get = function(self, variable, type = NULL, what = "value"){
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
  # Collect priors from ids
  collect = function(self, id){
    assertthat::assert_that(!missing(id) )
    id = as.character(id)
    assertthat::assert_that( all(id %in% as.character(self$ids())) )
    priors( self$priors[id] )
  },
  # Remove a set prior by id
  rm = function(self, id){
    assertthat::assert_that(is.Id(id) || is.character(id),
                            id %in% self$ids()
                            )
    self$priors[[id]] <- NULL
    invisible()
  },
  # Summary function that lists all priors
  summary = function(self){
    if(is.Waiver(self$priors)) return(NULL)
    # Loop through each id and summarize the hyper parameters
    ids <- self$ids()
    varn <- self$varnames()
    out <- data.frame(id = character(), class = character(), variable = character(),
                      type = character(), hyper = character())
    for(id in ids){
      pp <- self$collect(ids)
      co <- data.frame(id = as.character( id ),
                       class = self$classes()[[id]],
                       variable = varn[[id]],
                       type = self$types()[id],
                       hyper = paste0( self$get(varn[id]) ,collapse = " | ") )
      out <- rbind(out, co)
    }
    assertthat::assert_that(nrow(out)>0)
    return(out)
  },
  # Combining function to combine this PriorList with another
  combine = function(self, x){
    assertthat::assert_that(inherits(x, 'PriorList'))
    if(is.Waiver(self$priors)){
      self$priors <- x$priors
    } else {
      # Check whether priors on variable and type already exist. If yes, replace
      for(p in x$priors) self$add(p)
    }
  }
  # TODO: Plotting function. Plots the distribution of all priors
  # plot = function(self){
  #   # Queries included prior objects and iteratively calls their plot function
  # }
)
