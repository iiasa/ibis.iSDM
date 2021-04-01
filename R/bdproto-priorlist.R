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
#' @examples{
#'  priorlist(
#'     INLAPrior('var1',c(0,0.1)),
#'     INLAPrior('var2',c(0,0.1))
#'    )
#' }
#' @aliases PriorList
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
        message('Set priors: ',length(self))
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
    return(
      vapply(self$priors, function(x) x$id, character(1) )
            )
  },
  # Variable names
  varnames = function(self) {
    if(is.Waiver(self$priors)) return( character(0) )
    vapply(self$priors, function(x) x$variable, character(1) )
  },
  # Exists a specific prior for a variable and type combination?
  exists = function(self, variable, type){
    assertthat::assert_that(!missing(variable), !missing(type),
                            is.character(variable),
                            is.character(type),msg = 'Specify both variable and type to query a prior.')
    vn <- self$varnames()
    ty <- vapply(self$priors, function(x) x$type, character(1) )
    # Return the id of the combination
    if((variable %in% vn) && (type %in% ty) ) return( names(vn)[which(variable %in% vn)] ) else return(NULL)
  },
  # Get specific prior values from the list if set
  get = function(self, variable, type){
    assertthat::assert_that(!missing(variable), !missing(type),
                            is.character(variable),
                            is.character(type),msg = 'Specify both variable and type to query a prior.')
    # Check whether there is an id
    ex <- self$exists(variable, type)
    # Catch and return NULL in case if not set
    if(is.null(ex)) ex else self$priors[[ex]]$get()
  },
  # Add a new prior
  add = function(self, x) {
    assertthat::assert_that(inherits(x, "Prior"))

    # If variable and type already exist, replace
    ex <- self$exists(x$variable, x$type)
    # Catch and return NULL in case if not set
    if(is.null(ex)) {
      # Set prior
      self$priors[[as.character(x$id)]] <- x
    } else {
      # Otherwise delete previous prior and add new one
      self$rm(ex)
      self$priors[[as.character(x$id)]] <- x
    }
    invisible()
  },
  # Remove a set prior by id
  rm = function(self, id){
    assertthat::assert_that(is.Id(id) || is.character(id),
                            id %in% self$ids()
                            )
    x$priors[[id]] <- NULL
    invisible()
  },
  # Combining function to combine this PriorList with another
  combine = function(self, x){
    assertthat::assert_that(inherits(x, 'PriorList'))
    if(is.Waiver(self$priors)) self$priors <- x$priors
    # Check whether priors on variable and type already exist. If yes, replace
    for(p in x$priors) self$add(p)
  }
  # TODO: Plotting function. Plots the distribution of all priors
  # plot = function(self){
  #   # Queries included prior objects and iteratively calls their plot function
  # }
)
