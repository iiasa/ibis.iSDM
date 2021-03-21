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
  # Add a new prior
  add = function(self, x) {
    assertthat::assert_that(inherits(x, "Prior"))
    self$priors[[as.character(x$id)]] <- x
    invisible()
  },
  # Remove a set prior
  rm = function(self, x){
    assertthat::assert_that(is.Id(x) || is.character(x))
  }
  # TODO: Plotting function. Plots the distribution of all priors
  # plot = function(self){
  #   # Queries included prior objects and iteratively calls their plot function
  # }
)
