#' @include utils.R bdproto.R bdproto-priorlist.R
NULL

#' @export
if (!methods::isClass("Prior")) methods::setOldClass("Prior")
NULL

#' Base Prior class
#'
#' This class sets up the base class for priors and is the base class that is
#' inherited by any specific prior.
#'
#' @name Prior-class
#' @aliases Prior
#' @family bdproto
#' @keywords bdproto
NULL

#' @export
Prior <- bdproto(
  "Prior",
  id = new_waiver(),
  name = character(0),
  type = new_waiver(),
  variable = character(0),
  distribution = new_waiver(),
  value = vector(), # A named vector containing the hyper-parameters
  # Print
  print = function(self) {
    message(paste0(
      class(self)[1], ': ', self$type, ' prior for \'', self$variable,'\''
    ))
  },
  # Generic validation function
  validate = function(self, x){
    assertthat::assert_that(
      all( is.numeric(x) ),
      all( !is.nan(x) ),
      all( !is.infinite(x) )
    )
  },
  # Get prior values
  get = function(self, what = "value"){
    assertthat::assert_that(what %in% names(self))
    return(self[[what]])
  },
  # Set prior
  set = function(self, x){
    self$value <- x
  }
)
