#' @include bdproto-engine.R bdproto-distributionmodel.R
NULL
#' Use STAN as engine
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param ... Other variables
#' @name engine_stan
NULL
#' @rdname engine_stan
#' @export
engine_stan <- function(x,
                        ...) {
  # Check whether INLA package is available
  check_package('rstan')
  if(!isNamespaceLoaded("rstan")) { attachNamespace("rstan");requireNamespace('rstan') }
  # assert that arguments are valid
  assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                          inherits(x$background,'sf')
  )
  # Use stan directly or rstanarm for now (latter impossible to specify spatial latent effects yet, but see  https://github.com/stan-dev/rstanarm/issues/207)
  # https://github.com/ConnorDonegan/Stan-IAR
  # BRMS joint modelling?
  # rstanarm object to improve speed through approximation
  # Set engine in distribution object
  x$set_engine(
    bdproto(
      "STAN-Engine",
      Engine,
      name = "<STAN>",
      data = list(),
      # Setup a model
      setup = function(self){
        # Parallel processing
        options(mc.cores = parallel::detectCores())
      },
      train = function(self){

      }
    )
  ) # End of engine definition
}
