#' @include bdproto-engine.R bdproto-distributionmodel.R
NULL
#' Use STAN as engine
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param chains A positive integer specifying the number of Markov chains. The default is 4.
#' @param iter A positive integer specifying the number of iterations for each chain (including warmup). The default is 2000.
#' @param warmup positive integer specifying the number of warmup (aka burnin) iterations per chain. If step-size adaptation is on (which it is by default), this also controls the number of iterations for which adaptation is run (and hence these warmup samples should not be used for inference). The number of warmup iterations should be smaller than iter and the default is iter/2.
#' @param cores If set to NULL take values from specified ibis option getOption('ibis.nthread')
#' @param init Initial values for parameters. Default: 'random'. Can also be specified as list (see: rstan::stan)
#' @param control See rstan::stan for more details on specifying the controls
#' @param ... Other variables
#' @seealso rstan
#' @name engine_stan
NULL
#' @rdname engine_stan
#' @export
engine_stan <- function(x,
                        chains = 4,
                        iter = 2000,
                        warmup = 500,
                        init = "random",
                        cores = NULL,
                        control = NULL,
                        ...) {
  # Check whether INLA package is available
  check_package('rstan')
  if(!isNamespaceLoaded("rstan")) { attachNamespace("rstan");requireNamespace('rstan') }
  # assert that arguments are valid
  assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                          inherits(x$background,'sf')
  )
  # Other checks of parameters
  assertthat::assert_that(
    is.numeric(chains), is.numeric(iter), is.numeric(warmup),
    is.null(cores) || is.numeric(cores),
    is.character(init) || is.list(init),
    is.null(control) || is.list(control),
    msg = 'Input parameters wrongly specified!'
  )
  # CHECK:
  # Use stan directly or rstanarm for now (latter impossible to specify spatial latent effects yet, but see https://github.com/stan-dev/rstanarm/issues/207)
  # https://github.com/ConnorDonegan/Stan-IAR
  # rstanarm object to improve speed through approximation

  if(is.null(cores)) cores <- getOption('ibis.nthread')

  # Set engine in distribution object
  x$set_engine(
    bdproto(
      "STAN-Engine",
      Engine,
      name = "<STAN>",
      data = list(),
      # Stan options
      stan_param = list(
        chains = chains, iter = iter,
        warmup = warmup, init = init,
        cores = cores, control = control,
      ),
      # Function to respecify the control parameters
      set_control = function(self,
                             chains = 4,
                             iter = 2000,
                             warmup = 500,
                             init = "random",
                             cores = NULL,
                             control = NULL){

        # Overwrite existing
        self$stan_param <- list(
          chains = chains, iter = iter,
          warmup = warmup, init = init,
          cores = cores, control = control,
        )

      },
      # Setup a model
      setup = function(self){
        # Parallel processing

      },
      train = function(self){
      }
    )
  ) # End of engine definition
}
