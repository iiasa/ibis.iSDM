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
                        warmup = floor(iter/2),
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
  # Use rstanarm for now (latter impossible to specify spatial latent effects yet, but see https://github.com/stan-dev/rstanarm/issues/207)
  # https://github.com/ConnorDonegan/Stan-IAR
  # rstanarm object to improve speed through approximation
  check_package('rstanarm')
  if(!isNamespaceLoaded("rstanarm")) { attachNamespace("rstanarm");requireNamespace('rstanarm') }
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

  if(is.null(cores)) cores <- getOption('ibis.nthread')

  # Create a background raster
  if(is.Waiver(x$predictors)){
    # Create from background
    template <- raster::raster(
      ext = raster::extent(x$background),
      crs = raster::projection(x$background),
      res = c(diff( (sf::st_bbox(x$background)[c(1,3)]) ) / 100, # Simplified assumption for resolution
              diff( (sf::st_bbox(x$background)[c(1,3)]) ) / 100
      )
    )
  } else {
    # If predictor existing, use them
    template <- emptyraster(x$predictors$get_data() )
  }
  # Burn in the background
  template <- raster::rasterize(x$background, template, field = 0)

  # Set engine in distribution object
  x$set_engine(
    bdproto(
      "STAN-Engine",
      Engine,
      name = "<STAN>",
      data = list(
        'template' = template
      ),
      # Stan options
      stan_param = list(
        chains = chains, iter = iter,
        warmup = warmup, init = init,
        cores = cores, control = control
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
          cores = cores, control = control
        )

      },
      # Spatial latent effect
      get_equation_latent_spatial = function(self){
        stop("TBD")
      },
      # Setup a model
      setup = function(self, model, ...){
        # Parallel processing
        # Simple security checks
        assertthat::assert_that(
          assertthat::has_name(model, 'background'),
          assertthat::has_name(model, 'biodiversity'),
          # Check that all predictors are present
          all(all.vars(model$biodiversity[[1]]$equation)[-1] %in% names(model$biodiversity[[1]]$predictors)),
          nrow(model$predictors) == ncell(self$get_data('template')),
          length(model$biodiversity) == 1 # Only works with single likelihood. To be processed separately
        )
        # Detect and format the family
        fam <- model$biodiversity[[1]]$family
        if( fam == 'poisson' ) {fam <- poisson } else if( fam == 'binomial') fam <- binomial(link = "cloglog")
        self$data[['family']] <- fam

        # Check that all stan parameters are appropriately set
        assertthat::assert_that(
          is.numeric(self$stan_param$chains),
          is.numeric(self$stan_param$iter),
          is.numeric(self$stan_param$warmup),
          is.numeric(self$stan_param$init)
        )

        # Define an algorithm for MCMC sampling
        # an be "sampling" for MCMC (the default), "optimizing" for optimization,
        # "meanfield" for variational inference with independent normal distributions,
        # or "fullrank" for variational inference with a multivariate normal distribution.
        settings$set('algorithm', 'optimizing')

        # Set a model seed for reproducibility
        settings$set('seed', 31337)

        # Messenger
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Engine setup.')

        invisible()
      },
      train = function(self, model, settings, ...){
        # Check assumptions

        # Collect data
        equation <- model$biodiversity[[1]]$equation
        data <- cbind(model$biodiversity[[1]]$predictors,
                      data.frame(observed = model$biodiversity[[1]]$observations[,'observed']) )
        # Weights
        w <- model$biodiversity[[1]]$expect
        if(model$biodiversity[[1]]$family=='binomial') w <- NULL # Set weights to 0 when binomial
        # Offsets
        if('offset' %in% names(model$biodiversity[[1]]) ){
          # Add offset to full prediction and load vector
          stop("Needs work")
          off <- model$biodiversity[[1]]$offset[, names(model$offset)[3] ]
        } else { off = NULL }


        # Model estimation
        try({
          fit_stan <- rstanarm::stan_glm(
            equation, # The equation
            data = data, # The data for estimation
            weights = w, # The weights
            family = fam, # The family
            offset = off, # any set offset
            prior = rstanarm::lasso(1, 0), # Lasso prior
            prior_intercept = normal(0, 5),
            # MCMC options
            algorithm = settings$get('algorithm'),
            seed = settings$get('seed'),
            chains = self$stan_param$chains,
            iter = self$stan_param$iter,
            cores = self$stan_param$cores,
            warmup = self$stan_param$warmup
          )
        }, silent = TRUE)
        if(class(fit_stan) == 'try-error') {stop("Model did not converge")}

        prior_summary(fit_stan)
        bayesplot::color_scheme_set("brightblue")
        plot(fit_stan)
        pp_check(fit_stan, plotfun = "hist", nreps = 5) # ?bayesplot::ppc_hist


      }
    )
  ) # End of engine definition
}
