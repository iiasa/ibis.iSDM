#' Checks whether cmdstanr is available and otherwise tries to install it
#'
#' @param install A [`logical`] factor to indicate whether [cmdstanr] should be directly installed (Default: TRUE)
#' @param ask [`logical`] whether the cmdstanr package is to be installed (Default: FALSE)
#' @keywords stan, utils
stan_check_cmd <- function(install = TRUE, ask = FALSE){
  assertthat::assert_that(
    is.logical(install), is.logical(ask)
  )
  # Check if available
  if(!requireNamespace("cmdstanr", quietly = TRUE)){
    if(install){
      if(ask){ a <- askYesNo("Install cmdstanr?") } else { a <- TRUE}
      if(a){
        install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
        check_cmdstan_toolchain()
        install_cmdstan(cores = 2)
      }
    } else {
      check_package("cmdstanr")
    }
  } else {
    invisible()
  }
}

#' Write a cmdstanr model output to a specific file
#'
#' @description Write a [cmdstanr] model output to a specific destination
#' @param mod A supplied [cmdstanr] model
#' @param dir The model directory where the model chould be written. Should be a character / existing dir.
#' @keywords stan, utils
write_stanmodel <- function( mod, dir = tempdir() ) {
  assertthat::assert_that(
    dir.exists(dir)
  )
  fname <- file.path( dir , concat("rt_cmdstanr_", digest::digest(mod,"md5")) )
  file_stan <- concat( fname, ".stan" )
  fileConn <- file( file_stan )
  writeLines( mod , fileConn )
  close(fileConn)
  return(file_stan)
}

#' Fit [cmdstanr] model and convert to [rstan] object
#'
#' @description This function fits a stan model using the light-weight cmdstanr interface. Code
#' was adapted from McElreath rethinking package.
#' @param model_code A [`character`] pointing to the stan modelling code
#' @param data A [`list`] with all the parameters required to run the [model_code] in stan.
#' @param algorithm A [`character`] giving the algorithm to use. Either 'sampling' (Default), 'optimize' or 'variational' for penalized likelihood estimation.
#' @param chains A [`numeric`] indicating the number of chains to use for estimation.
#' @param cores Number of threads for sampling. Default set to 'getOption("ibis.nthread")'. See [ibis_options()]
#' @param threads [`numeric`] giving the number of threads to be run per chain. Has to be specified in accordance with cores.
#' @param iter A [`numeric`] value giving the number of MCMC samples to generate.
#' @param warmup [`numeric`] for the number of warm-up samples for MCMC. Default set to 1/2 of iter.
#' @param control A [`list`] with further control options for stan.
#' @param cpp_options A [`list`] with options for the Cpp compiling.
#' @param save_warmup A [`logical`] flag whether to save the warmup samples.
#' @param ... Other non-specified parameters
#' @seealso [rethinking] R package
#' @returns A [rstan] object
#' @keywords misc, stan
run_stan <- function( model_code, data = list(),
                      algorithm = "sampling",
                      chains = 4, cores = getOption("ibis.nthread"),
                      threads = 1,
                      iter = 1000, warmup = floor(iter / 2),
                      control = list(adapt_delta=0.95),
                      cpp_options = list(),
                      save_warmup = TRUE, ... ) {
  assertthat::assert_that(
    is.numeric(chains), is.numeric(cores),
    is.numeric(iter), is.numeric(warmup),
    is.numeric(threads),
    threads < cores,
    is.list(control), is.list(cpp_options),
    is.logical(save_warmup)
  )
  # Check that cmdstanr is available
  check_package("cmdstanr")
  cmdstanr::check_cmdstan_toolchain(quiet = TRUE)

  # Match the algorithm to be used
  algorithm <- match.arg(algorithm, c("sampling", "optimize", "variational"), several.ok = FALSE)

  if( threads > 1 ) cpp_options[['stan_threads']] <- TRUE

  # Check extension
  assertthat::assert_that(
    is.character(model_code),
    assertthat::has_extension(model_code, "stan")
  )

  # Now compile the model
  mod <- cmdstanr::cmdstan_model( model_code,
                                  compile = TRUE,
                                  cpp_options = cpp_options )

  # Final parameters for sampling
  samp <- iter - warmup
  warm <- warmup

  # pull out any control arguments
  carg_adapt_delta <- 0.95
  if ( !is.null( control[['adapt_delta']] ) )
    carg_adapt_delta <- as.numeric(control[['adapt_delta']])
  carg_max_treedepth <- 11
  if ( !is.null( control[['max_treedepth']] ) )
    carg_max_treedepth <- as.numeric(control[['max_treedepth']])

  if(algorithm == "sampling"){
    # Sample
    if ( threads > 1 ) {
      cmdstanfit <- mod$sample( data = data,
                                chains = chains,
                                parallel_chains = cores,
                                iter_sampling = samp, iter_warmup = warm,
                                adapt_delta = carg_adapt_delta,
                                max_treedepth = carg_max_treedepth,
                                threads_per_chain = threads,
                                save_warmup = save_warmup,
                                ... )
      # coerce to stanfit object
      stanfit <- rstan::read_stan_csv( cmdstanfit$output_files() )

    } else {
      cmdstanfit <- mod$sample( data = data,
                                chains = chains,
                                parallel_chains = cores,
                                iter_sampling = samp , iter_warmup = warm,
                                adapt_delta = carg_adapt_delta,
                                max_treedepth = carg_max_treedepth,
                                save_warmup = save_warmup,
                                ... )
    }
    # coerce to stanfit object
    stanfit <- rstan::read_stan_csv( cmdstanfit$output_files() )

  } else if(algorithm == "optimize"){
    # Optimize for getting point estimates
    stanfit <- mod$optimize(data = data,
                               #seed = seed, # This could be passed on
                               threads = threads
                               )
  } else if(algorithm == "variational") {
    # Variational for approximating the posterior
    stanfit <- mod$variational(data = data,
                                  # seed = seed,
                                  threads = threads
    )
  }

  return(stanfit)
}
