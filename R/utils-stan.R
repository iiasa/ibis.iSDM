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

#' Fit cmdstanr model and convert to rstan
#'
#' @param file
#' @seealso [rethinking] R package
#' @returns A [rstan] object
#' @keywords misc, stan
fit_stan <- function( file , model_code , data=list(),
                      chains = 1, cores = 1,
                      iter = 1000, warmup,
                      threads = 1,
                      control = list(adapt_delta=0.95),
                      cpp_options = list(),
                      save_warmup = TRUE, ... ) {
  assertthat::assert_that(
    is.numeric(chains), is.numeric(cores),
    is.numeric(iter), is.numeric(warmup),
    is.numeric(threads),
    is.list(control), is.list(cpp_options),
    is.logical(save_warmup)
  )
  # Check that cmdstanr is available
  check_package("cmdstanr")

  if ( threads>1 ) cpp_options[['stan_threads']] <- TRUE

  if ( missing(file) & !missing(model_code) ) {
    file <- cmdstanr_model_write( model_code )
  }



  mod <- cmdstan_model( file , compile=TRUE , cpp_options=cpp_options )

  if ( missing(warmup) ) {
    samp <- floor(iter/2)
    warm <- floor(iter/2)
  } else {
    samp <- iter - warmup
    warm <- warmup
  }

  # pull out any control arguments
  carg_adapt_delta <- 0.95
  if ( !is.null( control[['adapt_delta']] ) )
    carg_adapt_delta <- as.numeric(control[['adapt_delta']])
  carg_max_treedepth <- 11
  if ( !is.null( control[['max_treedepth']] ) )
    carg_max_treedepth <- as.numeric(control[['max_treedepth']])

  # sample
  if ( threads > 1 )
    cmdstanfit <- mod$sample( data=data ,
                              chains=chains ,
                              parallel_chains=cores ,
                              iter_sampling=samp , iter_warmup=warm ,
                              adapt_delta=carg_adapt_delta ,
                              max_treedepth=carg_max_treedepth ,
                              threads_per_chain=threads ,
                              save_warmup=save_warmup , ... )
  else
    cmdstanfit <- mod$sample( data=data ,
                              chains=chains ,
                              parallel_chains=cores ,
                              iter_sampling=samp , iter_warmup=warm ,
                              adapt_delta=carg_adapt_delta ,
                              max_treedepth=carg_max_treedepth ,
                              save_warmup=save_warmup , ... )

  # coerce to stanfit object
  stanfit <- rstan::read_stan_csv(cmdstanfit$output_files())

  return(stanfit)
}
