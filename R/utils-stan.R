#' Logistic (invlogit) transformation function
#' @param x A [`numeric`] value
#' @keywords utils
#' @noRd
logistic <- function(a){
  ilink(a, link = "logit")
}

#' Logit transformation function
#' @param x A [`numeric`] value
#' @keywords utils
#' @noRd
logit <- function(a){
  assertthat::assert_that(is.numeric(a),msg = 'Provided value not numeric')
  log(a/(1-a))
}

#' Inverse transformation function for the link function
#' @description back transforms a [numeric] vector using the appropriate link function
#' @param x A [`numeric`] vector generated by a model
#' @param link [`character`] indicating the link function to use. Default: "log"
#' @noRd
ilink <- function(x, link = "log"){
  assertthat::assert_that(is.numeric(x),
                          is.character(link)
                          )
  link <- match.arg(link, c("identity", "log","logm1","log1p",
                            "inverse", "sqrt", "logit", "probit",
                            "cauchit", "cloglog") , several.ok = FALSE)
  switch (link,
    identity = x,
    log = exp(x),
    logm1 = expp1(x),
    log1p = expm1(x),
    inverse = 1/x,
    sqrt = x^2,
    logit = ( 1 / (1 + exp(-x)) ),
    probit = pnorm(x),
    cauchit = pcauchy(x),
    cloglog = (1 - exp(-exp(x)))
  )
}

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

#' Wrap a list with stan model code
#'
#' @description [engine_stan] builds a list with stan model code. This function
#' concatenates them together
#' @param sm_code A [list] object with exactly 7 entries
#' @returns A [character] object
#' @keywords stan_utils
wrap_stanmodel <- function(sm_code){
  assertthat::assert_that(is.list(sm_code),
                          length(sm_code)==7)
  out <- character(0)

  # Functions
  out <- paste0("functions {")
  for(i in sm_code$functions) out <- paste0(out, i, "\n")
  out <- paste0(out, "\n}\n")
  # Data
  out <- paste(out, "data {")
  for(i in sm_code$data) out <- paste0(out, i, "\n")
  out <- paste0(out, "\n}\n")
  # Transformed data
  out <- paste(out, "transformed data {")
  for(i in sm_code$transformed_data) out <- paste0(out, i, "\n")
  out <- paste0(out, "\n}\n")
  # Parameters
  out <- paste(out, "parameters {")
  for(i in sm_code$parameters) out <- paste0(out, i, "\n")
  out <- paste0(out, "\n}\n")
  # Transformed parameters
  out <- paste(out, "transformed parameters {")
  for(i in sm_code$transformed_parameters) out <- paste0(out, i, "\n")
  out <- paste0(out, "\n}\n")
  # Model
  out <- paste(out, "model {")
  for(i in sm_code$model) out <- paste0(out, i, "\n")
  out <- paste0(out, "\n}\n")
  # Generated quantities
  out <- paste(out, "generated quantities {")
  for(i in sm_code$generated_quantities) out <- paste0(out, i, "\n")
  out <- paste0(out, "}")

  assertthat::assert_that(is.character(out), length(out)>0)
  return(out)
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
  fname <- file.path( dir , paste0("rt_cmdstanr_", digest::digest(mod,"md5")) )
  file_stan <- paste0( fname, ".stan" )
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
#' @param force [`logical`] indication whether to force recompile the model (default: False).
#' @param path [`character`] indicating a path to be made available to the stan compiler.
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
                      control = list(adapt_delta = 0.95),
                      cpp_options = list(),
                      force = FALSE,
                      path = base::getwd(),
                      save_warmup = TRUE, ... ) {
  assertthat::assert_that(
    is.numeric(chains), is.numeric(cores),
    is.numeric(iter), is.numeric(warmup),
    is.numeric(threads),
    threads < cores,
    is.list(control), is.list(cpp_options),
    is.logical(save_warmup),
    is.logical(force)
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
                                  force_recompile = force,
                                  cpp_options = cpp_options,
                                  include_paths = path
                                  )

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

#' Create a posterior prediction with rstanfit object
#'
#' @description This function does ...
#' @param obj A [`stanfit`] object (as used by [rstan]).
#' @param form A [`formula`] object created for the [ibis.iSDM::DistributionModel].
#' @param newdata A [data.frame] with new data to be used for prediction.
#' @param mode A [`character`] of whether the linear `predictor` or the `response` is to be summarized
#' @param family A [`character`] giving the family for simulating linear response values. Default: NULL
#' @param offset A [vector] with an optionally specified offset.
#' @param draws [numeric] indicating whether a specific number of draws should be taken.
#' @import posterior
#' @references [https://medium.com/@alex.pavlakis/making-predictions-from-stan-models-in-r-3e349dfac1ed](https://medium.com/@alex.pavlakis/making-predictions-from-stan-models-in-r-3e349dfac1ed)
#' @references The [brms] R-package
#' @export
posterior_predict_stanfit <- function(obj, form, newdata, mode = "predictor", family = NULL, offset = NULL, draws = NULL){
  assertthat::assert_that(
    inherits(obj, "stanfit") || inherits(obj, "CmdStanFit"),
    is.formula(form),
    is.data.frame(newdata),
    is.null(family) || is.character(family),
    is.null(draws) || is.numeric(draws),
    is.null(offset) || (length(offset) == nrow(newdata))
  )
  mode <- match.arg(mode, c("predictor", "response"), several.ok = FALSE)
  # Build model matrix
  # Note: This removes all NA cells from matrix
  A <- model.matrix(object = delete.response(terms(form)),
                     data = newdata)
  assertthat::assert_that(nrow(A)>0, inherits(A, "matrix") || inherits(A, "dgCMatrix"))
  # Remove intercept unless set
  # if(attr(terms(form),"intercept") != 1) {
    if(length(grep("Intercept", colnames(A), ignore.case = TRUE))>0){
      A <- A[,-(grep("Intercept", colnames(A), ignore.case = TRUE))]
    }
  # }

  # Draw from the posterior
  if(inherits(obj, "stanfit")) {
    pp <- posterior::as_draws_df(obj)
  } else {
    pp <- obj$draws() |> as.data.frame()
  }
  # Create a subset?
  if (!is.null(draws)) {
    pp <- pp[sample.int(nrow(pp), draws),]
  }
  # Get only beta coefficients and Intercept if set
  if("Intercept" %in% colnames(pp)) what <- "beta|Intercept" else what <- "beta"
  suppressWarnings( pp <- pp[ c(grep(what, colnames(pp), value = TRUE)) ] )

  # Prepare offset if set
  if(!is.null(offset)) {
    # Get only the rows in the A matrix (minus NA)
    offset <- offset[as.numeric(row.names(A))]
  } else { offset <- rep(0, nrow(A) ) }

  # Security checks
  assertthat::assert_that(
    nrow(A)>0, nrow(pp) > 0,
    ncol(pp) == ncol(A),
    is.numeric(offset)
  )

  if(mode == "predictor"){
    # Summarize the coefficients from the posterior
    pp <- posterior::summarise_draws(pp) |>
      subset(select = c("variable", "mean", "q5", "median", "q95", "sd"))  |>
      as.data.frame()
    # --- #
    pp$variable <- colnames(A)
    # Calculate b*X + offset if set
    preds <- cbind(
      A %*% pp[,"mean"] + ifelse(is.null(offset),0, offset),
      A %*% pp[,"q5"] + ifelse(is.null(offset),0, offset),
      A %*% pp[,"median"] + ifelse(is.null(offset),0, offset),
      A %*% pp[,"q95"] + ifelse(is.null(offset),0, offset),
      A %*% pp[,"sd"] + ifelse(is.null(offset),0, offset)
    )

    # Add random noise equivalent to the posterior length
    .rnorm_matrix <- function(mean, sd) {
      stopifnot(length(dim(mean)) == 2)
      error <- matrix(rnorm(length(mean), 0, sd), ncol=ncol(mean), byrow=TRUE)
      mean + error
    }
    preds <- .rnorm_matrix(preds, pp[,"sd"]) # FIXME: This only makes sense for mean. Apply mad to median?

    # Apply ilink
    if(!is.null(family)){
      preds <- switch (family,
        "poisson" = ilink(preds, link = "log"),
        "binomial" = ilink(preds, link = "logit"),
        ilink(preds, link = "log")
      )
    }

  } else {
    # Simulate linear response approximating poisson_rng in stan
    out <- vector("list", nrow(pp))
    # TODO: Parallelize over threads
    for(i in 1:nrow(pp)){
      # Build eta as additive beta with the A matrix row
      eta <- 0 + base::tcrossprod(as.matrix(pp)[i,] |> base::unname(), A) + offset
      out[[i]] <- base::unname(eta)
    }
    # ilink transform?
    # Get family and number of observations and
    # draw random variable for each draw and lambda value
    if(family == "poisson"){
      a <- lapply(out, function(lambda) rpois(nrow(A), ilink(lambda, link = "log")) )
    } else if(family == "binomial") {
      a <- lapply(out, function(mu) rbinom(nrow(A), size = 1, prob = ilink(mu, link = "logit")) )
    } else {
      stop("Not yet implemented method for linear response prediction")
    }
    # Combine link
    a <- do.call(cbind, a)
    colnames(a) <- rownames(a) <- NULL
    # Finally summarize
    preds <- cbind(
      matrixStats::rowMeans2(a, na.rm = TRUE),
      matrixStats::rowQuantiles(a, probs = c(.05,.5,.95), na.rm = TRUE),
      matrixStats::rowSds(a, na.rm = TRUE)
    )
  }
  # Create output with cellid
  out <- tibble::rowid_to_column(newdata,var = "cellid")["cellid"] |> as.data.frame()
  out$cv <- out$q95 <- out$q50 <- out$q05 <- out$sd <- out$mean <- NA
  out$mean[as.numeric(row.names(A))] <- preds[,1]
  out$sd[as.numeric(row.names(A))] <- preds[,5]
  out$q05[as.numeric(row.names(A))] <- preds[,2]
  out$q50[as.numeric(row.names(A))] <- preds[,3]
  out$q95[as.numeric(row.names(A))] <- preds[,4]
  out$cv[as.numeric(row.names(A))] <- preds[,5] / preds[,1]
  out$cellid <- NULL

  return(out)
}