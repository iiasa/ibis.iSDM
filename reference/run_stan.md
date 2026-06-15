# Fit a cmdstanr model

This function fits a stan model using the light-weight interface
provided by cmdstanr. The code was adapted from McElreath rethinking
package.

## Usage

``` r
run_stan(
  model_code,
  data = list(),
  algorithm = "sampling",
  chains = 4,
  cores = getOption("ibis.nthread"),
  threads = 1,
  iter = 1000,
  warmup = floor(iter/2),
  control = list(adapt_delta = 0.95),
  cpp_options = list(),
  force = FALSE,
  path = base::getwd(),
  save_warmup = TRUE,
  return_stanfit = FALSE,
  ...
)
```

## Arguments

- model_code:

  A [`character`](https://rdrr.io/r/base/character.html) pointing to the
  stan modelling code.

- data:

  A [`list`](https://rdrr.io/r/base/list.html) with all the parameters
  required to run the model_code in stan.

- algorithm:

  A [`character`](https://rdrr.io/r/base/character.html) giving the
  algorithm to use. Either `'sampling'` (Default), `'optimize'` or
  `'variational'` for penalized likelihood estimation.

- chains:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) indicating the
  number of chains to use for estimation.

- cores:

  Number of threads for sampling. Default set to
  `'getOption("ibis.nthread")'`. See
  [`ibis_options()`](https://iiasa.github.io/ibis.iSDM/reference/ibis_options.md).

- threads:

  [`numeric`](https://rdrr.io/r/base/numeric.html) giving the number of
  threads to be run per chain. Has to be specified in accordance with
  cores.

- iter:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) value giving the
  number of MCMC samples to generate.

- warmup:

  [`numeric`](https://rdrr.io/r/base/numeric.html) for the number of
  warm-up samples for MCMC. Default set to 1/2 of iter.

- control:

  A [`list`](https://rdrr.io/r/base/list.html) with further control
  options for stan.

- cpp_options:

  A [`list`](https://rdrr.io/r/base/list.html) with options for the Cpp
  compiling.

- force:

  [`logical`](https://rdrr.io/r/base/logical.html) indication whether to
  force recompile the model (Default: `FALSE`).

- path:

  [`character`](https://rdrr.io/r/base/character.html) indicating a path
  to be made available to the stan compiler.

- save_warmup:

  A [`logical`](https://rdrr.io/r/base/logical.html) flag whether to
  save the warmup samples.

- return_stanfit:

  A [`logical`](https://rdrr.io/r/base/logical.html) flag whether to
  convert sampling output to an
  [`rstan`](https://mc-stan.org/rstan/reference/rstan.html) stanfit
  object. Defaults to `FALSE`; the native cmdstanr CmdStanFit object is
  used otherwise.

- ...:

  Other non-specified parameters.

## Value

A cmdstanr object by default, or a rstan object when requested and
conversion succeeds.

## See also

rethinking R package
