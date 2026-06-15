# Use Stan as engine

Stan is probabilistic programming language that can be used to specify
linear Bayesian species distribution models with one or more point
presence-only and presence-absence biodiversity datasets. Multiple
datasets are fitted jointly with shared covariate effects and
dataset-specific intercepts when requested via the biodiversity dataset
settings. **Requires the `"cmdstanr"` package to be installed!**

## Usage

``` r
engine_stan(
  x,
  chains = 4,
  iter = 2000,
  warmup = floor(iter/2),
  init = "random",
  cores = getOption("ibis.nthread"),
  algorithm = "sampling",
  control = list(adapt_delta = 0.95),
  type = "response",
  ...
)
```

## Arguments

- x:

  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

- chains:

  A positive [`integer`](https://rdrr.io/r/base/integer.html) specifying
  the number of Markov chains (Default: `4` chains).

- iter:

  A positive [`integer`](https://rdrr.io/r/base/integer.html) specifying
  the number of iterations for each chain (including warmup). (Default:
  `2000`).

- warmup:

  A positive [`integer`](https://rdrr.io/r/base/integer.html) specifying
  the number of warmup iterations per chain. The default is `iter/2`.

- init:

  Initial values for parameters (Default: `'random'`).

- cores:

  If set to NULL take values from specified ibis option
  `getOption('ibis.nthread')`.

- algorithm:

  Mode used to sample from the posterior. Available options are
  `"sampling"`, `"optimize"`, or `"variational"`.

- control:

  See `"cmdstanr"` for more details on sampler controls.

- type:

  The mode used for creating posterior predictions. Either `"predictor"`
  or `"response"` (Default: `"response"`).

- ...:

  Other variables.

## Value

An
[Engine](https://iiasa.github.io/ibis.iSDM/reference/Engine-class.md).

## See also

rstan, cmdstanr

Other engine:
[`engine_bart()`](https://iiasa.github.io/ibis.iSDM/reference/engine_bart.md),
[`engine_breg()`](https://iiasa.github.io/ibis.iSDM/reference/engine_breg.md),
[`engine_gdb()`](https://iiasa.github.io/ibis.iSDM/reference/engine_gdb.md),
[`engine_glm()`](https://iiasa.github.io/ibis.iSDM/reference/engine_glm.md),
[`engine_glmnet()`](https://iiasa.github.io/ibis.iSDM/reference/engine_glmnet.md),
[`engine_inla()`](https://iiasa.github.io/ibis.iSDM/reference/engine_inla.md),
[`engine_inlabru()`](https://iiasa.github.io/ibis.iSDM/reference/engine_inlabru.md),
[`engine_scampr()`](https://iiasa.github.io/ibis.iSDM/reference/engine_scampr.md),
[`engine_xgboost()`](https://iiasa.github.io/ibis.iSDM/reference/engine_xgboost.md)

## Examples

``` r
if (FALSE) { # \dontrun{
x <- distribution(background) |> engine_stan(iter = 1000)
} # }
```
