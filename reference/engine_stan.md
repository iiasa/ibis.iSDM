# Use Stan as engine

Stan is probabilistic programming language that can be used to specify
most types of statistical linear and non-linear regression models. Stan
provides full Bayesian inference for continuous-variable models through
Markov chain Monte Carlo methods such as the No-U-Turn sampler, an
adaptive form of Hamiltonian Monte Carlo sampling. Stan code has to be
written separately and this function acts as compiler to build the
stan-model. **Requires the `"cmdstanr"` package to be installed!**

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
  the number of warmup (aka burnin) iterations per chain. If step-size
  adaptation is on (Default: `TRUE`), this also controls the number of
  iterations for which adaptation is run (and hence these warmup samples
  should not be used for inference). The number of warmup iterations
  should be smaller than `iter` and the default is `iter/2`.

- init:

  Initial values for parameters (Default: `'random'`). Can also be
  specified as [list](https://rdrr.io/r/base/list.html) (see:
  `"rstan::stan"`)

- cores:

  If set to NULL take values from specified ibis option
  `getOption('ibis.nthread')`.

- algorithm:

  Mode used to sample from the posterior. Available options are
  `"sampling"`, `"optimize"`, or `"variational"`. See `"cmdstanr"`
  package for more details. (Default: `"sampling"`).

- control:

  See `"rstan::stan"` for more details on specifying the controls.

- type:

  The mode used for creating posterior predictions. Either summarizing
  the linear `"predictor"` or `"response"` (Default: `"response"`).

- ...:

  Other variables

## Value

An
[Engine](https://iiasa.github.io/ibis.iSDM/reference/Engine-class.md).

## Details

By default the posterior is obtained through sampling, however stan also
supports approximate inference forms through penalized maximum
likelihood estimation (see Carpenter et al. 2017).

## Note

The function `obj$stancode()` can be used to print out the stancode of
the model.

## References

- Jonah Gabry and Rok Češnovar (2021). cmdstanr: R Interface to
  'CmdStan'. https://mc-stan.org/cmdstanr,
  https://discourse.mc-stan.org.

- Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B.,
  Betancourt, M., ... & Riddell, A. (2017). Stan: A probabilistic
  programming language. Journal of statistical software, 76(1), 1-32.

- Piironen, J., & Vehtari, A. (2017). Sparsity information and
  regularization in the horseshoe and other shrinkage priors. Electronic
  Journal of Statistics, 11(2), 5018-5051.

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
# Add Stan as an engine
x <- distribution(background) |> engine_stan(iter = 1000)
} # }
```
