# Engine for Bayesian regularized regression models

Efficient MCMC algorithm for linear regression models that makes use of
'spike-and-slab' priors for some modest regularization on the amount of
posterior probability for a subset of the coefficients.

## Usage

``` r
engine_breg(
  x,
  iter = 10000,
  nthread = getOption("ibis.nthread"),
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

- iter:

  [`numeric`](https://rdrr.io/r/base/numeric.html) on the number of MCMC
  iterations to run (Default: `10000`).

- nthread:

  [`numeric`](https://rdrr.io/r/base/numeric.html) on the number of
  CPU-threads to use for data augmentation.

- type:

  The mode used for creating posterior predictions. Either making
  `"link"` or `"response"` (Default: `"response"`).

- ...:

  Other none specified parameters passed on to the model.

## Value

An
[Engine](https://iiasa.github.io/ibis.iSDM/reference/Engine-class.md).

## Details

This engine provides efficient Bayesian predictions through the Boom
R-package. However note that not all link and models functions are
supported and certain functionalities such as offsets are generally not
available. This engines allows the estimation of linear and non-linear
effects via the `"only_linear"` option specified in
[train](https://iiasa.github.io/ibis.iSDM/reference/train.md).

## References

- Nguyen, K., Le, T., Nguyen, V., Nguyen, T., & Phung, D. (2016,
  November). Multiple kernel learning with data augmentation. In Asian
  Conference on Machine Learning (pp. 49-64). PMLR.

- Steven L. Scott (2021). BoomSpikeSlab: MCMC for Spike and Slab
  Regression. R package version 1.2.4.
  https://CRAN.R-project.org/package=BoomSpikeSlab

## See also

Other engine:
[`engine_bart()`](https://iiasa.github.io/ibis.iSDM/reference/engine_bart.md),
[`engine_gdb()`](https://iiasa.github.io/ibis.iSDM/reference/engine_gdb.md),
[`engine_glm()`](https://iiasa.github.io/ibis.iSDM/reference/engine_glm.md),
[`engine_glmnet()`](https://iiasa.github.io/ibis.iSDM/reference/engine_glmnet.md),
[`engine_inla()`](https://iiasa.github.io/ibis.iSDM/reference/engine_inla.md),
[`engine_inlabru()`](https://iiasa.github.io/ibis.iSDM/reference/engine_inlabru.md),
[`engine_scampr()`](https://iiasa.github.io/ibis.iSDM/reference/engine_scampr.md),
[`engine_stan()`](https://iiasa.github.io/ibis.iSDM/reference/engine_stan.md),
[`engine_xgboost()`](https://iiasa.github.io/ibis.iSDM/reference/engine_xgboost.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Add BREG as an engine
x <- distribution(background) |> engine_breg(iter = 1000)
} # }
```
