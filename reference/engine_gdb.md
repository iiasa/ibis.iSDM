# Use of Gradient Descent Boosting for model estimation

Gradient descent boosting is an efficient way to optimize any loss
function of a generalized linear or additive model (such as the GAMs
available through the `"mgcv"` R-package). It furthermore automatically
regularizes the fit, thus the resulting model only contains the
covariates whose baselearners have some influence on the response.
Depending on the type of the `add_biodiversity` data, either poisson
process models or logistic regressions are estimated. If the
`"only_linear"` term in
[train](https://iiasa.github.io/ibis.iSDM/reference/train.md) is set to
`FALSE`, splines are added to the estimation, thus providing a
non-linear additive inference.

## Usage

``` r
engine_gdb(
  x,
  iter = 2000,
  learning_rate = 0.1,
  empirical_risk = "inbag",
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

  An [`integer`](https://rdrr.io/r/base/integer.html) giving the number
  of boosting iterations (Default: `2e3L`).

- learning_rate:

  A bounded [`numeric`](https://rdrr.io/r/base/numeric.html) value
  between `0` and `1` defining the shrinkage parameter.

- empirical_risk:

  method for empirical risk calculation. Available options are
  `'inbag'`, `'oobag'` and `'none'`. (Default: `'inbag'`).

- type:

  The mode used for creating posterior predictions. Either making
  `"link"`, `"response"` or `"class"` (Default: `"response"`).

- ...:

  Other variables or control parameters

## Value

An engine.

## Details

This engine requires the `"mboost"` R-package to be installed. It is in
philosophy somewhat related to the
[engine_xgboost](https://iiasa.github.io/ibis.iSDM/reference/engine_xgboost.md)
and `"XGBoost"` R-package, however providing some additional desirable
features that make estimation quicker and particularly useful for
spatial projections. Such as for instance the ability to specifically
add spatial baselearners via
[add_latent_spatial](https://iiasa.github.io/ibis.iSDM/reference/add_latent_spatial.md)
or the specification of monotonically constrained priors via
[GDBPrior](https://iiasa.github.io/ibis.iSDM/reference/GDBPrior.md).

## Note

The coefficients resulting from gdb with poipa data (Binomial) are only
0.5 of the typical coefficients of a logit model obtained via glm (see
Binomial).

## References

- Hofner, B., Mayr, A., Robinzonov, N., & Schmid, M. (2014). Model-based
  boosting in R: a hands-on tutorial using the R package mboost.
  Computational statistics, 29(1-2), 3-35.

- Hofner, B., Müller, J., Hothorn, T., (2011). Monotonicity-constrained
  species distribution models. Ecology 92, 1895–901.

- Mayr, A., Hofner, B. and Schmid, M. (2012). The importance of knowing
  when to stop - a sequential stopping rule for component-wise gradient
  boosting. Methods of Information in Medicine, 51, 178–186.

## See also

Other engine:
[`engine_bart()`](https://iiasa.github.io/ibis.iSDM/reference/engine_bart.md),
[`engine_breg()`](https://iiasa.github.io/ibis.iSDM/reference/engine_breg.md),
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
# Add GDB as an engine
x <- distribution(background) |> engine_gdb(iter = 1000)
} # }
```
