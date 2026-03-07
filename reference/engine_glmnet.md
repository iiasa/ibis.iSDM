# Engine for regularized regression models

This engine allows the estimation of linear coefficients using either
ridge, lasso or elastic net regressions techniques. Backbone of this
engine is the glmnet R-package which is commonly used in SDMs, including
the popular `'maxnet'` (e.g. Maxent) package. Ultimately this engine is
an equivalent of
[engine_breg](https://iiasa.github.io/ibis.iSDM/reference/engine_breg.md),
but in a "frequentist" setting. If user aim to emulate a model that most
closely resembles maxent within the ibis.iSDM modelling framework, then
this package is the best way of doing so. Compared to the `'maxnet'`
R-package, a number of efficiency settings are implemented in particular
for cross-validation of alpha and lambda values.

Limited amount of prior information can be specified for this engine,
specifically via offsets or as
[`GLMNETPrior`](https://iiasa.github.io/ibis.iSDM/reference/GLMNETPrior.md),
which allow to specify priors as regularization constants.

## Usage

``` r
engine_glmnet(
  x,
  alpha = 0,
  nlambda = 100,
  lambda = NULL,
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

- alpha:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) giving the
  elasticnet mixing parameter, which has to be between `0` and `1`.
  `alpha=1` is the lasso penalty, and `alpha=0` the ridge penalty
  (Default: `0`).

- nlambda:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) giving the number
  of lambda values to be used (Default: `100`).

- lambda:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) with a user
  supplied estimate of lambda. Usually best to let this parameter be
  determined deterministically (Default: `NULL`).

- type:

  The mode used for creating posterior predictions. Either making
  `"link"` or `"response"` (Default: `"response"`).

- ...:

  Other parameters passed on to glmnet.

## Value

An
[Engine](https://iiasa.github.io/ibis.iSDM/reference/Engine-class.md).

## Details

Regularized regressions are effectively GLMs that are fitted with ridge,
lasso or elastic-net regularization. Which of them is chosen is
critically dependent on the alpha value:

- For `alpha` equal to `0` a ridge regularization is used. Ridge
  regularization has the property that it doesn't remove variables
  entirely, but instead sets their coefficients to `0`.

- For `alpha` equal to `1` a lasso regularization is used. Lassos tend
  to remove those coefficients fully from the final model that do not
  improve the loss function.

- For `alpha` values between `0` and `1` an elastic-net regularization
  is used, which is essentially a combination of the two. The optimal
  lambda parameter can be determined via cross-validation. For this
  option set `"varsel"` in
  [`train()`](https://iiasa.github.io/ibis.iSDM/reference/train.md) to
  `"reg"`.

## References

- Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
  Regularization Paths for Generalized Linear Models via Coordinate
  Descent. Journal of Statistical Software, 33(1), 1-22. URL
  https://www.jstatsoft.org/v33/i01/.

- Renner, I.W., Elith, J., Baddeley, A., Fithian, W., Hastie, T.,
  Phillips, S.J., Popovic, G. and Warton, D.I., 2015. Point process
  models for presence‐only analysis. Methods in Ecology and Evolution,
  6(4), pp.366-379.

- Fithian, W. & Hastie, T. (2013) Finite-sample equivalence in
  statistical models for presence-only data. The Annals of Applied
  Statistics 7, 1917–1939

## See also

Other engine:
[`engine_bart()`](https://iiasa.github.io/ibis.iSDM/reference/engine_bart.md),
[`engine_breg()`](https://iiasa.github.io/ibis.iSDM/reference/engine_breg.md),
[`engine_gdb()`](https://iiasa.github.io/ibis.iSDM/reference/engine_gdb.md),
[`engine_glm()`](https://iiasa.github.io/ibis.iSDM/reference/engine_glm.md),
[`engine_inla()`](https://iiasa.github.io/ibis.iSDM/reference/engine_inla.md),
[`engine_inlabru()`](https://iiasa.github.io/ibis.iSDM/reference/engine_inlabru.md),
[`engine_scampr()`](https://iiasa.github.io/ibis.iSDM/reference/engine_scampr.md),
[`engine_stan()`](https://iiasa.github.io/ibis.iSDM/reference/engine_stan.md),
[`engine_xgboost()`](https://iiasa.github.io/ibis.iSDM/reference/engine_xgboost.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Add GLMNET as an engine
x <- distribution(background) |> engine_glmnet(iter = 1000)
} # }
```
