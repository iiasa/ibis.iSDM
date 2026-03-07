# Engine for extreme gradient boosting (XGBoost)

Allows to estimate eXtreme gradient descent boosting for tree-based or
linear boosting regressions. The XGBoost engine is a flexible, yet
powerful engine with many customization options, supporting multiple
options to perform single and multi-class regression and classification
tasks. For a full list of options users are advised to have a look at
the [xgboost::xgb.train](https://rdrr.io/pkg/xgboost/man/xgb.train.html)
help file and <https://xgboost.readthedocs.io>.

## Usage

``` r
engine_xgboost(
  x,
  booster = "gbtree",
  iter = 8000L,
  learning_rate = 0.001,
  gamma = 6,
  reg_lambda = 0,
  reg_alpha = 0,
  max_depth = 2,
  subsample = 0.75,
  colsample_bytree = 0.4,
  min_child_weight = 3,
  nthread = getOption("ibis.nthread"),
  ...
)
```

## Arguments

- x:

  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

- booster:

  A [`character`](https://rdrr.io/r/base/character.html) of the booster
  to use. Either `"gbtree"` or `"gblinear"` (Default: `gblinear`)

- iter:

  [`numeric`](https://rdrr.io/r/base/numeric.html) value giving the the
  maximum number of boosting iterations for cross-validation (Default:
  `8e3L`).

- learning_rate:

  [`numeric`](https://rdrr.io/r/base/numeric.html) value indicating the
  learning rate (eta). Lower values generally being better but also
  computationally more costly. (Default: `1e-3`)

- gamma:

  [`numeric`](https://rdrr.io/r/base/numeric.html) A regularization
  parameter in the model. Lower values for better estimates (Default:
  `3`). Also see `"reg_lambda"` parameter for the L2 regularization on
  the weights

- reg_lambda:

  [`numeric`](https://rdrr.io/r/base/numeric.html) L2 regularization
  term on weights (Default: `0`).

- reg_alpha:

  [`numeric`](https://rdrr.io/r/base/numeric.html) L1 regularization
  term on weights (Default: `0`).

- max_depth:

  [`numeric`](https://rdrr.io/r/base/numeric.html) The Maximum depth of
  a tree (Default: `3`).

- subsample:

  [`numeric`](https://rdrr.io/r/base/numeric.html) The ratio used for
  subsampling to prevent overfitting. Also used for creating a random
  tresting dataset (Default: `0.75`).

- colsample_bytree:

  [`numeric`](https://rdrr.io/r/base/numeric.html) Sub-sample ratio of
  columns when constructing each tree (Default: `0.4`).

- min_child_weight:

  [`numeric`](https://rdrr.io/r/base/numeric.html) Broadly related to
  the number of instances necessary for each node (Default: `3`).

- nthread:

  [`numeric`](https://rdrr.io/r/base/numeric.html) on the number of
  CPU-threads to use.

- ...:

  Other none specified parameters.

## Value

An
[Engine](https://iiasa.github.io/ibis.iSDM/reference/Engine-class.md).

## Details

The default parameters have been set relatively conservative as to
reduce overfitting.

XGBoost supports the specification of monotonic constraints on certain
variables. Within ibis this is possible via
[`XGBPrior`](https://iiasa.github.io/ibis.iSDM/reference/XGBPrior.md).
However constraints are available only for the `"gbtree"` baselearners.

## Note

*'Machine learning is statistics minus any checking of models and
assumptions‘* ~ Brian D. Ripley, useR! 2004, Vienna

## References

- Tianqi Chen and Carlos Guestrin, "XGBoost: A Scalable Tree Boosting
  System", 22nd SIGKDD Conference on Knowledge Discovery and Data
  Mining, 2016, https://arxiv.org/abs/1603.02754

## See also

[xgboost::xgb.train](https://rdrr.io/pkg/xgboost/man/xgb.train.html)

Other engine:
[`engine_bart()`](https://iiasa.github.io/ibis.iSDM/reference/engine_bart.md),
[`engine_breg()`](https://iiasa.github.io/ibis.iSDM/reference/engine_breg.md),
[`engine_gdb()`](https://iiasa.github.io/ibis.iSDM/reference/engine_gdb.md),
[`engine_glm()`](https://iiasa.github.io/ibis.iSDM/reference/engine_glm.md),
[`engine_glmnet()`](https://iiasa.github.io/ibis.iSDM/reference/engine_glmnet.md),
[`engine_inla()`](https://iiasa.github.io/ibis.iSDM/reference/engine_inla.md),
[`engine_inlabru()`](https://iiasa.github.io/ibis.iSDM/reference/engine_inlabru.md),
[`engine_scampr()`](https://iiasa.github.io/ibis.iSDM/reference/engine_scampr.md),
[`engine_stan()`](https://iiasa.github.io/ibis.iSDM/reference/engine_stan.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Add xgboost as an engine
x <- distribution(background) |> engine_xgboost(iter = 4000)
} # }
```
