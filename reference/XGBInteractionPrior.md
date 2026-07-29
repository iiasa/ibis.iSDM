# Create a new interaction prior for XGBoost

Function to include prior information as interaction constraints in an
extreme gradient descent boosting model
[`engine_xgboost`](https://iiasa.github.io/ibis.iSDM/reference/engine_xgboost.md).
Interaction priors define groups of variables that are allowed to
interact in the same tree path. Variables outside the same group are not
allowed to interact.

## Usage

``` r
XGBInteractionPrior(variables, ...)

# S4 method for class 'character'
XGBInteractionPrior(variables, ...)
```

## Arguments

- variables:

  A [`character`](https://rdrr.io/r/base/character.html) vector matched
  against existing predictors or latent effects after XGBoost
  preprocessing.

- ...:

  Variables passed on to prior object.

## Value

A [`Prior`](https://iiasa.github.io/ibis.iSDM/reference/Prior-class.md)
object.

## Details

XGBoost interaction constraints are only supported by tree boosters.
They can be combined with monotonic constraints supplied through
[`XGBPrior`](https://iiasa.github.io/ibis.iSDM/reference/XGBPrior.md).

## See also

[`Prior`](https://iiasa.github.io/ibis.iSDM/reference/Prior-class.md),
[`XGBPrior`](https://iiasa.github.io/ibis.iSDM/reference/XGBPrior.md)
and
[`engine_xgboost`](https://iiasa.github.io/ibis.iSDM/reference/engine_xgboost.md).

Other prior:
[`BARTPrior()`](https://iiasa.github.io/ibis.iSDM/reference/BARTPrior.md),
[`BARTPriors()`](https://iiasa.github.io/ibis.iSDM/reference/BARTPriors.md),
[`BREGPrior()`](https://iiasa.github.io/ibis.iSDM/reference/BREGPrior.md),
[`BREGPriors()`](https://iiasa.github.io/ibis.iSDM/reference/BREGPriors.md),
[`GDBPrior()`](https://iiasa.github.io/ibis.iSDM/reference/GDBPrior.md),
[`GDBPriors()`](https://iiasa.github.io/ibis.iSDM/reference/GDBPriors.md),
[`GLMNETPrior()`](https://iiasa.github.io/ibis.iSDM/reference/GLMNETPrior.md),
[`GLMNETPriors()`](https://iiasa.github.io/ibis.iSDM/reference/GLMNETPriors.md),
[`INLAPrior()`](https://iiasa.github.io/ibis.iSDM/reference/INLAPrior.md),
[`INLAPriors()`](https://iiasa.github.io/ibis.iSDM/reference/INLAPriors.md),
[`STANPrior()`](https://iiasa.github.io/ibis.iSDM/reference/STANPrior.md),
[`STANPriors()`](https://iiasa.github.io/ibis.iSDM/reference/STANPriors.md),
[`XGBInteractionPriors()`](https://iiasa.github.io/ibis.iSDM/reference/XGBInteractionPriors.md),
[`XGBPrior()`](https://iiasa.github.io/ibis.iSDM/reference/XGBPrior.md),
[`XGBPriors()`](https://iiasa.github.io/ibis.iSDM/reference/XGBPriors.md),
[`add_priors()`](https://iiasa.github.io/ibis.iSDM/reference/add_priors.md),
[`get_priors()`](https://iiasa.github.io/ibis.iSDM/reference/get_priors.md),
[`priors()`](https://iiasa.github.io/ibis.iSDM/reference/priors.md),
[`rm_priors()`](https://iiasa.github.io/ibis.iSDM/reference/rm_priors.md)

## Examples

``` r
 pp <- XGBInteractionPrior(c("forest", "temperature"))
```
