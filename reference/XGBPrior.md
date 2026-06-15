# Create a new monotonic prior for boosted regressions

Function to include prior information as monotonic constrain to an
extreme gradient descent boosting model
[`engine_xgboost`](https://iiasa.github.io/ibis.iSDM/reference/engine_xgboost.md).
Monotonic priors enforce directionality in direction of certain
variables, however specifying a monotonic constrain does not guarantee
that the variable is not regularized out during model fitting.

## Usage

``` r
XGBPrior(variable, hyper = "increasing", ...)

# S4 method for class 'character,character'
XGBPrior(variable, hyper = "increasing", ...)
```

## Arguments

- variable:

  A [`character`](https://rdrr.io/r/base/character.html) matched against
  existing predictors or latent effects.

- hyper:

  A [`character`](https://rdrr.io/r/base/character.html) object
  describing the type of constrain. Available options are
  `'increasing'`, `'decreasing'`, `'positive'`, `'negative'`, `'none'`.

- ...:

  Variables passed on to prior object.

## References

- Chen, T., He, T., Benesty, M., Khotilovich, V., Tang, Y., & Cho, H.
  (2015). Xgboost: extreme gradient boosting. R package version 0.4-2,
  1(4), 1-4.

## See also

[`Prior`](https://iiasa.github.io/ibis.iSDM/reference/Prior-class.md)
and
[`GDBPrior`](https://iiasa.github.io/ibis.iSDM/reference/GDBPrior.md).

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
[`XGBInteractionPrior()`](https://iiasa.github.io/ibis.iSDM/reference/XGBInteractionPrior.md),
[`XGBInteractionPriors()`](https://iiasa.github.io/ibis.iSDM/reference/XGBInteractionPriors.md),
[`XGBPriors()`](https://iiasa.github.io/ibis.iSDM/reference/XGBPriors.md),
[`add_priors()`](https://iiasa.github.io/ibis.iSDM/reference/add_priors.md),
[`get_priors()`](https://iiasa.github.io/ibis.iSDM/reference/get_priors.md),
[`priors()`](https://iiasa.github.io/ibis.iSDM/reference/priors.md),
[`rm_priors()`](https://iiasa.github.io/ibis.iSDM/reference/rm_priors.md)

## Examples

``` r
if (FALSE) { # \dontrun{
 pp <- XGBPrior("forest", "increasing")
} # }
```
