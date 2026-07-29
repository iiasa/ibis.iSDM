# Helper function when multiple variables are supplied for XGBoost priors

This is a helper function to specify several
[XGBPrior](https://iiasa.github.io/ibis.iSDM/reference/XGBPrior.md) with
the same hyper-parameters, but different variables.

## Usage

``` r
XGBPriors(variable, hyper = "increasing", ...)

# S4 method for class 'character'
XGBPriors(variable, hyper = "increasing", ...)
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

## Value

A named [`list`](https://rdrr.io/r/base/list.html) of
[`Prior`](https://iiasa.github.io/ibis.iSDM/reference/Prior-class.md)
objects.

## See also

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
[`XGBPrior()`](https://iiasa.github.io/ibis.iSDM/reference/XGBPrior.md),
[`add_priors()`](https://iiasa.github.io/ibis.iSDM/reference/add_priors.md),
[`get_priors()`](https://iiasa.github.io/ibis.iSDM/reference/get_priors.md),
[`priors()`](https://iiasa.github.io/ibis.iSDM/reference/priors.md),
[`rm_priors()`](https://iiasa.github.io/ibis.iSDM/reference/rm_priors.md)

## Examples

``` r
priors <- XGBPriors(c("forest", "temperature"), hyper = "increasing")
length(priors)
#> [1] 2
```
