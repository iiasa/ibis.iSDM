# Helper function when multiple variables and types are supplied for INLA priors

This is a helper function to specify several
[INLAPrior](https://iiasa.github.io/ibis.iSDM/reference/INLAPrior.md)
objects with the same hyper-parameters, but different variables.

## Usage

``` r
INLAPriors(variables, type, hyper = c(0, 0.001), ...)

# S4 method for class 'vector,character'
INLAPriors(variables, type, hyper = c(0, 0.001), ...)
```

## Arguments

- variables:

  A [`vector`](https://rdrr.io/r/base/vector.html) of
  [`character`](https://rdrr.io/r/base/character.html) matched against
  existing predictors or latent effects.

- type:

  A [`character`](https://rdrr.io/r/base/character.html) specifying the
  type of prior to be set.

- hyper:

  A [`vector`](https://rdrr.io/r/base/vector.html) with
  [`numeric`](https://rdrr.io/r/base/numeric.html) values to be used as
  hyper-parameters.

- ...:

  Variables passed on to prior object.

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
[`STANPrior()`](https://iiasa.github.io/ibis.iSDM/reference/STANPrior.md),
[`STANPriors()`](https://iiasa.github.io/ibis.iSDM/reference/STANPriors.md),
[`XGBInteractionPrior()`](https://iiasa.github.io/ibis.iSDM/reference/XGBInteractionPrior.md),
[`XGBInteractionPriors()`](https://iiasa.github.io/ibis.iSDM/reference/XGBInteractionPriors.md),
[`XGBPrior()`](https://iiasa.github.io/ibis.iSDM/reference/XGBPrior.md),
[`XGBPriors()`](https://iiasa.github.io/ibis.iSDM/reference/XGBPriors.md),
[`add_priors()`](https://iiasa.github.io/ibis.iSDM/reference/add_priors.md),
[`get_priors()`](https://iiasa.github.io/ibis.iSDM/reference/get_priors.md),
[`priors()`](https://iiasa.github.io/ibis.iSDM/reference/priors.md),
[`rm_priors()`](https://iiasa.github.io/ibis.iSDM/reference/rm_priors.md)
