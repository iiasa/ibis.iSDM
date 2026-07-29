# Helper function when multiple variables are supplied for GLMNET priors

This is a helper function to specify several
[GLMNETPrior](https://iiasa.github.io/ibis.iSDM/reference/GLMNETPrior.md)
with the same hyper-parameters, but different variables.

## Usage

``` r
GLMNETPriors(variable, hyper = 0, lims = c(-Inf, Inf))

# S4 method for class 'character'
GLMNETPriors(variable, hyper = 0, lims = c(-Inf, Inf))
```

## Arguments

- variable:

  A [`character`](https://rdrr.io/r/base/character.html) variable passed
  on to the prior object.

- hyper:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) value between `0`
  and `1` that state the penalization factor. By default this is set to
  `0`, implying the `"variable"` provided is not regularized at all.

- lims:

  A [`numeric`](https://rdrr.io/r/base/numeric.html)
  [`vector`](https://rdrr.io/r/base/vector.html) of the lower and upper
  limits for each coefficient (Default: `c(-Inf, Inf)`).

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
[`INLAPrior()`](https://iiasa.github.io/ibis.iSDM/reference/INLAPrior.md),
[`INLAPriors()`](https://iiasa.github.io/ibis.iSDM/reference/INLAPriors.md),
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

## Examples

``` r
priors <- GLMNETPriors(c("forest", "temperature"), hyper = 0)
length(priors)
#> [1] 2
```
