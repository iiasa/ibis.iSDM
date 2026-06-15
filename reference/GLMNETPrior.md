# Regression penalty priors for GLMNET

The
[`engine_glmnet`](https://iiasa.github.io/ibis.iSDM/reference/engine_glmnet.md)
engine does not support priors in a typical sense, however it is
possible to specify so called penalty factors as well as lower and upper
limits on all variables in the model.

The default penalty multiplier is `1` for each coefficient X covariate,
i.e. coefficients are penalized equally and then informed by an
intersection of any absence information with the covariates. In contrast
a variable with penalty.factor equal to `0` is not penalized at all.

In addition, it is possible to specify a lower and upper limit for
specific coefficients, which constrain them to a certain range. By
default those ranges are set to `-Inf` and `Inf` respectively, but can
be reset to a specific value range by altering `"lims"` (see examples).

For a regularized regression that supports a few more options on the
priors, check out the Bayesian
[`engine_breg`](https://iiasa.github.io/ibis.iSDM/reference/engine_breg.md).

## Usage

``` r
GLMNETPrior(variable, hyper = 0, lims = c(-Inf, Inf), ...)

# S4 method for class 'character'
GLMNETPrior(variable, hyper = 0, lims = c(-Inf, Inf), ...)
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

- ...:

  Variables passed on to prior object.

## See also

[`Prior`](https://iiasa.github.io/ibis.iSDM/reference/Prior-class.md)

Other prior:
[`BARTPrior()`](https://iiasa.github.io/ibis.iSDM/reference/BARTPrior.md),
[`BARTPriors()`](https://iiasa.github.io/ibis.iSDM/reference/BARTPriors.md),
[`BREGPrior()`](https://iiasa.github.io/ibis.iSDM/reference/BREGPrior.md),
[`BREGPriors()`](https://iiasa.github.io/ibis.iSDM/reference/BREGPriors.md),
[`GDBPrior()`](https://iiasa.github.io/ibis.iSDM/reference/GDBPrior.md),
[`GDBPriors()`](https://iiasa.github.io/ibis.iSDM/reference/GDBPriors.md),
[`GLMNETPriors()`](https://iiasa.github.io/ibis.iSDM/reference/GLMNETPriors.md),
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
if (FALSE) { # \dontrun{
# Retain variable
p1 <- GLMNETPrior(variable = "forest", hyper = 0)
p1
# Smaller chance to be regularized
p2 <- GLMNETPrior(variable = "forest", hyper = 0.2, lims = c(0, Inf))
p2
} # }
```
