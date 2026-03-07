# Create a new spike and slab prior for Bayesian generalized linear models

Function to include prior information via Zellner-style spike and slab
prior for generalized linear models used in
[engine_breg](https://iiasa.github.io/ibis.iSDM/reference/engine_breg.md).
These priors are similar to the horseshoe priors used in regularized
[engine_stan](https://iiasa.github.io/ibis.iSDM/reference/engine_stan.md)
models and penalize regressions by assuming most predictors having an
effect of `0`.

## Usage

``` r
BREGPrior(variable, hyper = NULL, ip = NULL)

# S4 method for class 'character'
BREGPrior(variable, hyper = NULL, ip = NULL)
```

## Arguments

- variable:

  A [`character`](https://rdrr.io/r/base/character.html) matched against
  existing predictors.

- hyper:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) estimate of the
  mean regression coefficients.

- ip:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) estimate between 0
  and 1 of the inclusion probability of the target variable (Default:
  `NULL`).

## Details

The Zellner-style spike and slab prior for generalized linear models are
specified as described in the Boom R-package. Currently supported are
two options which work for models with `Poisson` and `binomial`
(`Bernoulli`) distributed errors. Two types of priors can be provided on
a variable:

- `"coefficient"` Allows to specify Gaussian priors on the mean
  coefficients of the model. Priors on the coefficients can be provided
  via the `"hyper"` parameter. Note that variables with such a prior can
  still be regularized out from the model.

- `"inclusion.probability"` A
  [`vector`](https://rdrr.io/r/base/vector.html) giving the prior
  probability of inclusion for the specified variable. This can be
  useful when prior information on preference is known but not the
  strength of it.

If coefficients are set, then the inclusion probability is also modified
by default. However even when not knowing a particular estimate of a
beta coefficients and their direction, one can still provide an estimate
of the inclusion probability. In other words: **The hyperparameters
'hyper' and 'ip' can't be both `NULL`.**

## References

- Hugh Chipman, Edward I. George, Robert E. McCulloch, M. Clyde, Dean P.
  Foster, Robert A. Stine (2001), "The Practical Implementation of
  Bayesian Model Selection" Lecture Notes-Monograph Series, Vol. 38, pp.
  65-134. Institute of Mathematical Statistics.

## See also

[`Prior`](https://iiasa.github.io/ibis.iSDM/reference/Prior-class.md)

Other prior:
[`BARTPrior()`](https://iiasa.github.io/ibis.iSDM/reference/BARTPrior.md),
[`BARTPriors()`](https://iiasa.github.io/ibis.iSDM/reference/BARTPriors.md),
[`BREGPriors()`](https://iiasa.github.io/ibis.iSDM/reference/BREGPriors.md),
[`GDBPrior()`](https://iiasa.github.io/ibis.iSDM/reference/GDBPrior.md),
[`GDBPriors()`](https://iiasa.github.io/ibis.iSDM/reference/GDBPriors.md),
[`GLMNETPrior()`](https://iiasa.github.io/ibis.iSDM/reference/GLMNETPrior.md),
[`GLMNETPriors()`](https://iiasa.github.io/ibis.iSDM/reference/GLMNETPriors.md),
[`INLAPrior()`](https://iiasa.github.io/ibis.iSDM/reference/INLAPrior.md),
[`INLAPriors()`](https://iiasa.github.io/ibis.iSDM/reference/INLAPriors.md),
[`STANPrior()`](https://iiasa.github.io/ibis.iSDM/reference/STANPrior.md),
[`STANPriors()`](https://iiasa.github.io/ibis.iSDM/reference/STANPriors.md),
[`XGBPrior()`](https://iiasa.github.io/ibis.iSDM/reference/XGBPrior.md),
[`XGBPriors()`](https://iiasa.github.io/ibis.iSDM/reference/XGBPriors.md),
[`add_priors()`](https://iiasa.github.io/ibis.iSDM/reference/add_priors.md),
[`get_priors()`](https://iiasa.github.io/ibis.iSDM/reference/get_priors.md),
[`priors()`](https://iiasa.github.io/ibis.iSDM/reference/priors.md),
[`rm_priors()`](https://iiasa.github.io/ibis.iSDM/reference/rm_priors.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Positive coefficient
p1 <- BREGPrior(variable = "forest", hyper = 2, ip = NULL)
p1
# Coefficient and direction unknown but variable def. important
p2 <- BREGPrior(variable = "forest", hyper = NULL, ip = 1)
p2
} # }
```
