# Create a new INLA prior

For any fixed and random effect INLA supports a range of different
priors of exponential distributions.

Currently supported for INLA in ibis.iSDM are the following priors that
can be specified via `"type"`:

- `"normal"` or `"gaussian"`: Priors on normal distributed and set to
  specified variable. Required parameters are a mean and a precision
  estimate provided to `"hyper"`. Note that precision is not equivalent
  (rather the inverse) to typical standard deviation specified in
  Gaussian priors. Defaults are set to a mean of `0` and a precision of
  `0.001`.

- `"clinear"`: Prior that places a constraint on the linear coefficients
  of a model so as that the coefficient is in a specified interval
  `"c(lower,upper)"`. Specified through hyper these values can be
  negative, positive or infinite.

- `"spde"`, specifically `'prior.range'` and `'prior.sigma'`:
  Specification of penalized complexity priors which can be added to a
  SPDE spatial random effect added via
  [`add_latent_spatial()`](https://iiasa.github.io/ibis.iSDM/reference/add_latent_spatial.md).
  Here the range of the penalized complexity prior can be specified
  through `'prior.range'` and the uncertainty via `'prior.sigma'` both
  supplied to the options 'type' and 'hyper'.

Other priors available in INLA ` names(INLA::inla.models()$prior) ) `
might also work, but have not been tested!

## Usage

``` r
INLAPrior(variable, type = "normal", hyper = c(0, 0.001), ...)

# S4 method for class 'character,character'
INLAPrior(variable, type = "normal", hyper = c(0, 0.001), ...)
```

## Arguments

- variable:

  A [`character`](https://rdrr.io/r/base/character.html) matched against
  existing predictors or latent effects.

- type:

  A [`character`](https://rdrr.io/r/base/character.html) specifying the
  type of prior to be set.

- hyper:

  A [`vector`](https://rdrr.io/r/base/vector.html) with
  [`numeric`](https://rdrr.io/r/base/numeric.html) values to be used as
  hyper-parameters. See description. The default values are set to a
  mean of `0` and a precision of `0.001`.

- ...:

  Variables passed on to prior object.

## Note

Compared to other engines, INLA does unfortunately does not support
priors related to more stringent parameter regularization such as
Laplace or Horseshoe priors, which limits the capability of
[`engine_inla`](https://iiasa.github.io/ibis.iSDM/reference/engine_inla.md)
for regularization. That being said many of the default uninformative
priors act already regularize the coefficients to some degree.

## References

- Rue, H., Riebler, A., Sørbye, S. H., Illian, J. B., Simpson, D. P., &
  Lindgren, F. K. (2017). Bayesian computing with INLA: a review. Annual
  Review of Statistics and Its Application, 4, 395-421.

- Simpson, D., Rue, H., Riebler, A., Martins, T. G., & Sørbye, S. H.
  (2017). Penalising model component complexity: A principled, practical
  approach to constructing priors. Statistical science, 32(1), 1-28.

## See also

[`Prior`](https://iiasa.github.io/ibis.iSDM/reference/Prior-class.md).

Other prior:
[`BARTPrior()`](https://iiasa.github.io/ibis.iSDM/reference/BARTPrior.md),
[`BARTPriors()`](https://iiasa.github.io/ibis.iSDM/reference/BARTPriors.md),
[`BREGPrior()`](https://iiasa.github.io/ibis.iSDM/reference/BREGPrior.md),
[`BREGPriors()`](https://iiasa.github.io/ibis.iSDM/reference/BREGPriors.md),
[`GDBPrior()`](https://iiasa.github.io/ibis.iSDM/reference/GDBPrior.md),
[`GDBPriors()`](https://iiasa.github.io/ibis.iSDM/reference/GDBPriors.md),
[`GLMNETPrior()`](https://iiasa.github.io/ibis.iSDM/reference/GLMNETPrior.md),
[`GLMNETPriors()`](https://iiasa.github.io/ibis.iSDM/reference/GLMNETPriors.md),
[`INLAPriors()`](https://iiasa.github.io/ibis.iSDM/reference/INLAPriors.md),
[`STANPrior()`](https://iiasa.github.io/ibis.iSDM/reference/STANPrior.md),
[`STANPriors()`](https://iiasa.github.io/ibis.iSDM/reference/STANPriors.md),
[`XGBPrior()`](https://iiasa.github.io/ibis.iSDM/reference/XGBPrior.md),
[`XGBPriors()`](https://iiasa.github.io/ibis.iSDM/reference/XGBPriors.md),
[`add_priors()`](https://iiasa.github.io/ibis.iSDM/reference/add_priors.md),
[`get_priors()`](https://iiasa.github.io/ibis.iSDM/reference/get_priors.md),
[`priors()`](https://iiasa.github.io/ibis.iSDM/reference/priors.md),
[`rm_priors()`](https://iiasa.github.io/ibis.iSDM/reference/rm_priors.md)
