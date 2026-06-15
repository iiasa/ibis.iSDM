# Create a new STAN prior

Function to create a new prior for
[engine_stan](https://iiasa.github.io/ibis.iSDM/reference/engine_stan.md)
models. Priors currently can be set on specific environmental
predictors.

## Usage

``` r
STANPrior(variable, type, hyper = c(0, 2), ...)

# S4 method for class 'character,character'
STANPrior(variable, type, hyper = c(0, 2), ...)
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
  hyper parameters. First entry is treated as mean (Default: `0`), the
  second as the standard variation (Default: `2`) of a Gaussian
  distribution on the respective coefficient.

- ...:

  Variables passed on to prior object.

## References

- Lemoine, N. P. (2019). Moving beyond noninformative priors: why and
  how to choose weakly informative priors in Bayesian analyses. Oikos,
  128(7), 912-928.

- Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B.,
  Betancourt, M., ... & Riddell, A. (2017). Stan: A probabilistic
  programming language. Journal of statistical software, 76(1), 1-32.

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
[`INLAPrior()`](https://iiasa.github.io/ibis.iSDM/reference/INLAPrior.md),
[`INLAPriors()`](https://iiasa.github.io/ibis.iSDM/reference/INLAPriors.md),
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
 pp <- STANPrior("forest", "normal", c(0,1))
} # }
```
