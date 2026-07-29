# Create a tree-based split probability prior for BART

Function to include prior information as split probability for the
Bayesian additive regression tree model added via
[engine_bart](https://iiasa.github.io/ibis.iSDM/reference/engine_bart.md).

Priors for
[engine_bart](https://iiasa.github.io/ibis.iSDM/reference/engine_bart.md)
have to be specified as transition probabilities of variables which are
internally used to generate splits in the regression tree. Specifying a
prior can thus help to 'enforce' a split with a given variable. These
can be numeric and coded as values between `0` and `1`.

## Usage

``` r
BARTPrior(variable, hyper = 0.75, ...)

# S4 method for class 'character'
BARTPrior(variable, hyper = 0.75, ...)
```

## Arguments

- variable:

  A [`character`](https://rdrr.io/r/base/character.html) matched against
  existing predictors or latent effects.

- hyper:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) object with a value
  between `0` and `1`. Defaults to `0.75`.

- ...:

  Variables passed on to prior object.

## Value

A [`Prior`](https://iiasa.github.io/ibis.iSDM/reference/Prior-class.md)
object.

## Note

Even if a given variable is included as split in the regression or
classification tree, this does not necessarily mean that the prediction
changes if the value is non-informative (as the split can occur early
on). It does however affect any variable importance estimates calculated
from the model.

## References

- Chipman, H., George, E., and McCulloch, R. (2009) BART: Bayesian
  Additive Regression Trees.

- Chipman, H., George, E., and McCulloch R. (2006) Bayesian Ensemble
  Learning. Advances in Neural Information Processing Systems 19,
  Scholkopf, Platt and Hoffman, Eds., MIT Press, Cambridge, MA, 265-272.

## See also

[`Prior`](https://iiasa.github.io/ibis.iSDM/reference/Prior-class.md).

Other prior:
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
[`XGBPriors()`](https://iiasa.github.io/ibis.iSDM/reference/XGBPriors.md),
[`add_priors()`](https://iiasa.github.io/ibis.iSDM/reference/add_priors.md),
[`get_priors()`](https://iiasa.github.io/ibis.iSDM/reference/get_priors.md),
[`priors()`](https://iiasa.github.io/ibis.iSDM/reference/priors.md),
[`rm_priors()`](https://iiasa.github.io/ibis.iSDM/reference/rm_priors.md)

## Examples

``` r
prior <- BARTPrior("forest", hyper = 0.75)
prior$get("value")
#> [1] 0.75
```
