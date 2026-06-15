# Create priors from an existing distribution model

Often it can make sense to fit an additional model to get a grasp on the
range of values that "beta" parameters can take. This function takes an
existing
[`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md)
object and creates
[`PriorList`](https://iiasa.github.io/ibis.iSDM/reference/PriorList-class.md)
object from them. The resulting object can be used to add for instance
[priors](https://iiasa.github.io/ibis.iSDM/reference/priors.md) to a new
model.

## Usage

``` r
get_priors(mod, target_engine, ...)

# S4 method for class 'ANY,character'
get_priors(mod, target_engine, ...)
```

## Arguments

- mod:

  A fitted
  [`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
  object. If instead a
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md)
  object is passed to this function, it simply returns the contained
  priors used for estimation (if any).

- target_engine:

  A [`character`](https://rdrr.io/r/base/character.html) for which the
  priors should be created.

- ...:

  Other parameters passed down.

## Note

Not all engines support priors in similar ways. See the vignettes and
help pages on that topic!

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
[`XGBPriors()`](https://iiasa.github.io/ibis.iSDM/reference/XGBPriors.md),
[`add_priors()`](https://iiasa.github.io/ibis.iSDM/reference/add_priors.md),
[`priors()`](https://iiasa.github.io/ibis.iSDM/reference/priors.md),
[`rm_priors()`](https://iiasa.github.io/ibis.iSDM/reference/rm_priors.md)

## Examples

``` r
if (FALSE) { # \dontrun{
 mod <- distribution(background) |>
    add_predictors(covariates) |>
    add_biodiversity_poipo(points) |>
    engine_inlabru() |>
    train()
 get_priors(mod, target_engine = "BART")
} # }
```
