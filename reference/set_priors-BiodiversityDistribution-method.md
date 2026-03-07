# Add priors to an existing distribution object

This function simply allows to add priors to an existing
[distribution](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
object. The supplied priors must be a
[`PriorList`](https://iiasa.github.io/ibis.iSDM/reference/PriorList-class.md)
object created through calling
[priors](https://iiasa.github.io/ibis.iSDM/reference/priors.md).

## Usage

``` r
# S4 method for class 'BiodiversityDistribution'
set_priors(x, priors = NULL, ...)
```

## Arguments

- x:

  [distribution](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

- priors:

  A
  [`PriorList`](https://iiasa.github.io/ibis.iSDM/reference/PriorList-class.md)
  object containing multiple priors.

- ...:

  Other parameters passed down.

## Note

Alternatively priors to environmental predictors can also directly added
as parameter via
[add_predictors](https://iiasa.github.io/ibis.iSDM/reference/add_predictors.md)

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
[`XGBPrior()`](https://iiasa.github.io/ibis.iSDM/reference/XGBPrior.md),
[`XGBPriors()`](https://iiasa.github.io/ibis.iSDM/reference/XGBPriors.md),
[`get_priors()`](https://iiasa.github.io/ibis.iSDM/reference/get_priors.md),
[`priors()`](https://iiasa.github.io/ibis.iSDM/reference/priors.md),
[`rm_priors()`](https://iiasa.github.io/ibis.iSDM/reference/rm_priors.md)

## Examples

``` r
if (FALSE) { # \dontrun{
 pp <-  GLMNETPrior("forest")
 x <- distribution(background) |>
  add_priors(pp)

} # }
```
