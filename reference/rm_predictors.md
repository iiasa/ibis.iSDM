# Remove specific predictors from a [distribution](https://iiasa.github.io/ibis.iSDM/reference/distribution.md) object

Remove a particular variable from an
[distribution](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
object with a
[`PredictorDataset`](https://iiasa.github.io/ibis.iSDM/reference/PredictorDataset-class.md).
See Examples.

## Usage

``` r
rm_predictors(x, names)

# S4 method for class 'BiodiversityDistribution,character'
rm_predictors(x, names)
```

## Arguments

- x:

  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

- names:

  [`vector`](https://rdrr.io/r/base/vector.html) A Vector of character
  names describing the environmental stack.

## Examples

``` r
if (FALSE) { # \dontrun{
distribution(background) |>
 add_predictors(my_covariates) |>
 rm_predictors(names = "Urban")
} # }
```
