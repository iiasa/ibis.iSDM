# Select specific predictors from a [distribution](https://iiasa.github.io/ibis.iSDM/reference/distribution.md) object

This function allows - out of a
[`character`](https://rdrr.io/r/base/character.html) vector with the
names of an already added
[`PredictorDataset`](https://iiasa.github.io/ibis.iSDM/reference/PredictorDataset-class.md)
object - to select a particular set of predictors. See Examples.

## Usage

``` r
sel_predictors(x, names)

# S4 method for class 'BiodiversityDistribution,character'
sel_predictors(x, names)
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
 sel_predictors(names = c("Forest", "Elevation"))
} # }
```
