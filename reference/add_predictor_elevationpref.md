# Create lower and upper limits for an elevational range and add them as separate predictors

Create lower and upper limits for an elevational range and add them as
separate predictors

## Usage

``` r
add_predictor_elevationpref(x, layer, lower, upper, transform = "none")

# S4 method for class 'BiodiversityDistribution,ANY,numeric,numeric'
add_predictor_elevationpref(x, layer, lower, upper, transform = "none")
```

## Arguments

- x:

  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

- layer:

  A [`character`](https://rdrr.io/r/base/character.html) stating the
  elevational layer in the Distribution object or
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object.

- lower:

  [`numeric`](https://rdrr.io/r/base/numeric.html) value for a lower
  elevational preference of a species.

- upper:

  [`numeric`](https://rdrr.io/r/base/numeric.html) value for a upper
  elevational preference of a species.

- transform:

  [`character`](https://rdrr.io/r/base/character.html) Any optional
  transformation to be applied. Usually not needed (Default: `"none"`).

## Examples

``` r
if (FALSE) { # \dontrun{
distribution(background) |>
  add_predictor_elevationpref(elevation, lower = 200, upper = 1000)
} # }
```
