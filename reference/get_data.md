# Small helper function to obtain predictions from an object

This function is a short helper function to return the fitted data from
a `DistributionModel` or `BiodiversityScenario` object. It can be used
to easily obtain for example the estimated prediction from a model or
the projected scenario from a
[`scenario()`](https://iiasa.github.io/ibis.iSDM/reference/scenario.md)
object.

## Usage

``` r
get_data(obj, what = NULL)

# S4 method for class 'ANY'
get_data(obj, what = NULL)
```

## Arguments

- obj:

  Provided
  [`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
  or
  [`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
  object.

- what:

  A [`character`](https://rdrr.io/r/base/character.html) of specific
  layer to be returned if existing (Default: `NULL`).

## Value

A `SpatRaster` or "stars" object depending on the input.

## Note

This function is essentially identical to querying the internal function
`x$get_data()` from the object. However it does attempt some lazy
character matching if what is supplied.

## Examples

``` r
if (FALSE) { # \dontrun{
 # Assumes previously computed model
 get_data(fit)
} # }
```
