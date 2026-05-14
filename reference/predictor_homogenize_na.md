# Homogenize NA values across a set of predictors.

This method allows the homogenization of missing data across a set of
environmental predictors. It is by default called when predictors are
added to
[`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md)
object. Only grid cells with NAs that contain values at some raster
layers are homogenized. Additional parameters allow instead of
homogenization to fill the missing data with neighbouring values

## Usage

``` r
predictor_homogenize_na(
  env,
  fill = FALSE,
  fill_method = "ngb",
  return_na_cells = FALSE
)
```

## Arguments

- env:

  A
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object with the predictors.

- fill:

  A [`logical`](https://rdrr.io/r/base/logical.html) value indicating
  whether missing data are to be filled (Default: `FALSE`).

- fill_method:

  A [`character`](https://rdrr.io/r/base/character.html) of the method
  for filling gaps to be used (Default: `'ngb'`).

- return_na_cells:

  A [`logical`](https://rdrr.io/r/base/logical.html) value of whether
  the ids of grid cells with NA values is to be returned instead
  (Default: `FALSE`).

## Value

A
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
object with the same number of layers as the input.

## Examples

``` r
if (FALSE) { # \dontrun{
 # Harmonize predictors
 env <- predictor_homogenize_na(env)
} # }
```
