# Approximate missing time steps between dates

This function linearly approximates shares between time steps, so that
gaps for instance between 2010 and 2020 are filled with data for 2010,
2011, 2012, etc.

## Usage

``` r
interpolate_gaps(env, date_interpolation = "annual", method = "linear")
```

## Arguments

- env:

  A [`stars`](https://rdrr.io/r/graphics/stars.html) object.

- date_interpolation:

  [`character`](https://rdrr.io/r/base/character.html) on how missing
  dates between events should be interpolated. See
  [`project()`](https://iiasa.github.io/ibis.iSDM/reference/project.md).

- method:

  A [`character`](https://rdrr.io/r/base/character.html) on the used
  method for approximation, either `"linear"` (Default) or `"constant"`
  through a step function.

## Value

[`logical`](https://rdrr.io/r/base/logical.html) indicating if the two
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
objects have the same

## Examples

``` r
if (FALSE) { # \dontrun{
  # Interpolate stars stack
  sc <- interpolate_gaps( stack, "annual")
} # }
```
