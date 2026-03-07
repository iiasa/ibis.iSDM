# Obtain spatial partial effects of trained model

Similar to
[partial](https://iiasa.github.io/ibis.iSDM/reference/partial.md), this
function calculates a partial response of a trained model for a given
variable. Differently from
[partial](https://iiasa.github.io/ibis.iSDM/reference/partial.md) in
space. However the result is a
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
showing the spatial magnitude of the partial response.

## Usage

``` r
spartial(mod, x.var, constant = NULL, newdata = NULL, plot = FALSE, ...)

# S4 method for class 'ANY,character'
spartial(mod, x.var, constant = NULL, newdata = NULL, plot = FALSE, ...)

spartial.DistributionModel(mod, ...)
```

## Arguments

- mod:

  A
  [`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
  object with trained model.

- x.var:

  A [character](https://rdrr.io/r/base/character.html) indicating the
  variable for which a partial effect is to be calculated.

- constant:

  A [numeric](https://rdrr.io/r/base/numeric.html) constant to be
  inserted for all other variables. Default calculates the
  [mean](https://rspatial.github.io/terra/reference/summarize-generics.html)
  per variable.

- newdata:

  A [`data.frame`](https://rdrr.io/r/base/data.frame.html) on which to
  calculate the spartial for. Can be for example created from a raster
  file (Default: `NULL`).

- plot:

  A [logical](https://rdrr.io/r/base/logical.html) indication of whether
  the result is to be plotted?

- ...:

  Other engine specific parameters.

## Value

A
[terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
containing the mapped partial response of the variable.

## Details

By default the
[mean](https://rspatial.github.io/terra/reference/summarize-generics.html)
is calculated across all parameters that are not `x.var`. Instead a
*constant* can be set (for instance `0`) to be applied to the output.

## See also

[partial](https://iiasa.github.io/ibis.iSDM/reference/partial.md)

## Examples

``` r
if (FALSE) { # \dontrun{
 # Create and visualize the spartial effect
 spartial(fit, x.var = "Forest.cover", plot = TRUE)
} # }
```
