# Specify a spatial explicit offset as bias

Including offsets is another option to integrate spatial prior
information in linear and additive regression models. Offsets shift the
intercept of the regression fit by a certain amount. Although only one
offset can be added to a regression model, it is possible to combine
several spatial-explicit estimates into one offset by calculating the
sum of all spatial-explicit layers.

## Usage

``` r
add_offset_bias(x, layer, add = TRUE, points = NULL)

# S4 method for class 'BiodiversityDistribution,SpatRaster'
add_offset_bias(x, layer, add = TRUE, points = NULL)
```

## Arguments

- x:

  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

- layer:

  A [`sf::sf`](https://r-spatial.github.io/sf/reference/sf.html) or
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object with the range for the target feature.

- add:

  [`logical`](https://rdrr.io/r/base/logical.html) specifying whether
  new offset is to be added. Setting this parameter to `FALSE` replaces
  the current offsets with the new one (Default: `TRUE`).

- points:

  An optional
  [`sf::sf`](https://r-spatial.github.io/sf/reference/sf.html) object
  with key points. The location of the points are then used to calculate
  the probability that a cell has been sampled while accounting for area
  differences. (Default: `NULL`).

## Value

Adds a bias offset to a
[`distribution`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
object.

## Details

This functions emulates the use of the
[`add_offset()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset.md)
function, however applies an inverse transformation to remove the
provided layer from the overall offset. So if for instance a offset is
already specified (such as area), this function removes the provided
`bias.layer` from it via `"offset(log(off.area)-log(bias.layer))"`

**Note that any transformation of the offset (such as `log`) has do be
done externally!**

If a generic offset is added, consider using the
[`add_offset()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset.md)
function. If the layer is a expert-based range and requires additional
parametrization, consider using the function
[`add_offset_range()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset_range.md)
or the `bossMaps` R-package.

## References

- Merow, C., Allen, J.M., Aiello-Lammens, M., Silander, J.A., 2016.
  Improving niche and range estimates with Maxent and point process
  models by integrating spatially explicit information. Glob. Ecol.
  Biogeogr. 25, 1022–1036. https://doi.org/10.1111/geb.12453

## See also

Other offset:
[`add_offset()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset.md),
[`add_offset_elevation()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset_elevation.md),
[`add_offset_range()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset_range.md),
[`rm_offset()`](https://iiasa.github.io/ibis.iSDM/reference/rm_offset.md)

## Examples

``` r
if (FALSE) { # \dontrun{
 x <- distribution(background) |>
   add_predictors(covariates) |>
   add_offset_bias(samplingBias)
} # }
```
