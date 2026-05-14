# Specify a spatial explicit offset

Including offsets is another option to integrate spatial prior
information in linear and additive regression models. Offsets shift the
intercept of the regression fit by a certain amount. Although only one
offset can be added to a regression model, it is possible to combine
several spatial-explicit estimates into one offset by calculating the
sum of all spatial-explicit layers.

## Usage

``` r
add_offset(x, layer, add = TRUE)

# S4 method for class 'BiodiversityDistribution,SpatRaster'
add_offset(x, layer, add = TRUE)

# S4 method for class 'BiodiversityDistribution,sf'
add_offset(x, layer, add = TRUE)
```

## Arguments

- x:

  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

- layer:

  A [`sf`](https://r-spatial.github.io/sf/reference/sf.html) or
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object with the range for the target feature.

- add:

  [`logical`](https://rdrr.io/r/base/logical.html) specifying whether
  new offset is to be added. Setting this parameter to `FALSE` replaces
  the current offsets with the new one (Default: `TRUE`).

## Value

Adds an offset to a
[`distribution`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
object.

## Details

This function allows to set any specific offset to a regression model.
The offset has to be provided as spatial
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
object. This function simply adds the layer to a
[`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
object. **Note that any transformation of the offset (such as `log`) has
do be done externally!**

If the layer is range and requires additional formatting, consider using
the function
[`add_offset_range()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset_range.md)
which has additional functionalities such such distance transformations.

## Note

Since offsets only make sense for linear regressions (and not for
instance regression tree based methods such as
[`engine_bart()`](https://iiasa.github.io/ibis.iSDM/reference/engine_bart.md)),
they do not work for all engines. Offsets specified for non-supported
engines are ignored during the estimation

## References

- Merow, C., Allen, J.M., Aiello-Lammens, M., Silander, J.A., 2016.
  Improving niche and range estimates with Maxent and point process
  models by integrating spatially explicit information. Glob. Ecol.
  Biogeogr. 25, 1022–1036. https://doi.org/10.1111/geb.12453

## See also

Other offset:
[`add_offset_bias()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset_bias.md),
[`add_offset_elevation()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset_elevation.md),
[`add_offset_range()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset_range.md),
[`rm_offset()`](https://iiasa.github.io/ibis.iSDM/reference/rm_offset.md)

## Examples

``` r
if (FALSE) { # \dontrun{
 x <- distribution(background) |>
   add_predictors(covariates) |>
   add_offset(nicheEstimate)
} # }
```
