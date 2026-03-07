# Specify a expert-based species range as offset

This function has additional options compared to the more generic
[`add_offset()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset.md),
allowing customized options specifically for expert-based ranges as
offsets or spatialized polygon information on species occurrences. If
even more control is needed, the user is informed of the `"bossMaps"`
package Merow et al. (2017). Some functionalities of that package
emulated through the `"distance_function"` set to `"log"`. This tries to
fit a 5-parameter logistic function to estimate the distance from the
range (Merow et al. 2017).

## Usage

``` r
add_offset_range(
  x,
  layer,
  distance_max = Inf,
  family = "poisson",
  presence_prop = 0.9,
  distance_clip = FALSE,
  distance_function = "negexp",
  field_occurrence = "observed",
  fraction = NULL,
  point = FALSE,
  add = TRUE
)

# S4 method for class 'BiodiversityDistribution,SpatRaster'
add_offset_range(x, layer, fraction = NULL, add = TRUE)

# S4 method for class 'BiodiversityDistribution,sf'
add_offset_range(
  x,
  layer,
  distance_max = Inf,
  family = "poisson",
  presence_prop = 0.9,
  distance_clip = FALSE,
  distance_function = "negexp",
  field_occurrence = "observed",
  fraction = NULL,
  point = FALSE,
  add = TRUE
)
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

- distance_max:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) threshold on the
  maximum distance beyond the range that should be considered to have a
  high likelihood of containing species occurrences (Default: `Inf`
  `"m"`). Can be set to `NULL` or `0` to indicate that no distance
  should be calculated.

- family:

  A [`character`](https://rdrr.io/r/base/character.html) denoting the
  type of model to which this offset is to be added. By default it
  assumes a `'poisson'` distributed model and as a result the output
  created by this function will be log-transformed. If however a
  `'binomial'` distribution is chosen, than the output will be
  `` `logit` `` transformed. For integrated models leave at default.

- presence_prop:

  [`numeric`](https://rdrr.io/r/base/numeric.html) giving the proportion
  of all records expected to be inside the range. By default this is set
  to `0.9` indicating that 10% of all records are likely outside the
  range.

- distance_clip:

  [`logical`](https://rdrr.io/r/base/logical.html) as to whether
  distance should be clipped after the maximum distance (Default:
  `FALSE`).

- distance_function:

  A [`character`](https://rdrr.io/r/base/character.html) specifying the
  distance function to be used. Available are linear (`"linear"`),
  negative exponential kernels (`"negexp"`, default) and a five
  parameters logistic curve (`"logcurve"`) as proposed by Merow et al.
  2017.

- field_occurrence:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) or
  [`character`](https://rdrr.io/r/base/character.html) location of
  biodiversity point records.

- fraction:

  An optional
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object that is multiplied with digitized raster layer. Can be used to
  for example to remove or reduce the expected value (Default: `NULL`).

- point:

  An optional
  [`sf::sf`](https://r-spatial.github.io/sf/reference/sf.html) layer
  with points or [`logical`](https://rdrr.io/r/base/logical.html)
  argument. In the case of the latter the point data is ignored
  (Default: `FALSE`).

- add:

  [`logical`](https://rdrr.io/r/base/logical.html) specifying whether
  new offset is to be added. Setting this parameter to `FALSE` replaces
  the current offsets with the new one (Default: `TRUE`).

## Value

Adds a range offset to a
[`distribution`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
object.

## Details

The output created by this function creates a
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
to be added to a provided distribution object. Offsets in regression
models are likelihood specific as they are added directly to the overall
estimate of `` `y^hat` ``.

Note that all offsets created by this function are by default
log-transformed before export. Background values (e.g. beyond
`"distance_max"`) are set to a very small constant (`1e-10`).

## References

- Merow, C., Wilson, A.M., Jetz, W., 2017. Integrating occurrence data
  and expert maps for improved species range predictions. Glob. Ecol.
  Biogeogr. 26, 243–258. https://doi.org/10.1111/geb.12539

- Merow, C., Allen, J.M., Aiello-Lammens, M., Silander, J.A., 2016.
  Improving niche and range estimates with Maxent and point process
  models by integrating spatially explicit information. Glob. Ecol.
  Biogeogr. 25, 1022–1036. https://doi.org/10.1111/geb.12453

## See also

`"bossMaps"`

Other offset:
[`add_offset()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset.md),
[`add_offset_bias()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset_bias.md),
[`add_offset_elevation()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset_elevation.md),
[`rm_offset()`](https://iiasa.github.io/ibis.iSDM/reference/rm_offset.md)

## Examples

``` r
if (FALSE) { # \dontrun{
 # Train a presence-only model with a simple offset
 fit <- distribution(background) |>
 add_biodiversity_poipo(virtual_points, field_occurrence = "Observed") |>
 add_predictors(predictors) |>
 add_offset_range(virtual_range, distance_max = 5,distance_function = "logcurve",
 distance_clip = TRUE ) |>
 engine_glm() |>
 train()
} # }
```
