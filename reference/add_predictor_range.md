# Add a range of a species as predictor to a distribution object

This function allows to add a species range which is usually drawn by
experts in a separate process as spatial explicit prior. Both
[`sf::sf`](https://r-spatial.github.io/sf/reference/sf.html) and
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)-objects
are supported as input.

Users are advised to look at the `"bossMaps"` R-package presented as
part of Merow et al. (2017), which allows flexible calculation of
non-linear distance transforms from the boundary of the range. Outputs
of this package could be added directly to this function. **Note that
this function adds the range as predictor and not as offset. For this
purpose a separate function
[`add_offset_range()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset_range.md)
exists.**

Additional options allow to include the range either as `"binary"` or as
`"distance"` transformed predictor. The difference being that the range
is either directly included as presence-only predictor or alternatively
with a linear distance transform from the range boundary. The parameter
`"distance_max"` can be specified to constrain this distance transform.

## Usage

``` r
add_predictor_range(
  x,
  layer,
  method = "distance",
  distance_max = NULL,
  fraction = NULL,
  priors = NULL
)

# S4 method for class 'BiodiversityDistribution,SpatRaster'
add_predictor_range(
  x,
  layer,
  method = "precomputed_range",
  fraction = NULL,
  priors = NULL
)

# S4 method for class 'BiodiversityDistribution,sf'
add_predictor_range(
  x,
  layer,
  method = "distance",
  distance_max = Inf,
  fraction = NULL,
  priors = NULL
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

- method:

  [`character`](https://rdrr.io/r/base/character.html) describing how
  the range should be included (`"binary"` \| `"distance"`).

- distance_max:

  Numeric threshold on the maximum distance (Default: `NULL`).

- fraction:

  An optional
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object that is multiplied with digitized raster layer. Can be used to
  for example to remove or reduce the expected value (Default: `NULL`).

- priors:

  A
  [`PriorList`](https://iiasa.github.io/ibis.iSDM/reference/PriorList-class.md)
  object. Default is set to NULL which uses default prior assumptions.

## References

- Merow, C., Wilson, A. M., & Jetz, W. (2017). Integrating occurrence
  data and expert maps for improved species range predictions. Global
  Ecology and Biogeography, 26(2), 243–258.
  https://doi.org/10.1111/geb.12539

## Examples

``` r
if (FALSE) { # \dontrun{
distribution(background) |>
  add_predictor_range(range, method = "distance", distance_max = 2)
} # }
```
