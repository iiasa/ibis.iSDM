# Functionality for geographic and environmental thinning

For most species distribution modelling approaches it is assumed that
occurrence records are unbiased, which is rarely the case. While
model-based control can alleviate some of the effects of sampling bias,
it can often be desirable to account for some sampling biases through
spatial thinning (Aiello‐Lammens et al. 2015). This is an approach based
on the assumption that over-sampled grid cells contribute little more
than bias, rather than strengthening any environmental responses. This
function provides some methods to apply spatial thinning approaches.
Note that this effectively removes data prior to any estimation and its
use should be considered with care (see also Steen et al. 2021).

## Usage

``` r
thin_observations(
  data,
  background,
  env = NULL,
  method = "random",
  remainpoints = 10,
  mindistance = NULL,
  zones = NULL,
  probs = 0.75,
  global = TRUE,
  centers = NULL,
  verbose = TRUE
)
```

## Arguments

- data:

  A [`sf::sf`](https://r-spatial.github.io/sf/reference/sf.html) object
  with observed occurrence points. All methods threat presence-only and
  presence-absence occurrence points equally.

- background:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object with the background of the study region. Use for assessing
  point density.

- env:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object with environmental covariates. Needed when method is set to
  `"environmental"` or `"bias"` (Default: `NULL`).

- method:

  A [`character`](https://rdrr.io/r/base/character.html) of the method
  to be applied (Default: `"random"`).

- remainpoints:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) giving the number
  of data points at minimum to remain (Default: `10`).

- mindistance:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) for the minimum
  distance of neighbouring observations (Default: `NULL`).

- zones:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  to be supplied when option `"zones"` is chosen (Default: `NULL`).

- probs:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) used as quantile
  threshold in `"bias"` method. (Default: `0.75`).

- global:

  A [`logical`](https://rdrr.io/r/base/logical.html) if during `"bias"`
  method global (entire `env` raster) or local (extracted at point
  locations) bias values are used as for quantile threshold. (Default:
  `TRUE`).

- centers:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) used as number of
  centers for `"environmental"` method. (Default: `NULL`). If not set,
  automatically set to three or nlayers - 1 (whatever is bigger).

- verbose:

  [`logical`](https://rdrr.io/r/base/logical.html) of whether to print
  some statistics about the thinning outcome (Default: `TRUE`).

## Details

All methods only remove points from "over-sampled" grid cells/areas.
These are defined as all cells/areas which either have more points than
`remainpoints` or more points than the global minimum point count per
cell/area (whichever is larger).

Currently implemented thinning methods:

- `"random"`: Samples at random across all over-sampled grid cells
  returning only `"remainpoints"` from over-sampled cells. Does not
  account for any spatial or environmental distance between
  observations.

- `"bias"`: This option removes explicitly points that are considered
  biased only (based on `"env"`). Points are only thinned from grid
  cells which are above the bias quantile (larger values equals greater
  bias). Thins the observations returning `"remainpoints"` from each
  over-sampled and biased cell.

- `"zones"`: Thins observations from each zone that is above the
  over-sampled threshold and returns `"remainpoints"` for each zone.
  Careful: If the zones are relatively wide this can remove quite a few
  observations.

- `"environmental"`: This approach creates an observation-wide
  clustering (k-means) under the assumption that the full environmental
  niche has been comprehensively sampled and is covered by the provided
  covariates `env`. For each over-sampled cluster, we then obtain
  (`"remainpoints"`) by thinning points.

- `"spatial"`: Calculates the spatial distance between all observations.
  Then points are removed iteratively until the minimum distance between
  points is crossed. The `"mindistance"` parameter has to be set for
  this function to work.

## References

- Aiello‐Lammens, M. E., Boria, R. A., Radosavljevic, A., Vilela, B., &
  Anderson, R. P. (2015). spThin: an R package for spatial thinning of
  species occurrence records for use in ecological niche models.
  Ecography, 38(5), 541-545.

- Steen, V. A., Tingley, M. W., Paton, P. W., & Elphick, C. S. (2021).
  Spatial thinning and class balancing: Key choices lead to variation in
  the performance of species distribution models with citizen science
  data. Methods in Ecology and Evolution, 12(2), 216-226.

## Examples

``` r
if (FALSE) { # \dontrun{
 # Thin a certain number of observations
 # At random
 thin_points <- thin_observations(points, background, method = "random")
 # using a bias layer
 thin_points <- thin_observations(points, background, method = "bias", env = bias)
} # }
```
