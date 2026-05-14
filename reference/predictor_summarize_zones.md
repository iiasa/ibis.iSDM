# Summarize and replace values in predictors with an aggregation

This function supports the aggregation and specification of values
inside a predictor set. Concretely, all predictor values underneath each
zonal polygon (or id) are summarized through a function and then set as
new values to the predictor set.

## Usage

``` r
predictor_summarize_zones(
  env,
  layer,
  fun = "mean",
  target_layer = NULL,
  value = NULL
)
```

## Arguments

- env:

  A
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object with the predictors.

- layer:

  Either a
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  or a [`sf`](https://r-spatial.github.io/sf/reference/sf.html) object.

- fun:

  A custom function or
  [`character`](https://rdrr.io/r/base/character.html). A range of
  methods are supported by default such as `"mean"`, `"sum"`, `"min"` or
  `"max"` (Default: `"mean"`).

- target_layer:

  (Optional)
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object with categorical zones for which the summarized value is to be
  set.

- value:

  (Optional) [`numeric`](https://rdrr.io/r/base/numeric.html) value
  instead of the layer.

## Value

A
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
object with the same number of layers as the input.

## Details

Potential uses of this function include for example the summary of
values within a set of regions (such as countries or super cells).
Optional target layers (for example another zone) and values can be
defined.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load example data
background <- terra::rast(system.file('extdata/europegrid_50km.tif',
package='ibis.iSDM',mustWork = TRUE))

co <- terra::rast(system.file('extdata/countries.tif',
package='ibis.iSDM',mustWork = TRUE))

# Get the UK
co <- co=="United Kingdom"
co[co==0] <- NA

# Get list of test predictors
ll <- list.files(system.file('extdata/predictors/', package = 'ibis.iSDM',
mustWork = TRUE),full.names = TRUE)
# Load them as rasters
env <- terra::rast(ll[3])

# Load the layer
env_r <- predictor_summarize_zones(env, layer = co, fun = 'mean')
terra::plot(c(env, env_r), main = c("original", "aggregated"))
} # }
```
