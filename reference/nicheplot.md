# Niche plot for distribution objects

The suitability of any given area for a biodiversity feature can in many
instances be complex and non-linear. Visualizing obtained suitability
predictions (e.g. from
[`train()`](https://iiasa.github.io/ibis.iSDM/reference/train.md))
against underlying predictors might help to explain the underlying
gradients of the niche.

Supported Inputs for this function are either single trained `ibis.iSDM`
[`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
objects or alternatively a set of three
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
objects. In both cases, users can specify `"xvar"` and `"yvar"`
explicitly or leave them empty. In the latter case a principal component
analysis (PCA) is conducted on the full environmental stack (loaded from
[`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
or supplied separately).

## Usage

``` r
nicheplot(
  mod,
  xvar = NULL,
  yvar = NULL,
  envvars = NULL,
  overlay_data = FALSE,
  plot = TRUE,
  fname = NULL,
  title = NULL,
  pal = NULL,
  ...
)

# S4 method for class 'ANY'
nicheplot(
  mod,
  xvar = NULL,
  yvar = NULL,
  envvars = NULL,
  overlay_data = FALSE,
  plot = TRUE,
  fname = NULL,
  title = NULL,
  pal = NULL,
  ...
)
```

## Arguments

- mod:

  A trained
  [`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
  or alternatively a
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object with `prediction` model within.

- xvar:

  A [`character`](https://rdrr.io/r/base/character.html) denoting the
  predictor on the x-axis. Alternatively a
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object can be provided.

- yvar:

  A [`character`](https://rdrr.io/r/base/character.html) denoting the
  predictor on the y-axis. Alternatively a
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object can be provided.

- envvars:

  A
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object containing all environmental variables. Only used if `xvar` and
  `yvar` is empty (Default: `NULL`).

- overlay_data:

  A [`logical`](https://rdrr.io/r/base/logical.html) on whether training
  data should be overlaid on the plot. Only used for
  [`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
  objects (Default: `FALSE`).

- plot:

  A [`logical`](https://rdrr.io/r/base/logical.html) indication of
  whether the result is to be plotted (Default: `TRUE`)?

- fname:

  A [`character`](https://rdrr.io/r/base/character.html) specifying the
  output file name a created figure should be written to.

- title:

  Allows to respecify the title through a
  [`character`](https://rdrr.io/r/base/character.html) (Default:
  `NULL`).

- pal:

  An optional [`vector`](https://rdrr.io/r/base/vector.html) with
  continuous custom colours (Default: `NULL`).

- ...:

  Other engine specific parameters.

## Value

Saved niche plot in `'fname'` if specified, otherwise plot.

## See also

[partial](https://iiasa.github.io/ibis.iSDM/reference/partial.md),
[plot.DistributionModel](https://iiasa.github.io/ibis.iSDM/reference/plot.md)

## Examples

``` r
# Make quick prediction
background <- terra::rast(system.file('extdata/europegrid_50km.tif',
package='ibis.iSDM',mustWork = TRUE))
virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM'), 'points',quiet = TRUE)
ll <- list.files(system.file('extdata/predictors/',package = 'ibis.iSDM',mustWork = TRUE),full.names = TRUE)

# Load them as rasters
predictors <- terra::rast(ll);names(predictors) <- tools::file_path_sans_ext(basename(ll))

# Add GLM as an engine and predict
fit <- distribution(background) |>
add_biodiversity_poipo(virtual_points, field_occurrence = 'Observed',
name = 'Virtual points',docheck = FALSE) |>
add_predictors(predictors, transform = 'none',derivates = 'none') |>
engine_glm() |>
train()
#> [Setup] 2026-05-14 20:28:36.762003 | Creating distribution object...
#> [Setup] 2026-05-14 20:28:36.763088 | Adding poipo dataset...
#> [Setup] 2026-05-14 20:28:36.768667 | Adding predictors...
#> [Estimation] 2026-05-14 20:28:36.944479 | Collecting input parameters.
#> [Estimation] 2026-05-14 20:28:37.11298 | Adding engine-specific parameters.
#> [Estimation] 2026-05-14 20:28:37.117973 | Engine setup.
#> [Estimation] 2026-05-14 20:28:37.373939 | Starting fitting: Virtual points
#> [Estimation] 2026-05-14 20:28:37.435185 | Starting prediction...
#> [Done] 2026-05-14 20:28:37.564561 | Completed after 0.62 secs

# Plot niche for prediction for temperature and forest cover
nicheplot(fit, xvar = "bio01_mean_50km", yvar = "CLC3_312_mean_50km" )
```
