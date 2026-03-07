# Generic function to write spatial outputs

The `write_output` function is a generic wrapper to writing any output
files (e.g. projections) created with the
[`ibis.iSDM-package`](https://iiasa.github.io/ibis.iSDM/reference/ibis.iSDM.md).
It is possible to write outputs of fitted
[`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md),
[`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
or individual
[`terra::terra`](https://rspatial.github.io/terra/reference/terra-package.html)
or [`stars`](https://rdrr.io/r/graphics/stars.html) objects. In case a
[`data.frame`](https://rdrr.io/r/base/data.frame.html) is supplied, the
output is written as csv file. **For creating summaries of distribution
and scenario parameters and performance, see
[`write_summary()`](https://iiasa.github.io/ibis.iSDM/reference/write_summary.md)**

## Usage

``` r
write_output(
  mod,
  fname,
  dt = "FLT4S",
  verbose = getOption("ibis.setupmessages", default = TRUE),
  ...
)

# S4 method for class 'ANY,character'
write_output(
  mod,
  fname,
  dt = "FLT4S",
  verbose = getOption("ibis.setupmessages", default = TRUE),
  ...
)

# S4 method for class 'BiodiversityScenario,character'
write_output(
  mod,
  fname,
  dt = "FLT4S",
  verbose = getOption("ibis.setupmessages", default = TRUE),
  ...
)

# S4 method for class 'SpatRaster,character'
write_output(
  mod,
  fname,
  dt = "FLT4S",
  verbose = getOption("ibis.setupmessages", default = TRUE),
  ...
)

# S4 method for class 'data.frame,character'
write_output(
  mod,
  fname,
  dt = "FLT4S",
  verbose = getOption("ibis.setupmessages", default = TRUE),
  ...
)

# S4 method for class 'stars,character'
write_output(
  mod,
  fname,
  dt = "FLT4S",
  verbose = getOption("ibis.setupmessages", default = TRUE),
  ...
)
```

## Arguments

- mod:

  Provided
  [`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md),
  [`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md),
  [`terra::terra`](https://rspatial.github.io/terra/reference/terra-package.html)
  or [`stars`](https://rdrr.io/r/graphics/stars.html) object.

- fname:

  A [`character`](https://rdrr.io/r/base/character.html) depicting an
  output filename.

- dt:

  A [`character`](https://rdrr.io/r/base/character.html) for the output
  datatype. Following the
  [`terra::writeRaster`](https://rspatial.github.io/terra/reference/writeRaster.html)
  options (Default: `'FLT4S'`).

- verbose:

  [`logical`](https://rdrr.io/r/base/logical.html) indicating whether
  messages should be shown. Overwrites `getOption("ibis.setupmessages")`
  (Default: `TRUE`).

- ...:

  Any other arguments passed on the individual functions.

## Value

No R-output is created. A file is written to the target direction.

## Note

By default output files will be overwritten if already existing!

## Examples

``` r
if (FALSE) { # \dontrun{
x <- distribution(background)  |>
 add_biodiversity_poipo(virtual_points, field_occurrence = 'observed', name = 'Virtual points') |>
 add_predictors(pred_current, transform = 'scale',derivates = 'none') |>
 engine_xgboost(nrounds = 2000) |> train(varsel = FALSE, only_linear = TRUE)
write_output(x, "testmodel.tif")
} # }
```
