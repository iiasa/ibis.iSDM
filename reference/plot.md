# Plot wrappers

Plots information from a given object where a plotting object is
available.

## Usage

``` r
# S3 method for class 'DistributionModel'
plot(x, what = "mean", ...)

# S3 method for class 'BiodiversityDatasetCollection'
plot(x, ...)

# S3 method for class 'PredictorDataset'
plot(x, ...)

# S3 method for class 'Engine'
plot(x, ...)

# S3 method for class 'BiodiversityScenario'
plot(x, ...)
```

## Arguments

- x:

  Any object belonging to
  [DistributionModel](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md),
  [BiodiversityDatasetCollection](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDatasetCollection-class.md),
  [PredictorDataset](https://iiasa.github.io/ibis.iSDM/reference/PredictorDataset-class.md)
  or
  [BiodiversityScenario](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md).

- what:

  In case a
  [terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  is supplied, this parameter specifies the layer to be shown (Default:
  `"mean"`).

- ...:

  Further arguments passed on to `x$plot`.

## Value

Graphical output

## Details

The plotted outputs vary depending on what object is being plotted. For
example for a fitted
[DistributionModel](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
the output is usually the fitted spatial prediction (Default: `'mean'`).

## Examples

``` r
if (FALSE) { # \dontrun{
# Build and train a model
mod <- distribution(background) |>
  add_biodiversity_poipo(species) |>
  add_predictors(predictors) |>
  engine_glmnet() |>
  train()
# Plot the resulting model
plot(mod)
} # }
```
