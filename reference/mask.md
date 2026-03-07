# Mask data with an external layer

This is a helper function that takes an existing object created by the
ibis.iSDM package and an external layer, then intersects both. It
currently takes either a
[DistributionModel](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md),
[BiodiversityDatasetCollection](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDatasetCollection-class.md),
[PredictorDataset](https://iiasa.github.io/ibis.iSDM/reference/PredictorDataset-class.md)
or
[BiodiversityScenario](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
as input.

As mask either a
[`sf::sf`](https://r-spatial.github.io/sf/reference/sf.html) or
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
object can be chosen. The mask will be converted internally depending on
the object.

## Usage

``` r
mask.DistributionModel(x, mask, inverse = FALSE, ...)

mask.BiodiversityDatasetCollection(x, mask, inverse = FALSE, ...)

mask.PredictorDataset(x, mask, inverse = FALSE, ...)

mask.BiodiversityScenario(x, mask, inverse = FALSE, ...)
```

## Arguments

- x:

  Any object belonging to
  [DistributionModel](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md),
  [BiodiversityDatasetCollection](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDatasetCollection-class.md),
  [PredictorDataset](https://iiasa.github.io/ibis.iSDM/reference/PredictorDataset-class.md)
  or
  [BiodiversityScenario](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md).

- mask:

  A [`sf::sf`](https://r-spatial.github.io/sf/reference/sf.html) or
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object.

- inverse:

  A [`logical`](https://rdrr.io/r/base/logical.html) flag whether to
  take inverse of the mask instead (Default: `FALSE`).

- ...:

  Passed on arguments

## Value

A respective object of the input type.

## See also

[`terra::mask()`](https://rspatial.github.io/terra/reference/mask.html)

## Examples

``` r
if (FALSE) { # \dontrun{
# Build and train a model
mod <- distribution(background) |>
  add_biodiversity_poipo(species) |>
  add_predictors(predictors) |>
  engine_glmnet() |>
  train()

# Constrain the prediction by another object
mod <- mask(mod, speciesrange)

} # }
```
