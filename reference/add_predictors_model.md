# Add predictions from a fitted model to a Biodiversity distribution object

This function is a convenience wrapper to add the output from a previous
fitted
[`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
to another
[`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md)
object. Obviously only works if a prediction was fitted in the model.
Options to instead add thresholds, or to transform / derivate the model
outputs are also supported.

## Usage

``` r
add_predictors_model(
  x,
  model,
  transform = "scale",
  derivates = "none",
  threshold_only = FALSE,
  priors = NULL,
  ...
)

# S4 method for class 'BiodiversityDistribution'
add_predictors_model(
  x,
  model,
  transform = "scale",
  derivates = "none",
  threshold_only = FALSE,
  priors = NULL,
  ...
)
```

## Arguments

- x:

  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

- model:

  A
  [`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
  object.

- transform:

  A [`vector`](https://rdrr.io/r/base/vector.html) stating whether
  predictors should be preprocessed in any way (Options:
  `'none'`,`'pca'`, `'scale'`, `'norm'`)

- derivates:

  A Boolean check whether derivate features should be considered
  (Options: `'none'`, `'thresh'`, `'hinge'`, `'quad'`) )

- threshold_only:

  A [`logical`](https://rdrr.io/r/base/logical.html) flag indicating
  whether to add thresholded layers from the fitted model (if existing)
  instead (Default: `FALSE`).

- priors:

  A
  [`PriorList`](https://iiasa.github.io/ibis.iSDM/reference/PriorList-class.md)
  object. Default is set to `NULL` which uses default prior assumptions.

- ...:

  Other parameters passed down

## Details

A transformation takes the provided rasters and for instance rescales
them or transforms them through a principal component analysis
([prcomp](https://rdrr.io/r/stats/prcomp.html)). In contrast, derivates
leave the original provided predictors alone, but instead create new
ones, for instance by transforming their values through a quadratic or
hinge transformation. Note that this effectively increases the number of
predictors in the object, generally requiring stronger regularization by
the used
[`Engine`](https://iiasa.github.io/ibis.iSDM/reference/Engine-class.md).
Both transformations and derivates can also be combined. Available
options for transformation are:

- `'none'` - Leaves the provided predictors in the original scale.

- `'pca'` - Converts the predictors to principal components. Note that
  this results in a renaming of the variables to principal component
  axes!

- `'scale'` - Transforms all predictors by applying
  [scale](https://rdrr.io/r/base/scale.html) on them.

- `'norm'` - Normalizes all predictors by transforming them to a scale
  from 0 to 1.

- `'windsor'` - Applies a windsorization to the target predictors. By
  default this effectively cuts the predictors to the 0.05 and 0.95,
  thus helping to remove extreme outliers.

Available options for creating derivates are:

- `'none'` - No additional predictor derivates are created.

- `'quad'` - Adds quadratic transformed predictors.

- `'interaction'` - Add interacting predictors. Interactions need to be
  specified (`"int_variables"`)!

- `'thresh'` - Add threshold transformed predictors.

- `'hinge'` - Add hinge transformed predictors.

- `'bin'` - Add predictors binned by their percentiles.

## Examples

``` r
if (FALSE) { # \dontrun{
 # Fit first model
 fit <- distribution(background) |>
        add_predictors(covariates) |>
        add_biodiversity_poipa(species) |>
        engine_glmnet() |>
        train()

 # New model object
 obj <- distribution(background) |>
        add_predictors_model(fit)
 obj
} # }
```
