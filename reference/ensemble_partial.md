# Function to create an ensemble of partial effects from multiple models

Similar to the
[`ensemble()`](https://iiasa.github.io/ibis.iSDM/reference/ensemble.md)
function, this function creates an ensemble of partial responses of
provided distribution models fitted with the
[`ibis.iSDM-package`](https://iiasa.github.io/ibis.iSDM/reference/ibis.iSDM.md).
Through the `layer` parameter it can be specified which part of the
partial prediction should be averaged in an ensemble (if given). This
can be for instance the *mean* prediction and/or the standard deviation
*sd*. Ensemble partial is also being called if more than one input
[`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
object is provided to `partial`.

By default the ensemble of partial responses is created as average
across all models with the uncertainty being the standard deviation of
responses.

## Usage

``` r
ensemble_partial(
  ...,
  x.var,
  method = "mean",
  layer = "mean",
  newdata = NULL,
  normalize = TRUE
)

# S4 method for class 'ANY'
ensemble_partial(
  ...,
  x.var,
  method = "mean",
  layer = "mean",
  newdata = NULL,
  normalize = TRUE
)
```

## Arguments

- ...:

  Provided
  [`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
  objects from which partial responses can be called. In the future
  provided data.frames might be supported as well.

- x.var:

  A [`character`](https://rdrr.io/r/base/character.html) of the variable
  from which an ensemble is to be created.

- method:

  Approach on how the ensemble is to be created. See details for options
  (Default: `'mean'`).

- layer:

  A [`character`](https://rdrr.io/r/base/character.html) of the layer to
  be taken from each prediction (Default: `'mean'`). If set to `NULL`
  ignore any of the layer names in ensembles of `SpatRaster` objects.

- newdata:

  A optional [`data.frame`](https://rdrr.io/r/base/data.frame.html) or
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object supplied to the model (DefaultL `NULL`). This object needs to
  have identical names as the original predictors.

- normalize:

  [`logical`](https://rdrr.io/r/base/logical.html) on whether the inputs
  of the ensemble should be normalized to a scale of 0-1 (Default:
  `TRUE`).

## Value

A [data.frame](https://rdrr.io/r/base/data.frame.html) with the combined
partial effects of the supplied models.

## Details

Possible options for creating an ensemble includes:

- `'mean'` - Calculates the mean of several predictions.

- `'median'` - Calculates the median of several predictions.

## Note

If a list is supplied, then it is assumed that each entry in the list is
a fitted
[`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
object. Take care not to create an ensemble of models constructed with
different link functions, e.g. logistic vs
[log](https://rspatial.github.io/terra/reference/math-generics.html). By
default the response functions of each model are normalized.

## Examples

``` r
if (FALSE) { # \dontrun{
 # Assumes previously computed models
 ex <- ensemble_partial(mod1, mod2, mod3, method = "mean")
} # }
```
