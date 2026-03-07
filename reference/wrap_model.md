# Wrap a model for later use

The `wrap_model` function uses
[`terra::wrap()`](https://rspatial.github.io/terra/reference/wrap.html)
to easier ship a `DistributionModel` object.

## Usage

``` r
wrap_model(mod, verbose = getOption("ibis.setupmessages", default = TRUE))

# S4 method for class 'ANY'
wrap_model(mod, verbose = getOption("ibis.setupmessages", default = TRUE))
```

## Arguments

- mod:

  Provided
  [`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
  object.

- verbose:

  [`logical`](https://rdrr.io/r/base/logical.html) indicating whether
  messages should be shown. Overwrites `getOption("ibis.setupmessages")`
  (Default: `TRUE`).

## Value

DistributionModel with wrapped raster layers

## See also

unwrap_model

## Examples

``` r
if (FALSE) { # \dontrun{
x <- distribution(background) |>
 add_biodiversity_poipo(virtual_points, field_occurrence = 'observed', name = 'Virtual points') |>
 add_predictors(pred_current, transform = 'scale',derivates = 'none') |>
 engine_xgboost(nrounds = 2000) |>
 train(varsel = FALSE, only_linear = TRUE)
wrap_model(x, "testmodel.rds")
} # }
```
