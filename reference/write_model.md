# Save a model for later use

The `write_model` function (opposed to the `write_output`) is a generic
wrapper to writing a
[`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
to disk. It is essentially a wrapper to
[`saveRDS`](https://rspatial.github.io/terra/reference/serialize.html).
Models can be loaded again via the `load_model` function.

## Usage

``` r
write_model(
  mod,
  fname,
  slim = FALSE,
  verbose = getOption("ibis.setupmessages", default = TRUE)
)

# S4 method for class 'ANY'
write_model(
  mod,
  fname,
  slim = FALSE,
  verbose = getOption("ibis.setupmessages", default = TRUE)
)
```

## Arguments

- mod:

  Provided
  [`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
  object.

- fname:

  A [`character`](https://rdrr.io/r/base/character.html) depicting an
  output filename.

- slim:

  A [`logical`](https://rdrr.io/r/base/logical.html) option to whether
  unnecessary entries in the model object should be deleted. This
  deletes for example predictions or any other non-model content from
  the object (Default: `FALSE`).

- verbose:

  [`logical`](https://rdrr.io/r/base/logical.html) indicating whether
  messages should be shown. Overwrites `getOption("ibis.setupmessages")`
  (Default: `TRUE`).

## Value

No R-output is created. A file is written to the target direction.

## Note

By default output files will be overwritten if already existing!

## See also

load_model

## Examples

``` r
if (FALSE) { # \dontrun{
x <- distribution(background) |>
 add_biodiversity_poipo(virtual_points, field_occurrence = 'observed', name = 'Virtual points') |>
 add_predictors(pred_current, transform = 'scale',derivates = 'none') |>
 engine_xgboost(nrounds = 2000) |> train(varsel = FALSE, only_linear = TRUE)
write_model(x, "testmodel.rds")
} # }
```
