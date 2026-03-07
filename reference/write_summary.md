# Generic function to write summary outputs from created models.

The `write_summary` function is a wrapper function to create summaries
from fitted
[`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
or
[`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
objects. This function will extract parameters and statistics about the
used data from the input object and writes the output as either `'rds'`
or `'rdata'` file. Alternative, more open file formats are under
consideration.

## Usage

``` r
write_summary(
  mod,
  fname,
  partial = FALSE,
  verbose = getOption("ibis.setupmessages", default = TRUE),
  ...
)

# S4 method for class 'ANY,character'
write_summary(
  mod,
  fname,
  partial = FALSE,
  verbose = getOption("ibis.setupmessages", default = TRUE),
  ...
)
```

## Arguments

- mod:

  Provided
  [`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
  or
  [`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
  object.

- fname:

  A [`character`](https://rdrr.io/r/base/character.html) depicting an
  output filename. The suffix determines the file type of the output
  (Options: `'rds'`, `'rdata'`).

- partial:

  A [`logical`](https://rdrr.io/r/base/logical.html) value determining
  whether partial variable contributions should be calculated and added
  to the model summary. **Note** that this can be rather slow (Default:
  `FALSE`).

- verbose:

  [`logical`](https://rdrr.io/r/base/logical.html) indicating whether
  messages should be shown. Overwrites `getOption("ibis.setupmessages")`
  (Default: `TRUE`).

- ...:

  Any other arguments passed on the individual functions.

## Value

No R-output is created. A file is written to the target direction.

## Note

No predictions or tabular data is saved through this function. Use
[`write_output()`](https://iiasa.github.io/ibis.iSDM/reference/write_output.md)
to save those.

## Examples

``` r
if (FALSE) { # \dontrun{
x <- distribution(background) |>
 add_biodiversity_poipo(virtual_points, field_occurrence = 'observed', name = 'Virtual points')  |>
 add_predictors(pred_current, transform = 'scale',derivates = 'none') |>
 engine_xgboost(nrounds = 2000) |> train(varsel = FALSE, only_linear = TRUE)
write_summary(x, "testmodel.rds")
} # }
```
