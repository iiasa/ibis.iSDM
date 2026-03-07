# Obtain partial effects of trained model

Create a partial response or effect plot of a trained model.

## Usage

``` r
partial(
  mod,
  x.var = NULL,
  constant = NULL,
  variable_length = 100,
  values = NULL,
  newdata = NULL,
  plot = FALSE,
  type = "response",
  ...
)

# S4 method for class 'ANY'
partial(
  mod,
  x.var = NULL,
  constant = NULL,
  variable_length = 100,
  values = NULL,
  newdata = NULL,
  plot = FALSE,
  type = "response",
  ...
)

partial.DistributionModel(mod, ...)
```

## Arguments

- mod:

  A trained `DistributionModel` object with `fit_best` model within.

- x.var:

  A [character](https://rdrr.io/r/base/character.html) indicating the
  variable for which a partial effect is to be calculated.

- constant:

  A [numeric](https://rdrr.io/r/base/numeric.html) constant to be
  inserted for all other variables. Default calculates a mean per
  variable.

- variable_length:

  [numeric](https://rdrr.io/r/base/numeric.html) The interpolation depth
  (nr. of points) to be used (Default: `100`).

- values:

  [numeric](https://rdrr.io/r/base/numeric.html) Directly specified
  values to compute partial effects for. If this parameter is set to
  anything other than `NULL`, the parameter `"variable_length"` is
  ignored (Default: `NULL`).

- newdata:

  An optional [data.frame](https://rdrr.io/r/base/data.frame.html) with
  provided data for partial estimation (Default: `NULL`).

- plot:

  A [`logical`](https://rdrr.io/r/base/logical.html) indication of
  whether the result is to be plotted?

- type:

  A specified type, either `'response'` or `'predictor'`. Can be
  missing.

- ...:

  Other engine specific parameters.

## Value

A [data.frame](https://rdrr.io/r/base/data.frame.html) with the created
partial response.

## Details

By default the mean is calculated across all parameters that are not
`x.var`. Instead a *constant* can be set (for instance `0`) to be
applied to the output.

## See also

partial

## Examples

``` r
if (FALSE) { # \dontrun{
 # Do a partial calculation of a trained model
 partial(fit, x.var = "Forest.cover", plot = TRUE)
} # }
```
