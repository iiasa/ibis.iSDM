# Plot effects of trained model

This functions is handy wrapper that calls the default plotting
functions for the model of a specific engine. Equivalent to calling
`effects` of a fitted
[distribution](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
function.

## Usage

``` r
# S3 method for class 'DistributionModel'
effects(object, ...)
```

## Arguments

- object:

  Any fitted
  [distribution](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  object.

- ...:

  Not used.

## Value

None.

## Note

For some models, where default coefficients plots are not available,
this function will attempt to generate
[partial](https://iiasa.github.io/ibis.iSDM/reference/partial.md)
dependency plots instead.

## Examples

``` r
if (FALSE) { # \dontrun{
# Where mod is an estimated distribution model
mod$effects()
} # }
```
