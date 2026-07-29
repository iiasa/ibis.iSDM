# Obtains the coefficients of a trained model

Similar as
[`summary`](https://iiasa.github.io/ibis.iSDM/reference/summary.md),
this helper function obtains the coefficients from a given
[DistributionModel](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
object.

## Usage

``` r
# S3 method for class 'DistributionModel'
coef(object, ...)
```

## Arguments

- object:

  Any prepared object.

- ...:

  not used.

## Note

For models trained with machine-learning approaches (e.g.
[`engine_bart`](https://iiasa.github.io/ibis.iSDM/reference/engine_bart.md)
etc) this function will return variable importance estimates rather than
linear coefficients. Similar can be said for trained non-linear models.

## See also

[`stats::coef()`](https://rdrr.io/r/stats/coef.html).

## Examples

``` r
model <- DistributionModel$new("Example model")
model$get_coefficients <- function() {
  data.frame(variable = "temperature", coefficient = 1)
}
coef(model)
#>      variable coefficient
#> 1 temperature           1
```
