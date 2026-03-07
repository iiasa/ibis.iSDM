# Remove control from an existing distribution object

This function allows to remove set control options from an existing
[distribution](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
object.

## Usage

``` r
rm_control(x, type)

# S4 method for class 'BiodiversityDistribution'
rm_control(x, type)
```

## Arguments

- x:

  [distribution](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

- type:

  A [`character`](https://rdrr.io/r/base/character.html) vector
  describing the type of control to be removed. Can be missing.

## See also

[`add_control_bias()`](https://iiasa.github.io/ibis.iSDM/reference/add_control_bias.md)

Other control:
[`rm_limits()`](https://iiasa.github.io/ibis.iSDM/reference/rm_limits.md)

## Examples

``` r
if (FALSE) { # \dontrun{
 x <- distribution(background) |>
   add_predictors(covariates) |>
   add_control_bias(method = "proximity")
 x <- x |> rm_control()
 x
} # }
```
