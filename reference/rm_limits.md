# Remove limits from an existing distribution object

This function allows to remove set limits from an existing
[distribution](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
object.

## Usage

``` r
rm_limits(x)

# S4 method for class 'BiodiversityDistribution'
rm_limits(x)
```

## Arguments

- x:

  [distribution](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

## Value

A
[`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md)
object with extrapolation limits removed.

## See also

[`add_limits_extrapolation()`](https://iiasa.github.io/ibis.iSDM/reference/add_limits_extrapolation.md)

Other control:
[`rm_control()`](https://iiasa.github.io/ibis.iSDM/reference/rm_control.md)

## Examples

``` r
if (FALSE) { # \dontrun{
 x <- distribution(background) |>
   add_predictors(covariates) |>
   add_limits_extrapolation(method = "zones", layer = zones)
 x <- x |> rm_limits()
 x
} # }
```
