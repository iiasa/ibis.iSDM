# Function to add GLOBIOM-DownScalr derived predictors to a Biodiversity distribution object

**This function is defunct! Use the BNRTools package for formatting
input data!**

## Usage

``` r
add_predictors_globiom(x, ...)

# S4 method for class 'BiodiversityDistribution'
add_predictors_globiom(x, ...)

# S4 method for class 'BiodiversityScenario'
add_predictors_globiom(x, ...)
```

## Arguments

- x:

  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

- ...:

  Other parameters passed down

## Value

No return value; this function is defunct and always errors.

## Details

See
[`add_predictors()`](https://iiasa.github.io/ibis.iSDM/reference/add_predictors.md)

## References

<https://github.com/iiasa/BNRTools>

## See also

[add_predictors](https://iiasa.github.io/ibis.iSDM/reference/add_predictors.md)

## Examples

``` r
if (FALSE) { # \dontrun{
 obj <- distribution(background) |>
        add_predictors_globiom(fname = "", transform = 'none')
 obj
} # }
```
