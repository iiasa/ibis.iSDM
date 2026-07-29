# Tests if an input is a stars object.

Tests if an input is a stars object.

## Usage

``` r
is.stars(x)
```

## Arguments

- x:

  an R Object.

## Value

Boolean evaluation with [logical](https://rdrr.io/r/base/logical.html)
output.

## Examples

``` r
x <- stars::st_as_stars(matrix(1, nrow = 1, ncol = 1))
is.stars(x)
#> [1] TRUE
```
