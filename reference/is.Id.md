# Check whether a provided object is truly of a specific type

Check whether a provided object is truly of a specific type

## Usage

``` r
is.Id(x)
```

## Arguments

- x:

  A provided Id object

## Value

Boolean evaluation with [logical](https://rdrr.io/r/base/logical.html)
output.

## Examples

``` r
id <- as.Id("example-id")
is.Id(id)
#> [1] TRUE
```
