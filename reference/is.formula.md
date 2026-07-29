# Check whether a formula is valid

Check whether a formula is valid

## Usage

``` r
is.formula(x)
```

## Arguments

- x:

  A [`character`](https://rdrr.io/r/base/character.html) object

## Value

Boolean evaluation with [logical](https://rdrr.io/r/base/logical.html)
output.

## Examples

``` r
f <- y ~ x
is.formula(f)
#> [1] TRUE
```
