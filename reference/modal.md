# Calculate the mode of a provided vector

Calculate the mode of a provided vector

## Usage

``` r
modal(x, na.rm = TRUE)
```

## Arguments

- x:

  A [`vector`](https://rdrr.io/r/base/vector.html) of values or
  characters.

- na.rm:

  [`logical`](https://rdrr.io/r/base/logical.html) whether `NA` values
  are to be removed (Default: `TRUE`)

## Value

The most common (mode) estimate.

## Examples

``` r
# Example
modal(trees$Girth)
#> Error: unable to find an inherited method for function ‘modal’ for signature ‘x = "numeric"’
```
