# Tests if an input is a SpatRaster object.

Tests if an input is a SpatRaster object.

## Usage

``` r
is.Raster(x)
```

## Arguments

- x:

  an R Object.

## Value

Boolean evaluation with [logical](https://rdrr.io/r/base/logical.html)
output.

## Examples

``` r
r <- terra::rast(nrows = 1, ncols = 1, vals = 1)
is.Raster(r)
#> [1] TRUE
```
