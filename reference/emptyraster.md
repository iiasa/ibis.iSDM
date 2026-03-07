# Create an empty `SpatRaster` based on a template

This function creates an empty copy of a provided `SpatRaster` object.
It is primarily used in the package to create the outputs for the
predictions.

## Usage

``` r
emptyraster(x, res = NULL, ...)
```

## Arguments

- x:

  A `SpatRaster*`, `stars` or a `sf` object from which coordinates can
  be obtained. Note that for `sf` objects the parameter `res` needs to
  be supplied.

- res:

  (Optional) [`numeric`](https://rdrr.io/r/base/numeric.html) estimate
  on the resolution of the output (Default: `NULL`).

- ...:

  other arguments that can be passed to
  [`terra`](https://rspatial.github.io/terra/reference/terra-package.html)

## Value

an empty
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html),
i.e. all cells are `NA`.

## Examples

``` r
require(terra)
#> Loading required package: terra
#> terra 1.8.93
#> 
#> Attaching package: ‘terra’
#> The following object is masked from ‘package:ibis.iSDM’:
#> 
#>     modal
r <- rast(matrix(1:100, 5, 20))
emptyraster(r)
#> class       : SpatRaster 
#> size        : 5, 20, 1  (nrow, ncol, nlyr)
#> resolution  : 1, 1  (x, y)
#> extent      : 0, 20, 0, 5  (xmin, xmax, ymin, ymax)
#> coord. ref. :  
```
