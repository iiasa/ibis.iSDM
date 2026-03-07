# Align a [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html) object to another by harmonizing geometry and extend.

If the data is not in the same projection as the template, the alignment
will be computed by reprojection only. If the data has already the same
projection, the data set will be cropped and aggregated prior to
resampling in order to reduce computation time.

## Usage

``` r
alignRasters(data, template, method = "bilinear", func = mean, cl = TRUE)
```

## Arguments

- data:

  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object to be resampled.

- template:

  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  or [`sf::sf`](https://r-spatial.github.io/sf/reference/sf.html) object
  from which geometry can be extracted.

- method:

  method for resampling (Options: `"near"` or `"bilinear"`).

- func:

  function for resampling (Default:
  [mean](https://rdrr.io/r/base/mean.html)).

- cl:

  [`logical`](https://rdrr.io/r/base/logical.html) value if multicore
  computation should be used (Default: `TRUE`).

## Value

New
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
object aligned to the supplied template layer.

## Details

Nearest Neighbour resampling (near) is recommended for discrete and
bilinear resampling recommended for continuous data. See also help from
[terra::resample](https://rspatial.github.io/terra/reference/resample.html)
for other options.

## Examples

``` r
if (FALSE) { # \dontrun{
 # Align one raster to another
 ras1 <- alignRasters( ras1, ras2, method = "near", cl = FALSE)
} # }
```
