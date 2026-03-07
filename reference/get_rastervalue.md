# Function to extract point values directly from a SpatRaster

This function simply extracts the values from a provided
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html),
[`terra::SpatRasterDataset`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
or
[`terra::SpatRasterCollection`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
object. For points where or `NA` values were extracted a small buffer is
applied to try and obtain the remaining values.

## Usage

``` r
get_rastervalue(coords, env, ngb_fill = TRUE, rm.na = FALSE)
```

## Arguments

- coords:

  A [`data.frame`](https://rdrr.io/r/base/data.frame.html),
  [`matrix`](https://rdrr.io/r/base/matrix.html) or
  [`sf::sf`](https://r-spatial.github.io/sf/reference/sf.html) object.

- env:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object with the provided predictors.

- ngb_fill:

  [`logical`](https://rdrr.io/r/base/logical.html) on whether cells
  should be interpolated from neighbouring values.

- rm.na:

  [`logical`](https://rdrr.io/r/base/logical.html) parameter which - if
  set - removes all rows with a missing data point (`NA`) from the
  result.

## Value

A [`data.frame`](https://rdrr.io/r/base/data.frame.html) with the
extracted covariate data from each provided data point.

## Details

It is essentially a wrapper for
[`terra::extract`](https://rspatial.github.io/terra/reference/extract.html).

## Examples

``` r
# Dummy raster:
r <- terra::rast(nrows = 10, ncols = 10, res = 0.05, xmin = -1.5, xmax = 1.5,
ymin = -1.5, ymax = 1.5, vals = rnorm(3600,mean = .5,sd = .1))
# (dummy points)
pp <- terra::spatSample(r,20,as.points = TRUE) |> sf::st_as_sf()

# Extract values
vals <- get_rastervalue(pp, r)
head(vals)
#>   ID     lyr.1      x      y
#> 1  1 0.6782599  1.225 -1.175
#> 2  2 0.6420574 -0.225  0.925
#> 3  3 0.5397701  0.075 -1.425
#> 4  4 0.4993814  0.975  0.775
#> 5  5 0.5545246  0.175  1.275
#> 6  6 0.4226162  0.625  0.625
```
