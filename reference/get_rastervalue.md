# Function to extract point values directly from a SpatRaster

This function simply extracts the values from a provided
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html),
[`SpatRasterDataset`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
or
[`SpatRasterCollection`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
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
  [`sf`](https://r-spatial.github.io/sf/reference/sf.html) object.

- env:

  A
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
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
#> 1  1 0.5154511  1.475 -1.275
#> 2  2 0.3948053  1.275  1.175
#> 3  3 0.5363174 -1.025  0.225
#> 4  4 0.3373516  0.625 -1.275
#> 5  5 0.5024094  0.225 -0.775
#> 6  6 0.4530912  1.025  1.275
```
