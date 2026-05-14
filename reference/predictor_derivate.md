# Create spatial derivative of raster stacks

This function creates derivatives of existing covariates and returns
them in Raster format. Derivative variables can in the machine learning
literature commonly be understood as one aspect of feature engineering.
They can be particularly powerful in introducing non-linearities in
otherwise linear models, for example is often done in the popular Maxent
framework.

## Usage

``` r
predictor_derivate(
  env,
  option,
  nknots = 4,
  deriv = NULL,
  int_variables = NULL,
  method = NULL,
  ...
)
```

## Arguments

- env:

  A
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object.

- option:

  A [`vector`](https://rdrr.io/r/base/vector.html) stating whether
  predictors should be preprocessed in any way (Options: `'none'`,
  `'quadratic'`, `'hinge'`, `'kmeans'`, `'thresh'`, `'bin'`).

- nknots:

  The number of knots to be used for the transformation (Default: `4`).

- deriv:

  A [`vector`](https://rdrr.io/r/base/vector.html) with
  [`character`](https://rdrr.io/r/base/character.html) of specific
  derivates to create (Default: `NULL`).

- int_variables:

  A [`vector`](https://rdrr.io/r/base/vector.html) with length greater
  or equal than `2` specifying the covariates (Default: `NULL`).

- method:

  As `'option'` for more intuitive method setting. Can be left empty (in
  this case option has to be set).

- ...:

  other options (Non specified).

## Value

Returns the derived adjusted
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
objects of identical resolution.

## Details

Available options are:

- `'none'` - The original layer(s) are returned.

- `'quadratic'` - A quadratic transformation (\\x^{2}\\) is created of
  the provided layers.

- `'hinge'` - Creates hinge transformation of covariates, which set all
  values lower than a set threshold to `0` and all others to a range of
  \\\[0,1\]\\. The number of thresholds and thus new derivates is
  specified via the parameter `'nknots'` (Default: `4`).

- `'interaction'` - Creates interactions between variables. Target
  variables have to be specified via `"int_variables"`.

- `'thresh'` - A threshold transformation of covariates, which sets all
  values lower than a set threshold at `0` and those larger to `1`. The
  number of thresholds and thus new derivates is specified via the
  parameter `'nknots'` (Default: `4`).

- `'bin'` - Creates a factor representation of a covariates by cutting
  the range of covariates by their percentiles. The number of percentile
  cuts and thus new derivates is specified via the parameter `'nknots'`
  (Default: `4`).

- `'kmeans'` Creates a factor representation of a covariates through a
  [`kmeans()`](https://rdrr.io/r/stats/kmeans.html) clustering. The
  number of clusters are specified via the parameter `'nknots'`.

## See also

predictor_transform

## Examples

``` r
# Dummy raster
r_ori <- terra::rast(nrows = 10, ncols = 10, res = 0.05, xmin = -1.5, xmax = 1.5, ymin = -1.5, ymax = 1.5, vals = rpois(3600, 10))

# Create a hinge transformation with 4 knots of one or multiple SpatRaster.
new <- predictor_derivate(r_ori, option = "hinge", knots = 4)
terra::plot(new)


# Or a quadratic transformation
new2 <- predictor_derivate(r_ori, option = "quad", knots = 4)
terra::plot(new2)

```
