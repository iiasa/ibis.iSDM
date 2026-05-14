# Spatial adjustment of environmental predictors and raster stacks

This function allows the transformation of provided environmental
predictors (in
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
format). A common use case is for instance the standardization (or
scaling) of all predictors prior to model fitting. This function works
both with
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
as well as with [`stars`](https://rdrr.io/r/graphics/stars.html)
objects.

## Usage

``` r
predictor_transform(
  env,
  option,
  windsor_props = c(0.05, 0.95),
  pca.var = 0.8,
  state = NULL,
  method = NULL,
  ...
)
```

## Arguments

- env:

  A
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  or [`stars`](https://rdrr.io/r/graphics/stars.html) object.

- option:

  A [`vector`](https://rdrr.io/r/base/vector.html) stating whether
  predictors should be preprocessed in any way (Options: `'none'`,
  `'scale'`, `'norm'`, `'windsor'`, `'windsor_thresh'`, `'percentile'`
  `'pca'`, `'revjack'`). See Details.

- windsor_props:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) vector specifying
  the proportions to be clipped for windsorization (Default:
  `c(.05,.95)`).

- pca.var:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) value between `>0`
  and `1` stating the minimum amount of variance to be covered (Default:
  `0.8`).

- state:

  A [`matrix`](https://rdrr.io/r/base/matrix.html) with one value per
  variable (column) providing either a ( `stats::mean()`,
  [`stats::sd()`](https://rdrr.io/r/stats/sd.html) ) for each variable
  in `env` for option `'scale'` or a range of minimum and maximum values
  for option `'norm'`. Effectively applies their value range for
  rescaling. (Default: `NULL`).

- method:

  As `'option'` for more intuitive method setting. Can be left empty (in
  this case option has to be set).

- ...:

  Currrently not implemented (Non specified).

## Value

Returns a adjusted
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
object of identical resolution.

## Details

Available options are:

- `'none'` The original layer(s) are returned.

- `'scale'` This run the
  [`scale()`](https://rspatial.github.io/terra/reference/scale.html)
  function with default settings (1 Standard deviation) across all
  predictors. A sensible default to for most model fitting.

- `'norm'` This normalizes all predictors to a range from `0-1`.

- `'windsor'` This applies a 'windsorization' to an existing raster
  layer by setting the lowest, respectively largest values to the value
  at a certain percentage level (e.g. 95%). Those can be set via the
  parameter `"windsor_props"`.

- `'windsor_thresh'` Same as option 'windsor', however in this case
  values are clamped to a thresholds rather than certain percentages
  calculated on the data.

- `'percentile'` This converts and bins all values into percentiles,
  e.g. the top 10% or lowest 10% of values and so on.

- `'pca'` This option runs a principal component decomposition of all
  predictors (via
  [`prcomp()`](https://rspatial.github.io/terra/reference/prcomp.html)).
  It returns new predictors resembling all components in order of the
  most important ones. Can be useful to reduce collinearity, however
  note that this changes all predictor names to 'PCX', where X is the
  number of the component. The parameter `'pca.var'` can be modified to
  specify the minimum variance to be covered by the axes.

- `'revjack'` Removes outliers from the supplied stack via a reverse
  jackknife procedure. Identified outliers are by default set to `NA`.

## Note

If future covariates are rescaled or normalized, it is highly
recommended to use the statistical moments on which the models were
trained for any variable transformations, also to ensure that variable
ranges are consistent among relative values.

## See also

predictor_derivate

## Examples

``` r
# Dummy raster
r_ori <- terra::rast(nrows = 10, ncols = 10, res = 0.05, xmin = -1.5, xmax = 1.5, ymin = -1.5, ymax = 1.5, vals = rnorm(3600,mean = .01,sd = .1))

# Normalize
r_norm <- predictor_transform(r_ori, option = 'norm')
new <- c(r_ori, r_norm)
names(new) <- c("original scale", "normalized units")
terra::plot(new)

```
