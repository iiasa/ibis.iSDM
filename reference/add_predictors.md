# Add predictors to a Biodiversity distribution object

This function allows to add predictors to
[distribution](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
or
[BiodiversityScenario](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
objects. Predictors are covariates that in spatial projection have to
match the geographic projection of the background layer in the
[distribution](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
object. This function furthermore allows to transform or create
derivates of provided predictors.

## Usage

``` r
add_predictors(
  x,
  env,
  names = NULL,
  transform = "none",
  derivates = "none",
  derivate_knots = 4,
  int_variables = NULL,
  bgmask = TRUE,
  harmonize_na = FALSE,
  explode_factors = FALSE,
  priors = NULL,
  state = NULL,
  ...
)

# S4 method for class 'BiodiversityDistribution,SpatRasterCollection'
add_predictors(
  x,
  env,
  names = NULL,
  transform = "none",
  derivates = "none",
  derivate_knots = 4,
  int_variables = NULL,
  bgmask = TRUE,
  harmonize_na = FALSE,
  explode_factors = FALSE,
  priors = NULL,
  state = NULL,
  ...
)

# S4 method for class 'BiodiversityDistribution,SpatRaster'
add_predictors(
  x,
  env,
  names = NULL,
  transform = "none",
  derivates = "none",
  derivate_knots = 4,
  int_variables = NULL,
  bgmask = TRUE,
  harmonize_na = FALSE,
  explode_factors = FALSE,
  priors = NULL,
  state = NULL,
  ...
)

# S4 method for class 'BiodiversityDistribution,stars'
add_predictors(
  x,
  env,
  names = NULL,
  transform = "none",
  derivates = "none",
  derivate_knots = 4,
  int_variables = NULL,
  bgmask = TRUE,
  harmonize_na = FALSE,
  explode_factors = FALSE,
  priors = NULL,
  state = NULL,
  ...
)

# S4 method for class 'BiodiversityDistribution,data.frame'
add_predictors(
  x,
  env,
  names = NULL,
  transform = "none",
  derivates = "none",
  derivate_knots = 4,
  int_variables = NULL,
  bgmask = TRUE,
  harmonize_na = FALSE,
  explode_factors = FALSE,
  priors = NULL,
  state = NULL,
  ...
)

# S4 method for class 'BiodiversityScenario,SpatRaster'
add_predictors(
  x,
  env,
  names = NULL,
  transform = "none",
  derivates = "none",
  derivate_knots = 4,
  int_variables = NULL,
  bgmask = TRUE,
  harmonize_na = FALSE,
  explode_factors = FALSE,
  priors = NULL,
  state = NULL,
  ...
)

# S4 method for class 'BiodiversityScenario,stars'
add_predictors(
  x,
  env,
  names = NULL,
  transform = "none",
  derivates = "none",
  derivate_knots = 4,
  int_variables = NULL,
  bgmask = TRUE,
  harmonize_na = FALSE,
  explode_factors = FALSE,
  priors = NULL,
  state = NULL,
  ...
)
```

## Arguments

- x:

  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

- env:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html),
  [`stars`](https://rdrr.io/r/graphics/stars.html) or
  [`data.frame`](https://rdrr.io/r/base/data.frame.html) object.

- names:

  A [`vector`](https://rdrr.io/r/base/vector.html) of character names
  describing the environmental stack in case they should be renamed.

- transform:

  A [`vector`](https://rdrr.io/r/base/vector.html) stating whether
  predictors should be preprocessed in any way (Options:
  `'none'`,`'pca'`, `'scale'`, `'norm'`)

- derivates:

  A Boolean check whether derivate features should be considered
  (Options: `'none'`, `'thresh'`, `'hinge'`, `'quad'`) )

- derivate_knots:

  A single [`numeric`](https://rdrr.io/r/base/numeric.html) or
  [`vector`](https://rdrr.io/r/base/vector.html) giving the number of
  knots for derivate creation if relevant (Default: `4`).

- int_variables:

  A [`vector`](https://rdrr.io/r/base/vector.html) with length greater
  or equal than `2` specifying the covariates (Default: `NULL`).

- bgmask:

  Check whether the environmental data should be masked with the
  background layer (Default: `TRUE`).

- harmonize_na:

  A [`logical`](https://rdrr.io/r/base/logical.html) value indicating of
  whether NA values should be harmonized among predictors (Default:
  `FALSE`).

- explode_factors:

  [`logical`](https://rdrr.io/r/base/logical.html) of whether any factor
  variables should be split up into binary variables (one per class).
  (Default: `FALSE`).

- priors:

  A
  [`PriorList`](https://iiasa.github.io/ibis.iSDM/reference/PriorList-class.md)
  object. Default is set to `NULL` which uses default prior assumptions.

- state:

  A [`matrix`](https://rdrr.io/r/base/matrix.html) with one value per
  variable (column) providing either a ( `stats::mean()`,
  [`stats::sd()`](https://rdrr.io/r/stats/sd.html) ) for each variable
  in `env` for option `'scale'` or a range of minimum and maximum values
  for option `'norm'`. Effectively applies their value range for
  rescaling. In the case of provided `stars` data to a
  [BiodiversityScenario](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
  object, the state variables are attempted to be compiled from the
  predictor ranges used for model inferrence (Default: `NULL`).

- ...:

  Other parameters passed down

## Details

A transformation takes the provided rasters and for instance rescales
them or transforms them through a principal component analysis
([prcomp](https://rdrr.io/r/stats/prcomp.html)). In contrast, derivates
leave the original provided predictors alone, but instead create new
ones, for instance by transforming their values through a quadratic or
hinge transformation. Note that this effectively increases the number of
predictors in the object, generally requiring stronger regularization by
the used
[`Engine`](https://iiasa.github.io/ibis.iSDM/reference/Engine-class.md).
Both transformations and derivates can also be combined. Available
options for transformation are:

- `'none'` - Leaves the provided predictors in the original scale.

- `'pca'` - Converts the predictors to principal components. Note that
  this results in a renaming of the variables to principal component
  axes!

- `'scale'` - Transforms all predictors by applying
  [scale](https://rdrr.io/r/base/scale.html) on them.

- `'norm'` - Normalizes all predictors by transforming them to a scale
  from `0` to `1`.

- `'windsor'` - Applies a windsorization to the target predictors. By
  default this effectively cuts the predictors to the `0.05` and `0.95`,
  thus helping to remove extreme outliers.

Available options for creating derivates are:

- `'none'` - No additional predictor derivates are created.

- `'quad'` - Adds quadratic derivate predictors.

- `'interaction'` - Add interacting predictors. Interactions need to be
  specified (`"int_variables"`)!

- `'thresh'` - Add threshold derivate predictors.

- `'hinge'` - Add hinge derivate predictors.

- `'kmeans'` - Add k-means derived factors.

- `'bin'` - Add predictors binned by their percentiles.

## Note

**Important:** Not every
[`Engine`](https://iiasa.github.io/ibis.iSDM/reference/Engine-class.md)
supported by the ibis.iSDM R-package allows missing data points among
extracted covariates. Thus any observation with missing data is
generally removed prior from model fitting. Thus ensure that covariates
have appropriate no-data settings (for instance setting `NA` values to
`0` or another out of range constant).

Not every engine does actually need covariates. For instance it is
perfectly legit to fit a model with only occurrence data and a spatial
latent effect
([add_latent_spatial](https://iiasa.github.io/ibis.iSDM/reference/add_latent_spatial.md)).
This correspondents to a spatial kernel density estimate.

Certain names such `"offset"` are forbidden as predictor variable names.
The function will return an error message if these are used.

Some engines use binary variables regardless of the parameter
`"explode_factors"` set here.

## Examples

``` r
if (FALSE) { # \dontrun{
 obj <- distribution(background) |>
        add_predictors(covariates, transform = 'scale')
 obj
} # }
```
