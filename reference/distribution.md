# Create distribution modelling procedure

This function creates an object that contains all the data, parameters
and settings for building an (integrated) species distribution model.
Key functions to add data are
[`add_biodiversity_poipo`](https://iiasa.github.io/ibis.iSDM/reference/add_biodiversity_poipo.md)
and the like,
[`add_predictors`](https://iiasa.github.io/ibis.iSDM/reference/add_predictors.md),
[`add_latent_spatial`](https://iiasa.github.io/ibis.iSDM/reference/add_latent_spatial.md),
[`engine_glmnet`](https://iiasa.github.io/ibis.iSDM/reference/engine_glmnet.md)
or similar,
[`add_priors`](https://iiasa.github.io/ibis.iSDM/reference/add_priors.md)
and
[`add_offset`](https://iiasa.github.io/ibis.iSDM/reference/add_offset.md).
It creates a prototype
[`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md)
object with its own functions. After setting input data and parameters,
model predictions can then be created via the
[train](https://iiasa.github.io/ibis.iSDM/reference/train.md) function
and predictions be created.

Additionally, it is possible to specify a `"limit"` to any predictions
conducted on the background. This can be for instance a buffered layer
by a certain dispersal distance (Cooper and Soberon, 2018) or a
categorical layer representing biomes or soil conditions. Another option
is to create a constraint by constructing a minimum convex polygon (MCP)
using the supplied biodiversity data. This option can be enabled by
setting `"limits_method"` to `"mcp"`. It is also possible to provide a
small buffer to constructed MCP that way. See the frequently asked
question (FAQ) section on the homepage for more information.

See **Details** for a description of the internal functions available to
modify or summarize data within the created object.

**Note that any model requires at minimum a single added biodiversity
dataset as well as a specified engine.**

## Usage

``` r
distribution(
  background,
  limits = NULL,
  limits_method = "none",
  mcp_buffer = 0,
  limits_clip = FALSE
)

# S4 method for class 'SpatRaster'
distribution(
  background,
  limits = NULL,
  limits_method = "none",
  mcp_buffer = 0,
  limits_clip = FALSE
)

# S4 method for class 'sf'
distribution(
  background,
  limits = NULL,
  limits_method = "none",
  mcp_buffer = 0,
  limits_clip = FALSE
)
```

## Arguments

- background:

  Specification of the modelling background. Must be a
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  or [`sf::sf`](https://r-spatial.github.io/sf/reference/sf.html)
  object.

- limits:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html),
  [`sf::sf`](https://r-spatial.github.io/sf/reference/sf.html) or
  [`stars`](https://rdrr.io/r/graphics/stars.html) object that limits
  the prediction surface when intersected with input data (Default:
  `NULL`). In case of a [`stars`](https://rdrr.io/r/graphics/stars.html)
  object the first factorized time entry is taken.

- limits_method:

  A [`character`](https://rdrr.io/r/base/character.html) of the method
  used for hard limiting a projection. Available options are `"none"`
  (Default), `"zones"` or `"mcp"`. See also
  [`add_limits_extrapolation()`](https://iiasa.github.io/ibis.iSDM/reference/add_limits_extrapolation.md).

- mcp_buffer:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) distance to buffer
  the mcp (Default `0`). Only used if `"mcp"` is used.

- limits_clip:

  [`logical`](https://rdrr.io/r/base/logical.html) Should the limits
  clip all predictors before fitting a model (`TRUE`) or just the
  prediction (`FALSE`, default).

## Value

[`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md)
object containing data for building a biodiversity distribution
modelling problem.

## Details

This function creates a
[`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md)
object that in itself contains other functions and stores parameters and
(pre-)processed data. A full list of functions available can be queried
via `"names(object)"`. Some of the functions are not intended to be
manipulated directly, but rather through convenience functions (e.g.
`"object$set_predictors()"`). Similarly other objects are stored in the
[`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md)
object that have their own functions as well and can be queried (e.g.
`"names(object)"`). For a list of functions see the reference
documentation. By default, if some datasets are not set, then a
`"Waiver"` object is returned instead.

The following objects can be stored:

- `object$biodiversity` A
  [`BiodiversityDatasetCollection`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDatasetCollection-class.md)
  object with the added biodiversity data.

- `object$engine` An `"engine"` object (e.g.
  [`engine_inlabru()`](https://iiasa.github.io/ibis.iSDM/reference/engine_inlabru.md))
  with function depended on the added engine.

- `object$predictors` A
  [`PredictorDataset`](https://iiasa.github.io/ibis.iSDM/reference/PredictorDataset-class.md)
  object with all set predictors.

- `object$priors` A
  [`PriorList`](https://iiasa.github.io/ibis.iSDM/reference/PriorList-class.md)
  object with all specified priors.

- `object$log` A
  [`Log`](https://iiasa.github.io/ibis.iSDM/reference/Log-class.md)
  object that captures messages.

Useful high-level functions to address those objects are for instance:

- `object$show()` A generic summary of the
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md)
  object contents. Can also be called via
  [print](https://iiasa.github.io/ibis.iSDM/reference/print.md).

- `object$get_biodiversity_equations()` Lists the equations used for
  each biodiversity dataset with given id. Defaults to all predictors.

- `object$get_biodiversity_types()` Lists the type of each specified
  biodiversity dataset with given id.

- `object$get_extent()` Outputs the
  [terra::ext](https://rspatial.github.io/terra/reference/ext.html) of
  the modelling region.

- `object$show_background_info()` Returns a
  [`list`](https://rdrr.io/r/base/list.html) with the
  [terra::ext](https://rspatial.github.io/terra/reference/ext.html) and
  the [terra::crs](https://rspatial.github.io/terra/reference/crs.html).

- `object$get_extent_dimensions()` Outputs the
  [terra::ext](https://rspatial.github.io/terra/reference/ext.html)
  dimension by calling the `"extent_dimensions()"` function.

- `object$get_predictor_names()` Returns a
  [character](https://rdrr.io/r/base/character.html) vector with the
  names of all added predictors.

- `object$get_prior_variables()` Returns a description of
  [`priors`](https://iiasa.github.io/ibis.iSDM/reference/priors.md)
  added.

There are other functions as well but those are better accessed through
their respective wrapper functions.

## References

- Fletcher, R.J., Hefley, T.J., Robertson, E.P., Zuckerberg, B.,
  McCleery, R.A., Dorazio, R.M., (2019) A practical guide for combining
  data to model species distributions. Ecology 100, e02710.
  https://doi.org/10.1002/ecy.2710

- Cooper, Jacob C., and Jorge Soberón. "Creating individual accessible
  area hypotheses improves stacked species distribution model
  performance." Global Ecology and Biogeography 27, no. 1 (2018):
  156-165.

## See also

[`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md)
and other classes.

## Examples

``` r
# Load background raster
background <- terra::rast(system.file("extdata/europegrid_50km.tif",package = "ibis.iSDM"))
# Define model
x <- distribution(background)
#> [Setup] 2026-03-07 11:56:08.360655 | Creating distribution object...
x
#> <Biodiversity distribution model>
#> Background extent: 
#>      xmin: -16.064, xmax: 34.95,
#>      ymin: 36.322, ymax: 71.535
#>    projection: +proj=longlat +datum=WGS84 +no_defs
#>  --------- 
#> Biodiversity data:
#>    None
#>  --------- 
#>   predictors:     None
#>   priors:         <Default>
#>   latent:         None
#>   log:            <Console>
#>   engine:         <NONE>
```
