# Add a control to a BiodiversityModel object to limit extrapolation

One of the main aims of species distribution models (SDMs) is to project
in space and time. For projections a common issue is extrapolation as -
unconstrained - SDMs can indicate areas as suitable which are unlikely
to be occupied by species or habitats (often due to historic or biotic
factors). To some extent this can be related to an insufficient
quantification of the niche (e.g. niche truncation by considering only a
subset of observations within the actual distribution), in other cases
there can also be general barriers or constraints that limit any
projections (e.g. islands). This limit method adds some of those options
to a model distribution object. Currently supported methods are:

- `"zones"` - This is a wrapper to allow the addition of zones to a
  distribution model object, similar to what is also possible via
  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md).
  Required is a spatial layer that describes an environmental zoning.

- `"mcp"` - Rather than using an external or additional layer, this
  option constrains predictions by a certain distance of points in its
  vicinity. Buffer distances have to be in the unit of the projection
  used and can be configured via `"mcp_buffer"`.

- `"nt2"` - Constrains the predictions using the multivariate
  combination novelty index (NT2) following Mesgaran et al. (2014). This
  method is also available in the
  [`similarity()`](https://iiasa.github.io/ibis.iSDM/reference/similarity.md)
  function.

- `"mess"` - Constrains the predictions using the Multivariate
  Environmental Similarity Surfaces (MESS) following Mesgaran et al.
  (2014). This method is also available in the
  [`similarity()`](https://iiasa.github.io/ibis.iSDM/reference/similarity.md)
  function.

- `"shape"` - This is an implementation of the 'shape' method introduced
  by Velazco et al. (2023). Through a user defined threshold it
  effectively limits model extrapolation so that no projections are made
  beyond the extent judged as defensible and informed by the training
  observations. **Not yet implemented!**

See also details for further explanations.

## Usage

``` r
add_limits_extrapolation(
  x,
  layer,
  method = "mcp",
  mcp_buffer = 0,
  novel = "within",
  limits_clip = FALSE
)

# S4 method for class 'BiodiversityDistribution'
add_limits_extrapolation(
  x,
  layer,
  method = "mcp",
  mcp_buffer = 0,
  novel = "within",
  limits_clip = FALSE
)
```

## Arguments

- x:

  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

- layer:

  A
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  or [`sf`](https://r-spatial.github.io/sf/reference/sf.html) object
  that limits the prediction surface when intersected with input data
  (Default: `NULL`).

- method:

  A [`character`](https://rdrr.io/r/base/character.html) vector
  describing the method used for controlling extrapolation. Available
  options are `"zones"`, `"mcp"` (Default), or `"nt2"`, `"mess"` or
  `"shape"`.

- mcp_buffer:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) distance to buffer
  the mcp (Default `0`). Only used if `"mcp"` is used.

- novel:

  Which conditions are to be masked out respectively, either the novel
  conditions within only `"within"` (Default) or also including outside
  reference conditions `"outside"`. Only use for `method = "nt2"`, for
  `method = "mess"` this variable is always `"within"`.

- limits_clip:

  [`logical`](https://rdrr.io/r/base/logical.html) Should the limits
  clip all predictors before fitting a model (`TRUE`) or just the
  prediction (`FALSE`, default).

## Value

Adds extrapolation limit option to a
[`distribution`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
object.

## Details

For method `"zones"` a zoning layer can be supplied which is then used
to intersect the provided training points with. Any projections made
with the model can then be constrained so as to not project into areas
that do not consider any training points and are unlikely to have any.
Examples for zones are for the separation of islands and mainlands,
biomes, or lithological soil conditions.

If no layer is available, it is also possible to constraint predictions
by the distance to a minimum convex polygon surrounding the training
points with method `"mcp"` (optionally buffered). This can make sense
particular for rare species or those fully sampled across their niche.

For the `"NT2"` and `"MESS"` index it is possible to constrain the
prediction to conditions within (`novel = "within"`) or also include
outside (`novel = "outside"`) conditions.

## Note

The method `"zones"` is also possible directly within
[`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md).

## References

- Randin, C. F., Dirnböck, T., Dullinger, S., Zimmermann, N. E., Zappa,
  M., & Guisan, A. (2006). Are niche‐based species distribution models
  transferable in space?. Journal of biogeography, 33(10), 1689-1703.
  https://doi.org/10.1111/j.1365-2699.2006.01466.x

- Chevalier, M., Broennimann, O., Cornuault, J., & Guisan, A. (2021).
  Data integration methods to account for spatial niche truncation
  effects in regional projections of species distribution. Ecological
  Applications, 31(7), e02427. https://doi.org/10.1002/eap.2427

- Velazco, S. J. E., Brooke, M. R., De Marco Jr., P., Regan, H. M., &
  Franklin, J. (2023). How far can I extrapolate my species distribution
  model? Exploring Shape, a novel method. Ecography, 11, e06992.
  https://doi.org/10.1111/ecog.06992

- Mesgaran, M. B., R. D. Cousens, B. L. Webber, and J. Franklin. (2014)
  Here be dragons: a tool for quantifying novelty due to covariate range
  and correlation change when projecting species distribution models.
  Diversity and Distributions 20:1147-1159.

## Examples

``` r
if (FALSE) { # \dontrun{
 # To add a zone layer for extrapolation constraints.
 x <- distribution(background) |>
   add_predictors(covariates) |>
   add_limits_extrapolation(method = "zones", layer = zones)
} # }
```
