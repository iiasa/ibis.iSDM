# Add a control to a BiodiversityModel object to control biases

Sampling and other biases are pervasive drivers of the spatial location
of biodiversity datasets. While the integration of other, presumably
less biased data can be one way of controlling for sampling biases,
another way is to control directly for the bias in the model. Currently
supported methods are:

- `"partial"` - An approach described by Warton et al. (2013) to control
  the biases in a model, by including a specified variable ("layer") in
  the model, but "partialling" it out during the projection phase.
  Specifically the variable is set to a specified value ("bias_value"),
  which is by default the minimum value observed across the background.

- `"offset"` - Dummy method that points to the
  [`add_offset_bias()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset_bias.md)
  functionality (see note). Makes use of offsets to factor out a
  specified bias variable.

- `"proximity"` - Use the proximity or distance between points as a
  weight in the model. This option effectively places greater weight on
  points farther away. *Note:* In the best case this can control for
  spatial bias and aggregation, in the worst case it can place a lot of
  emphasis on points that likely outliers or misidentification (in terms
  of species).

See also details for some explanations.

## Usage

``` r
add_control_bias(
  x,
  layer,
  method = "partial",
  bias_value = NULL,
  maxdist = NULL,
  alpha = 1,
  add = TRUE
)

# S4 method for class 'BiodiversityDistribution'
add_control_bias(
  x,
  layer,
  method = "partial",
  bias_value = NULL,
  maxdist = NULL,
  alpha = 1,
  add = TRUE
)
```

## Arguments

- x:

  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

- layer:

  A [`sf::sf`](https://r-spatial.github.io/sf/reference/sf.html) or
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object with the range for the target feature. Specify a variable that
  is not already added to `"x"` to avoid issues with duplications.

- method:

  A [`character`](https://rdrr.io/r/base/character.html) vector
  describing the method used for bias control. Available options are
  `"partial"` (Default), `"offset"` or `"proximity"`.

- bias_value:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) with a value for
  `"layer"`. Specifying a
  [`numeric`](https://rdrr.io/r/base/numeric.html) value here sets
  `layer` to the target value during projection. By default the value is
  set to the minimum value found in the layer (Default: `NULL`).

- maxdist:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) giving the maximum
  distance if method `"proximity"` is used. If unset it uses by default
  the distance to the centroid of a minimum convex polygon encircling
  all points.

- alpha:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) given the initial
  weight to points if method `"proximity"` is used (Default: `1`). For
  example, if set to values smaller than `1` neighbouring points will be
  weighted less.

- add:

  [`logical`](https://rdrr.io/r/base/logical.html) specifying whether a
  new offset is to be added. Setting this parameter to `FALSE` replaces
  the current offsets with the new one (Default: `TRUE`).

## Value

Adds bias control option to a
[`distribution`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
object.

## Details

In the case of `"proximity"` weights are assigned to each point, placing
higher weight on points further away and with less overlap. Weights are
assigned up to a maximum of distance which can be provided by the user
(parameter `"maxdist"`). This distance is ideally informed by some
knowledge of the species to be modelled (e.g., maximum dispersal
distance). If not provided, it is set to the distance of the centroid of
a minimum convex polygon encircling all observations. The parameter
`"alpha"` is a weighting factor which can be used to diminish the effect
of neighboring points.  
For a given observation \\i\\, the weight \\w\\ is defined as \$\$w_i =
1 / (1 + \epsilon)\$\$ where \$\$\epsilon = \sum\_{n=1}^{N}((1 -
d_n)/d_sac)^\alpha\$\$ in which \\N\\ is the total number of points
closer than the maximum distance (\\d_sac\\) of point \\i\\, and \\d_n\\
the distance between focal point \\i\\ and point \\n\\.

## Note

**Covariate transformations applied to other predictors need to be
applied to bias too.** Another option to consider biases particular in
Poisson-point process models is to remove them through an offset.
Functionality to do so is available through the
[`add_offset_bias()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset_bias.md)
method. Setting the method to `"offset"` will automatically point to
this option.

## References

- Warton, D.I., Renner, I.W. and Ramp, D., 2013. Model-based control of
  observer bias for the analysis of presence-only data in ecology. PloS
  one, 8(11), p.e79168.

- Merow, C., Allen, J.M., Aiello-Lammens, M., Silander, J.A., 2016.
  Improving niche and range estimates with Maxent and point process
  models by integrating spatially explicit information. Glob. Ecol.
  Biogeogr. 25, 1022–1036. https://doi.org/10.1111/geb.12453

- Botella, C., Joly, A., Bonnet, P., Munoz, F., & Monestiez, P. (2021).
  Jointly estimating spatial sampling effort and habitat suitability for
  multiple species from opportunistic presence‐only data. Methods in
  Ecology and Evolution, 12(5), 933-945.

## See also

[`add_limits_extrapolation()`](https://iiasa.github.io/ibis.iSDM/reference/add_limits_extrapolation.md)

## Examples

``` r
if (FALSE) { # \dontrun{
 x <- distribution(background) |>
   add_predictors(covariates) |>
   add_control_bias(biasvariable, bias_value = NULL)
} # }
```
