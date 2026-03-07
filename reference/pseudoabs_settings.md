# Settings for specifying pseudo-absence points within the model background

This function defines the settings for pseudo-absence sampling of the
background. For many engines such points are necessary to model Poisson
(or Binomial) distributed point process data. Specifically we call
absence points for Binomial (Bernoulli really) distributed responses
`'pseudo-absence'` and absence data for Poisson responses `'background'`
points. For more details read Renner et al. (2015).

The function `'add_pseudoabsence'` allows to add absence points to any
[`sf::sf`](https://r-spatial.github.io/sf/reference/sf.html) object. See
**Details** for additional parameter description and examples on how to
'turn' a presence-only dataset into a presence-(pseudo-)absence.

## Usage

``` r
pseudoabs_settings(
  background = NULL,
  nrpoints = 10000,
  min_ratio = 0.25,
  method = "random",
  buffer_distance = 10000,
  inside = FALSE,
  layer = NULL,
  bias = NULL,
  ...
)

# S4 method for class 'ANY'
pseudoabs_settings(
  background = NULL,
  nrpoints = 10000,
  min_ratio = 0.25,
  method = "random",
  buffer_distance = 10000,
  inside = FALSE,
  layer = NULL,
  bias = NULL,
  ...
)
```

## Arguments

- background:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  or [`sf::sf`](https://r-spatial.github.io/sf/reference/sf.html) object
  over which background points can be sampled. Default is `NULL`
  (Default) and the background is then added when the sampling is first
  called.

- nrpoints:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) given the number of
  absence points to be created. Has to be larger than `0` and normally
  points are not created in excess of the number of cells of the
  background (Default: `10000`).

- min_ratio:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) with the minimum
  ratio of background points relative to the presence points. Setting
  this value to `1` generates an equal amount of absence points relative
  to the presence points. Usually ignored unless the ratio exceeds the
  `nrpoints` parameters (Default: `0.25`).

- method:

  [`character`](https://rdrr.io/r/base/character.html) denoting how the
  sampling should be done. See details for options (Default:
  `"random"`).

- buffer_distance:

  [`numeric`](https://rdrr.io/r/base/numeric.html) A distance from the
  observations in which pseudo-absence points are not to be generated.
  Note that units follow the units of the projection (e.g. `m` or `°`).
  Only used when `method = "buffer"`.

- inside:

  A [`logical`](https://rdrr.io/r/base/logical.html) value of whether
  absence points should be sampled outside (Default) or inside a minimum
  convex polygon or range provided the respective method is chosen
  (parameter `method = "mcp"` or `method = "range"`).

- layer:

  A [`sf::sf`](https://r-spatial.github.io/sf/reference/sf.html) or
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  (in the case of method `'zones'`) object indicating the range of a
  species. Only used with `method = "range"` or `method = "zones"`
  (Default: `NULL`).

- bias:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  with the same extent and projection and background. Absence points
  will be preferentially sampled in areas with higher (!) bias.
  (Default: `NULL`).

- ...:

  Any other settings to be added to the pseudoabs settings.

## Details

There are multiple methods available for sampling a biased background
layer. Possible parameters for `method` are:

- `'random'` Absence points are generated randomly over the background
  (Default),

- `'buffer'` Absence points are generated only within a buffered
  distance of existing points. This option requires the specification of
  the parameter `buffer_distance`.

- `'mcp'` Can be used to only generate absence points within or outside
  a minimum convex polygon of the presence points. The parameter
  `inside` specifies whether points should be sampled inside or outside
  (Default) the minimum convex polygon.

- `'range'` Absence points are created either inside or outside a
  provided additional layer that indicates for example a range of
  species (controlled through parameter `inside`).

- `'zones'` A ratified (e.g. of type
  [factor](https://rdrr.io/r/base/factor.html))
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  layer depicting zones from which absence points are to be sampled.
  This method checks which points fall within which zones and then
  samples absence points either within or outside these zones
  exclusively. Both `'layer'` and `'inside'` have to be set for this
  option.

- `'target'` Make use of a target background for sampling absence
  points. Here a
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object has to be provided through the parameter `'layer'`. Absence
  points are then sampled exclusively within the target areas for grid
  cells with non-zero values.

## References

- Renner IW, Elith J, Baddeley A, Fithian W, Hastie T, Phillips SJ,
  Popovic G, Warton DI. 2015. Point process models for presence-only
  analysis. Methods in Ecology and Evolution 6:366–379. DOI:
  10.1111/2041-210X.12352.

- Renner, I. W., & Warton, D. I. (2013). Equivalence of MAXENT and
  Poisson point process models for species distribution modeling in
  ecology. Biometrics, 69(1), 274-281.

## Examples

``` r
if (FALSE) { # \dontrun{
# This setting generates 10000 pseudo-absence points outside the
# minimum convex polygon of presence points
ass1 <- pseudoabs_settings(nrpoints = 10000, method = 'mcp', inside = FALSE)

# This setting would match the number of presence-absence points directly.
ass2 <- pseudoabs_settings(nrpoints = 0, min_ratio = 1)

# These settings can then be used to add pseudo-absence data to a
# presence-only dataset. This effectively adds these simulated absence
# points to the resulting model
all_my_points <- add_pseudoabsence(
                     df = virtual_points,
                      field_occurrence = 'observed',
                      template = background,
                      settings = ass1)
} # }
```
