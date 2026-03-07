# Add biodiversity polygon dataset to a distribution object (presence-absence)

This function can be used to add a
[`sf::sf`](https://r-spatial.github.io/sf/reference/sf.html) polygon
dataset to an existing distribution object. Presence-absence polygon
data assumes that each area within the polygon can be treated as
'presence' for the species, while each area outside the polygon is where
the species is absent.

## Usage

``` r
add_biodiversity_polpa(
  x,
  polpa,
  name = NULL,
  field_occurrence = "observed",
  formula = NULL,
  family = "binomial",
  link = NULL,
  weight = 1,
  simulate = FALSE,
  simulate_points = 100,
  simulate_bias = NULL,
  simulate_strategy = "random",
  separate_intercept = TRUE,
  docheck = TRUE,
  pseudoabsence_settings = NULL,
  ...
)

# S4 method for class 'BiodiversityDistribution,sf'
add_biodiversity_polpa(
  x,
  polpa,
  name = NULL,
  field_occurrence = "observed",
  formula = NULL,
  family = "binomial",
  link = NULL,
  weight = 1,
  simulate = FALSE,
  simulate_points = 100,
  simulate_bias = NULL,
  simulate_strategy = "random",
  separate_intercept = TRUE,
  docheck = TRUE,
  pseudoabsence_settings = NULL,
  ...
)
```

## Arguments

- x:

  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

- polpa:

  A [`sf::sf`](https://r-spatial.github.io/sf/reference/sf.html) polygon
  object of presence-absence occurrences.

- name:

  The name of the biodiversity dataset used as internal identifier.

- field_occurrence:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) or
  [`character`](https://rdrr.io/r/base/character.html) location of
  biodiversity point records.

- formula:

  A [`character`](https://rdrr.io/r/base/character.html) or
  [`formula`](https://rdrr.io/r/stats/formula.html) object to be passed.
  Default (`NULL`) is to use all covariates .

- family:

  A [`character`](https://rdrr.io/r/base/character.html) stating the
  family to be used (Default: `binomial`).

- link:

  A [`character`](https://rdrr.io/r/base/character.html) to overwrite
  the default link function (Default: `NULL`).

- weight:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) value acting as a
  multiplier with regards to any weights used in the modelling. Larger
  weights indicate higher weighting relative to any other datasets. By
  default set to `1` if only one dataset is added. A
  [`vector`](https://rdrr.io/r/base/vector.html) is also supported but
  must be of the same length as `"polpa"`.

- simulate:

  Simulate poipa points within its boundaries. Results are passed to
  [`add_biodiversity_poipa`](https://iiasa.github.io/ibis.iSDM/reference/add_biodiversity_poipa.md)
  (Default: `FALSE`).

- simulate_points:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) number of points to
  be created by simulation.

- simulate_bias:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  layer describing an eventual preference for simulation (Default:
  `NULL`).

- simulate_strategy:

  A [`character`](https://rdrr.io/r/base/character.html) stating the
  strategy for sampling. Can be set to either `'random'` or `'regular'`,
  the latter requiring a raster supplied in the `'simulate_weights'`
  parameter.

- separate_intercept:

  A [`logical`](https://rdrr.io/r/base/logical.html) value stating
  whether a separate intercept is to be added in shared likelihood
  models for engines
  [engine_inla](https://iiasa.github.io/ibis.iSDM/reference/engine_inla.md),
  [engine_inlabru](https://iiasa.github.io/ibis.iSDM/reference/engine_inlabru.md)
  and
  [engine_stan](https://iiasa.github.io/ibis.iSDM/reference/engine_stan.md).

- docheck:

  [`logical`](https://rdrr.io/r/base/logical.html) on whether additional
  checks should be performed (e.g. intersection tests) (Default:
  `TRUE`).

- pseudoabsence_settings:

  Either `NULL` or a
  [`pseudoabs_settings()`](https://iiasa.github.io/ibis.iSDM/reference/pseudoabs_settings.md)
  created settings object.

- ...:

  Other parameters passed down.

## Value

Adds biodiversity data to
[distribution](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
object.

## Details

The default approach for polygon data is to sample presence-absence
points across the region of the polygons. This function thus acts as a
wrapper to
[`add_biodiversity_poipa()`](https://iiasa.github.io/ibis.iSDM/reference/add_biodiversity_poipa.md)
as presence-absence points are created by the model. Note if the polygon
is used directly in the modelling the link between covariates and
polygonal data is established by regular sampling of points within the
polygon and is thus equivalent to simulating the points directly.

For an integration of range data as predictor or offset, see
[`add_predictor_range()`](https://iiasa.github.io/ibis.iSDM/reference/add_predictor_range.md)
and
[`add_offset_range()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset_range.md)
instead.

## See also

Other add_biodiversity:
[`add_biodiversity_poipa()`](https://iiasa.github.io/ibis.iSDM/reference/add_biodiversity_poipa.md),
[`add_biodiversity_poipo()`](https://iiasa.github.io/ibis.iSDM/reference/add_biodiversity_poipo.md),
[`add_biodiversity_polpo()`](https://iiasa.github.io/ibis.iSDM/reference/add_biodiversity_polpo.md)

## Examples

``` r
if (FALSE) { # \dontrun{
 x <- distribution(background) |>
   add_biodiversity_polpa(protectedArea)
} # }
```
