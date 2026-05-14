# Add biodiversity point dataset to a distribution object (presence-only)

This function adds a presence-only biodiversity dataset to a
distribution object.

## Usage

``` r
add_biodiversity_poipo(
  x,
  poipo,
  name = NULL,
  field_occurrence = "observed",
  formula = NULL,
  family = "poisson",
  link = NULL,
  weight = 1,
  separate_intercept = TRUE,
  docheck = TRUE,
  pseudoabsence_settings = NULL,
  ...
)

# S4 method for class 'BiodiversityDistribution,sf'
add_biodiversity_poipo(
  x,
  poipo,
  name = NULL,
  field_occurrence = "observed",
  formula = NULL,
  family = "poisson",
  link = NULL,
  weight = 1,
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

- poipo:

  A [`data.frame`](https://rdrr.io/r/base/data.frame.html) or
  [`sf`](https://r-spatial.github.io/sf/reference/sf.html) object of
  presence-only point occurrences.

- name:

  The name of the biodiversity dataset used as internal identifier.

- field_occurrence:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) or
  [`character`](https://rdrr.io/r/base/character.html) location of
  biodiversity point records.

- formula:

  A [`character`](https://rdrr.io/r/base/character.html) or
  [`formula`](https://rdrr.io/r/stats/formula.html) object to be passed.
  Default is to use all covariates (if specified).

- family:

  A [`character`](https://rdrr.io/r/base/character.html) stating the
  family to be used (Default: `'Poisson'`).

- link:

  A [`character`](https://rdrr.io/r/base/character.html) to overwrite
  the default link function (Default: `NULL`).

- weight:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) value acting as a
  multiplier with regards to any weights used in the modelling. Larger
  weights indicate higher weighting relative to any other datasets. By
  default set to `1` if only one dataset is added. A
  [`vector`](https://rdrr.io/r/base/vector.html) is also supported but
  must be of the same length as `"poipo"`. **Note: Weights are
  reformated to the inverse for models with area offsets (e.g. 5 is
  converted to 1/5).**

- separate_intercept:

  A [`logical`](https://rdrr.io/r/base/logical.html) value stating
  whether a separate intercept is to be added in shared likelihood
  models for engines
  [engine_inla](https://iiasa.github.io/ibis.iSDM/reference/engine_inla.md),
  [engine_inlabru](https://iiasa.github.io/ibis.iSDM/reference/engine_inlabru.md)
  and
  [engine_stan](https://iiasa.github.io/ibis.iSDM/reference/engine_stan.md).
  Otherwise ignored.

- docheck:

  [`logical`](https://rdrr.io/r/base/logical.html) on whether additional
  checks should be performed (e.g. intersection tests) (Default:
  `TRUE`).

- pseudoabsence_settings:

  Either `NULL` or a
  [`pseudoabs_settings()`](https://iiasa.github.io/ibis.iSDM/reference/pseudoabs_settings.md)
  created settings object.

- ...:

  Other parameters passed down to the object. Normally not used unless
  described in details.

## Value

Adds biodiversity data to
[distribution](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
object.

## Details

This function allows to add presence-only biodiversity records to a
[distribution](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
ibis.iSDM Presence-only data are usually modelled through an inferential
model (see Guisan and Zimmerman, 2000) that relate their occurrence in
relation to environmental covariates to a selected sample of
'background' points. The most common approach for estimation and the one
supported by this type of dataset are poisson-process models (PPM) in
which presence-only points are fitted through a down-weighted Poisson
regression. See Renner et al. 2015 for an overview.

## References

- Guisan A. and Zimmerman N. 2000. Predictive habitat distribution
  models in ecology. Ecol. Model. 135: 147–186.

- Renner, I. W., J. Elith, A. Baddeley, W. Fithian, T. Hastie, S. J.
  Phillips, G. Popovic, and D. I. Warton. 2015. Point process models for
  presence-only analysis. Methods in Ecology and Evolution 6:366–379.

## See also

Other add_biodiversity:
[`add_biodiversity_poipa()`](https://iiasa.github.io/ibis.iSDM/reference/add_biodiversity_poipa.md),
[`add_biodiversity_polpa()`](https://iiasa.github.io/ibis.iSDM/reference/add_biodiversity_polpa.md),
[`add_biodiversity_polpo()`](https://iiasa.github.io/ibis.iSDM/reference/add_biodiversity_polpo.md)

## Examples

``` r
# Load background
background <- terra::rast(system.file('extdata/europegrid_50km.tif',
package='ibis.iSDM',mustWork = TRUE))
# Load virtual species
virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg',
package='ibis.iSDM',mustWork = TRUE),'points',quiet = TRUE)
# Define model
x <- distribution(background) |>
add_biodiversity_poipo(virtual_points, field_occurrence = "Observed")
#> [Setup] 2026-05-14 20:28:22.591026 | Creating distribution object...
#> [Setup] 2026-05-14 20:28:22.592192 | Adding poipo dataset...
```
