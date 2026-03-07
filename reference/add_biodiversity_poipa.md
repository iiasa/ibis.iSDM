# Add biodiversity point dataset to a distribution object (presence-absence).

This function adds a presence-absence biodiversity dataset to a
distribution object. Opposed to presence-only data, presence-absence
biodiversity records usually originate from structured biodiversity
surveys where the absence of a species in a given region was
specifically assessed.

If it is the analysts choice it is also possible to format presence-only
biodiversity data into a presence-absence form, by adding pseudo-absence
through
[`add_pseudoabsence`](https://iiasa.github.io/ibis.iSDM/reference/add_pseudoabsence.md).
See the help file for more information.

## Usage

``` r
add_biodiversity_poipa(
  x,
  poipa,
  name = NULL,
  field_occurrence = "observed",
  formula = NULL,
  family = "binomial",
  link = NULL,
  weight = 1,
  separate_intercept = TRUE,
  docheck = TRUE,
  ...
)

# S4 method for class 'BiodiversityDistribution,sf'
add_biodiversity_poipa(
  x,
  poipa,
  name = NULL,
  field_occurrence = "observed",
  formula = NULL,
  family = "binomial",
  link = NULL,
  weight = 1,
  separate_intercept = TRUE,
  docheck = TRUE,
  ...
)
```

## Arguments

- x:

  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

- poipa:

  A [`data.frame`](https://rdrr.io/r/base/data.frame.html) or
  [`sf::sf`](https://r-spatial.github.io/sf/reference/sf.html) object of
  presence-absence point occurrences.

- name:

  The name of the biodiversity dataset used as internal identifier.

- field_occurrence:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) or
  [`character`](https://rdrr.io/r/base/character.html) location of
  biodiversity point records indicating presence/absence. By default set
  to `"observed"` and an error will be thrown if a
  [`numeric`](https://rdrr.io/r/base/numeric.html) column with that name
  does not exist.

- formula:

  A [`character`](https://rdrr.io/r/base/character.html) or
  [`formula`](https://rdrr.io/r/stats/formula.html) object to be passed.
  Default (`NULL`) is to use all covariates.

- family:

  A [`character`](https://rdrr.io/r/base/character.html) stating the
  family to be used (Default: `'binomial'`).

- link:

  A [`character`](https://rdrr.io/r/base/character.html) to overwrite
  the default link function (Default: `NULL`).

- weight:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) value acting as a
  multiplier with regards to any weights used in the modelling. Larger
  weights indicate higher weighting relative to any other datasets. By
  default set to `1` if only one dataset is added. A
  [`vector`](https://rdrr.io/r/base/vector.html) is also supported but
  must be of the same length as parameter `"poipa"`.

- separate_intercept:

  A [`logical`](https://rdrr.io/r/base/logical.html) value stating
  whether a separate intercept is to be added in. shared likelihood
  models for engines
  [engine_inla](https://iiasa.github.io/ibis.iSDM/reference/engine_inla.md),
  [engine_inlabru](https://iiasa.github.io/ibis.iSDM/reference/engine_inlabru.md)
  and
  [engine_stan](https://iiasa.github.io/ibis.iSDM/reference/engine_stan.md).

- docheck:

  [`logical`](https://rdrr.io/r/base/logical.html) on whether additional
  checks should be performed (e.g. intersection tests) (Default:
  `TRUE`).

- ...:

  Other parameters passed down.

## Value

Adds biodiversity data to
[distribution](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
object.

## Details

By default, the logit link function is used in a logistic regression
setting unless the specific engine does not support generalised linear
regressions (e.g.
[engine_bart](https://iiasa.github.io/ibis.iSDM/reference/engine_bart.md)).

## References

- Renner, I. W., J. Elith, A. Baddeley, W. Fithian, T. Hastie, S. J.
  Phillips, G. Popovic, and D. I. Warton. 2015. Point process models for
  presence-only analysis. Methods in Ecology and Evolution 6:366–379.

- Guisan A. and Zimmerman N. 2000. Predictive habitat distribution
  models in ecology. Ecol. Model. 135: 147–186.

## See also

Other add_biodiversity:
[`add_biodiversity_poipo()`](https://iiasa.github.io/ibis.iSDM/reference/add_biodiversity_poipo.md),
[`add_biodiversity_polpa()`](https://iiasa.github.io/ibis.iSDM/reference/add_biodiversity_polpa.md),
[`add_biodiversity_polpo()`](https://iiasa.github.io/ibis.iSDM/reference/add_biodiversity_polpo.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Define model
x <- distribution(background) |> add_biodiversity_poipa(virtual_species)
} # }
```
