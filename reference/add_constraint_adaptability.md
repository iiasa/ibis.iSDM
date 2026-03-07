# Adds an adaptability constraint to a scenario object

Adaptability constraints assume that suitable habitat for species in
(future) projections might be unsuitable if it is outside the range of
conditions currently observed for the species.

## Usage

``` r
add_constraint_adaptability(
  mod,
  method = "nichelimit",
  names = NULL,
  approach = "thresh",
  value = 1,
  value_min = NULL,
  value_max = NULL,
  increment = 0,
  ...
)

# S4 method for class 'BiodiversityScenario'
add_constraint_adaptability(
  mod,
  method = "nichelimit",
  names = NULL,
  approach = "thresh",
  value = 1,
  value_min = NULL,
  value_max = NULL,
  increment = 0,
  ...
)
```

## Arguments

- mod:

  A
  [`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
  object with specified predictors.

- method:

  A [`character`](https://rdrr.io/r/base/character.html) indicating the
  type of constraints to be added to the scenario. See details for more
  information.

- names:

  A [`character`](https://rdrr.io/r/base/character.html) vector with
  names of the predictors for which an adaptability threshold should be
  set (Default: `NULL` for all).

- approach:

  [`character`](https://rdrr.io/r/base/character.html) on whether
  thresholds or hinges are to be calculated (Default: `'thresh'`). For
  `'fixedlimit'` this controls how strongly the limits are enforced
  (e.g. abrupt or linearly).

- value:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) value in units of
  standard deviation (Default: `1`) or alternatively as prefered value
  for `"fixedlimit"`.

- value_min:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) minimum value used
  for method `"fixedlimit"`.

- value_max:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) maximum value used
  for method `"fixedlimit"`.

- increment:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) constant that is
  added to value at every time step (Default: `0`). Allows incremental
  widening of the niche space, thus opening constraints.

- ...:

  passed on parameters. See also the specific methods for adding
  constraints.

## Details

Currently implemented are the following approaches:

[\*](https://rdrr.io/r/base/Arithmetic.html) `'nichelimit'` = This adds
a simple constrain on the predictor parameter space, which can be
defined through the `"value"` parameter. For example by setting it to
`1` (Default), any projections are constrained to be within the range of
at maximum 1 standard deviation from the range of covariates used for
model training. The parameter `"increment"` furthermore allows to
step-wise increase the value range by a certain amount per time step.

[\*](https://rdrr.io/r/base/Arithmetic.html) `'fixedlimit'` = Here we
can supply a fixed limit for a given variable, provided as a minimum
(`"value_min"`), maximum (`"value_max"`) and preferred (`"value"`) range
in which a biodiversity feature exist. Common applications include for
example known thermal limits with regards to temperature. Internally
this function applies a normalized hinge transform based on the supplied
values.

## See also

Other constraint:
[`add_constraint()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint.md),
[`add_constraint_MigClim()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_MigClim.md),
[`add_constraint_boundary()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_boundary.md),
[`add_constraint_connectivity()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_connectivity.md),
[`add_constraint_dispersal()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_dispersal.md),
[`add_constraint_minsize()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_minsize.md),
[`add_constraint_threshold()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_threshold.md),
[`simulate_population_steps()`](https://iiasa.github.io/ibis.iSDM/reference/simulate_population_steps.md)

## Examples

``` r
if (FALSE) { # \dontrun{
scenario(fit) |>
 add_constraint_adaptability(value = 1)
} # }
```
