# Adds a size constraint on a scenario

This function applies a minimum size constraint on a
[`scenario()`](https://iiasa.github.io/ibis.iSDM/reference/scenario.md)
created object. The rationale here is that for a given species isolated
habitat patches smaller than a given size might not be viable /
unrealistic for a species to establish a (long-term) presence.

The idea thus is to apply a constraint in that only patches bigger than
a certain size are retained between timesteps. It has thus the potential
to reduce subsequent colonizations of neighbouring patches.

## Usage

``` r
add_constraint_minsize(
  mod,
  value,
  unit = "km2",
  establishment_step = FALSE,
  ...
)

# S4 method for class 'BiodiversityScenario,numeric'
add_constraint_minsize(
  mod,
  value,
  unit = "km2",
  establishment_step = FALSE,
  ...
)
```

## Arguments

- mod:

  A
  [`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
  object with specified predictors.

- value:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) value describing
  the minimum amount of area of a given patch

- unit:

  A [`character`](https://rdrr.io/r/base/character.html) of the unit of
  area. Options available are `km2` (Default), `ha` and `pixel`.

- establishment_step:

  A [`logical`](https://rdrr.io/r/base/logical.html) flag indicating
  whether a given patch is only to be removed if wasn't small in a
  previous time step (not yet implemented!)

- ...:

  passed on parameters. See also the specific methods for adding
  constraints.

## Value

A
[`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
object with the minimum-size constraint added.

## Details

Area values in a specific unit need to be supplied.

## Note

*This function requires that a scenario has a set
[`threshold()`](https://iiasa.github.io/ibis.iSDM/reference/threshold.md)!*

## See also

Other constraint:
[`add_constraint()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint.md),
[`add_constraint_MigClim()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_MigClim.md),
[`add_constraint_adaptability()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_adaptability.md),
[`add_constraint_boundary()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_boundary.md),
[`add_constraint_connectivity()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_connectivity.md),
[`add_constraint_dispersal()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_dispersal.md),
[`add_constraint_threshold()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_threshold.md),
[`simulate_population_steps()`](https://iiasa.github.io/ibis.iSDM/reference/simulate_population_steps.md)

## Examples

``` r
if (FALSE) { # \dontrun{
scenario(fit) |>
 add_predictors(future_covariates) |>
 threshold() |>
 add_constraint_minsize(value = 1000, unit = "km2") |>
 project()
} # }
```
