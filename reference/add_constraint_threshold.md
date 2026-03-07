# Adds a threshold constraint to a scenario object

This option adds a
[`threshold()`](https://iiasa.github.io/ibis.iSDM/reference/threshold.md)
constraint to a scenario projection, thus effectively applying the
threshold as mask to each projection step made during the scenario
projection.

Applying this constraint thus means that the `"suitability"` projection
is clipped to the threshold. This method requires the
[`threshold()`](https://iiasa.github.io/ibis.iSDM/reference/threshold.md)
set for a scenario object.

It could be in theory possible to re calculate the threshold for each
time step based on supplied parameters or even observation records. So
far this option has not been necessary to implement.

## Usage

``` r
add_constraint_threshold(mod, updatevalue = NA, ...)

# S4 method for class 'BiodiversityScenario'
add_constraint_threshold(mod, updatevalue = NA, ...)
```

## Arguments

- mod:

  A
  [`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
  object with specified predictors.

- updatevalue:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) indicating to what
  the masked out values (those outside) the threshold should become
  (Default: `NA`).

- ...:

  passed on parameters. See also the specific methods for adding
  constraints.

## Note

Threshold values are taken from the original fitted model.

## See also

Other constraint:
[`add_constraint()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint.md),
[`add_constraint_MigClim()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_MigClim.md),
[`add_constraint_adaptability()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_adaptability.md),
[`add_constraint_boundary()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_boundary.md),
[`add_constraint_connectivity()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_connectivity.md),
[`add_constraint_dispersal()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_dispersal.md),
[`add_constraint_minsize()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_minsize.md),
[`simulate_population_steps()`](https://iiasa.github.io/ibis.iSDM/reference/simulate_population_steps.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Add scenario constraint
scenario(fit) |> threshold() |>
add_constraint_threshold()
} # }
```
