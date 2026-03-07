# Adds a connectivity constraint to a scenario object.

Adds a connectivity constraint to a scenario object.

## Usage

``` r
add_constraint_connectivity(mod, method, value = NULL, resistance = NULL, ...)

# S4 method for class 'BiodiversityScenario'
add_constraint_connectivity(mod, method, value = NULL, resistance = NULL, ...)
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

- value:

  For many dispersal `"constrain"` this is set as
  [`numeric`](https://rdrr.io/r/base/numeric.html) value specifying a
  fixed constrain or constant in units `"m"` (Default: `NULL`). For
  kissmig the value needs to give the number of iteration steps (or
  within year migration steps). For adaptability constraints this
  parameter specifies the extent (in units of standard deviation) to
  which extrapolations should be performed.

- resistance:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object describing a resistance surface or barrier for use in
  connectivity constrains (Default: `NULL`).

- ...:

  passed on parameters. See also the specific methods for adding
  constraints.

## Details

- `hardbarrier` - Defines a hard barrier to any dispersal events. By
  definition this sets all values larger than `0` in the barrier layer
  to `0` in the projection. Barrier has to be provided through the
  `"resistance"` parameter.

- `resistance` - Allows the provision of a static or dynamic layer that
  is multiplied with the projection at each time step. Can for example
  be used to reduce the suitability of any given area (using pressures
  not included in the model). The respective layer(s) have to be
  provided through the `"resistance"` parameter. Provided layers are
  incorporated as `abs(resistance - 1)` and multiplied with the
  prediction.

## See also

Other constraint:
[`add_constraint()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint.md),
[`add_constraint_MigClim()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_MigClim.md),
[`add_constraint_adaptability()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_adaptability.md),
[`add_constraint_boundary()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_boundary.md),
[`add_constraint_dispersal()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_dispersal.md),
[`add_constraint_minsize()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_minsize.md),
[`add_constraint_threshold()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_threshold.md),
[`simulate_population_steps()`](https://iiasa.github.io/ibis.iSDM/reference/simulate_population_steps.md)
