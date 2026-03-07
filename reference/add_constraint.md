# Add a constraint to an existing `scenario`

This function adds a constraint to a
[`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
object to constrain (future) projections. These constraints can for
instance be constraints on a possible dispersal distance, connectivity
between identified patches or limitations on species adaptability.

**Most constraints require pre-calculated thresholds to be present in
the
[`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
object!**

## Usage

``` r
add_constraint(mod, method, ...)

# S4 method for class 'BiodiversityScenario'
add_constraint(mod, method, ...)
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

- ...:

  passed on parameters. See also the specific methods for adding
  constraints.

## Value

Adds constraints data to a
[`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
object.

## Details

Constraints can be added to scenario objects to increase or decrease the
suitability of a given area for the target feature. This function acts
as a wrapper to add these constraints. Currently supported are the
following options:

**Dispersal**:

- `sdd_fixed` - Applies a fixed uniform dispersal distance per modelling
  timestep.

- `sdd_nexpkernel` - Applies a dispersal distance using a negative
  exponential kernel from its origin.

- `kissmig` - Applies the kissmig stochastic dispersal model. Requires
  `` `kissmig` `` package. Applied at each modelling time step.

- `migclim` - Applies the dispersal algorithm MigClim to the modelled
  objects. Requires `"MigClim"` package.

A comprehensive overview of the benefits of including dispersal
constraints in species distribution models can be found in Bateman et
al. (2013).

**Connectivity**:

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

**Adaptability**:

- `nichelimit` - Specifies a limit on the environmental niche to only
  allow a modest amount of extrapolation beyond the known occurrences.
  This can be particular useful to limit the influence of increasing
  marginal responses and avoid biologically unrealistic projections.

- `fixedlimit` - Sets limit on a given variable based on minimum,
  maximum and preferred range of the variable in question. Common use
  case is for example to constrain future climatically-forced
  projections by the thermal limit of a species.

**Boundary and size**:

- `boundary` - Applies a hard boundary constraint on the projection,
  thus disallowing an expansion of a range outside the provide layer.
  Similar as specifying projection limits (see
  [`distribution`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)),
  but can be used to specifically constrain a projection within a
  certain area (e.g. a species range or an island).

- `minsize` - Allows to specify a certain size that must be satisfied in
  order for a thresholded patch to be occupied. Can be thought of as a
  minimum size requirement. See
  [`add_constraint_minsize()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_minsize.md)
  for the required parameters.

- `threshold` - Applies the set threshold as a constrain directly on the
  suitability projections. Requires a threshold to be set.

## References

- Bateman, B. L., Murphy, H. T., Reside, A. E., Mokany, K., &
  VanDerWal, J. (2013). Appropriateness of full‐, partial‐and
  no‐dispersal scenarios in climate change impact modelling. Diversity
  and Distributions, 19(10), 1224-1234.

- Nobis MP and Normand S (2014) KISSMig - a simple model for R to
  account for limited migration in analyses of species distributions.
  Ecography 37: 1282-1287.

- Mendes, P., Velazco, S. J. E., de Andrade, A. F. A., &
  Júnior, P. D. M. (2020). Dealing with overprediction in species
  distribution models: How adding distance constraints can improve model
  accuracy. Ecological Modelling, 431, 109180.

## See also

Other constraint:
[`add_constraint_MigClim()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_MigClim.md),
[`add_constraint_adaptability()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_adaptability.md),
[`add_constraint_boundary()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_boundary.md),
[`add_constraint_connectivity()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_connectivity.md),
[`add_constraint_dispersal()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_dispersal.md),
[`add_constraint_minsize()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_minsize.md),
[`add_constraint_threshold()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_threshold.md),
[`simulate_population_steps()`](https://iiasa.github.io/ibis.iSDM/reference/simulate_population_steps.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Assumes that a trained 'model' object exists
 mod <- scenario(model) |>
  add_predictors(env = predictors, transform = 'scale', derivates = "none") |>
  add_constraint_dispersal(method = "kissmig", value = 2, pext = 0.1) |>
  project()
} # }
```
