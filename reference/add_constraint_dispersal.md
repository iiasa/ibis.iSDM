# Add dispersal constraint to an existing `scenario`

Add dispersal constraint to an existing `scenario`

## Usage

``` r
add_constraint_dispersal(
  mod,
  method,
  value = NULL,
  unit = "m",
  type = NULL,
  ...
)

# S4 method for class 'BiodiversityScenario'
add_constraint_dispersal(
  mod,
  method,
  value = NULL,
  unit = "m",
  type = NULL,
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

- value:

  For many dispersal `"constrain"` this is set as
  [`numeric`](https://rdrr.io/r/base/numeric.html) value specifying a
  fixed constrain or constant in units `"m"` (Default: `NULL`). For
  kissmig the value needs to give the number of iteration steps (or
  within year migration steps). For adaptability constraints this
  parameter specifies the extent (in units of standard deviation) to
  which extrapolations should be performed.

- unit:

  A [`character`](https://rdrr.io/r/base/character.html) indicating the
  unit of the value parameter. Available are meter (`"m"`) and kilometre
  (`"km"`) (Default: `"m"`).

- type:

  A [`character`](https://rdrr.io/r/base/character.html) indicating the
  type used in the method. See for instance `` `kissmig` ``.

- ...:

  passed on parameters. See also the specific methods for adding
  constraints.

## Details

**Dispersal**: Parameters for `'method'`:

- `sdd_fixed` - Applies a fixed uniform dispersal distance per modelling
  timestep.

- `sdd_nexpkernel` - Applies a dispersal distance using a negative
  exponential kernel from its origin. \#' The negative exponential
  kernel is defined as: \$\$f(x) = \frac{1}{2 \pi a^2}
  e^{-\frac{x}{a}}\$\$ where \\a\\ is the mean dispersal distance (in m)
  divided by 2.

- `kissmig` - Applies the kissmig stochastic dispersal model. Requires
  `` `kissmig` `` package. Applied at each modelling time step.

- `migclim` - Applies the dispersal algorithm MigClim to the modelled
  objects. Requires `"MigClim"` package.

A comprehensive overview of the benefits of including dispersal
constraints in species distribution models can be found in Bateman et
al. (2013).

The following additional parameters can be set:

- `pext`: [`numeric`](https://rdrr.io/r/base/numeric.html) indicator for
  `` `kissmig` `` of the probability a colonized cell becomes
  uncolonised, i.e., the species gets locally extinct (Default: `0.1`).

- `pcor`: [`numeric`](https://rdrr.io/r/base/numeric.html) probability
  that corner cells are considered in the 3x3 neighbourhood (Default:
  `0.2`).

## Note

Unless otherwise stated, the default unit of supplied distance values
(e.g. average dispersal distance) should be in `"m"`.

## References

- Bateman, B. L., Murphy, H. T., Reside, A. E., Mokany, K., &
  VanDerWal, J. (2013). Appropriateness of full‐, partial‐and
  no‐dispersal scenarios in climate change impact modelling. Diversity
  and Distributions, 19(10), 1224-1234.

## See also

Other constraint:
[`add_constraint()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint.md),
[`add_constraint_MigClim()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_MigClim.md),
[`add_constraint_adaptability()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_adaptability.md),
[`add_constraint_boundary()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_boundary.md),
[`add_constraint_connectivity()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_connectivity.md),
[`add_constraint_minsize()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_minsize.md),
[`add_constraint_threshold()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_threshold.md),
[`simulate_population_steps()`](https://iiasa.github.io/ibis.iSDM/reference/simulate_population_steps.md)
