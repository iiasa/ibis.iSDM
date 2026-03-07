# Adds a boundary constraint to a scenario object

The purpose of boundary constraints is to limit a future projection
within a specified area (such as for example a range or ecoregion). This
can help to limit unreasonable projections into geographic space.

Similar to boundary constraints it is also possible to define a `"zone"`
for the scenario projections, similar as was done for model training.
The difference to a boundary constraint is that the boundary constraint
is applied posthoc as a hard cut on any projection, while the zones
would allow any projection (and other constraints) to be applied within
the zone. **Note: Setting a boundary constraint for future projections
effectively potentially suitable areas!**

## Usage

``` r
add_constraint_boundary(mod, layer, ...)

# S4 method for class 'BiodiversityScenario,sf'
add_constraint_boundary(mod, layer, method = "boundary", ...)

# S4 method for class 'BiodiversityScenario,ANY'
add_constraint_boundary(mod, layer, method = "boundary", ...)
```

## Arguments

- mod:

  A
  [`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
  object with specified predictors.

- layer:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  or [`sf::sf`](https://r-spatial.github.io/sf/reference/sf.html) object
  with the same extent as the model background. Has to be binary and is
  used for a posthoc masking of projected grid cells.

- ...:

  passed on parameters. See also the specific methods for adding
  constraints.

- method:

  A [`character`](https://rdrr.io/r/base/character.html) indicating the
  type of constraints to be added to the scenario. See details for more
  information.

## See also

Other constraint:
[`add_constraint()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint.md),
[`add_constraint_MigClim()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_MigClim.md),
[`add_constraint_adaptability()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_adaptability.md),
[`add_constraint_connectivity()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_connectivity.md),
[`add_constraint_dispersal()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_dispersal.md),
[`add_constraint_minsize()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_minsize.md),
[`add_constraint_threshold()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_threshold.md),
[`simulate_population_steps()`](https://iiasa.github.io/ibis.iSDM/reference/simulate_population_steps.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Add scenario constraint
scenario(fit) |> add_constraint_boundary(range)
} # }
```
