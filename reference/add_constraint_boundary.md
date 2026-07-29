# Adds a boundary or zone constraint to a scenario object

The purpose of boundary constraints is to limit a future projection
within a specified area (such as for example a range or ecoregion). This
can help to limit unreasonable projections into geographic space.

When `method = "zone"` is used, the constraint is applied per projection
timestep as a mask on the suitability output (before any threshold is
computed), rather than posthoc on the full stacked result. This allows
other per-timestep constraints (e.g. dispersal) to interact with the
zone at each step. Zones can be either static (one layer applied
identically at every timestep) or time-series (one layer per timestep,
matched by nearest date).

**Note: Setting a boundary constraint for future projections effectively
potentially clips suitable areas!**

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

  For `method = "boundary"`: a single-layer
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  or [`sf`](https://r-spatial.github.io/sf/reference/sf.html) object.
  Has to be binary. For `method = "zone"`: must be a
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html).
  A single-layer raster is used as a static zone (same mask at every
  timestep). A multi-layer raster is treated as a time-series zone;
  layers are matched to projection timesteps either by the raster's time
  attribute (`terra::has.time(layer) == TRUE`, nearest date selected)
  or, when no time attribute is present, by layer order
  (`terra::nlyr(layer)` must equal the number of covariate timesteps).

- ...:

  passed on parameters. See also the specific methods for adding
  constraints.

- method:

  A [`character`](https://rdrr.io/r/base/character.html) specifying the
  constraint type. Either `"boundary"` (default, posthoc mask on full
  stacked projection) or `"zone"` (per-timestep mask applied within the
  projection loop).

## Value

A
[`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
object with the boundary or zone constraint added.

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
# Static boundary constraint (posthoc)
scenario(fit) |> add_constraint_boundary(range)

# Static zone constraint (per-timestep)
scenario(fit) |> add_constraint_boundary(range, method = "zone")

# Time-series zone: SpatRaster with terra::time() set
scenario(fit) |> add_constraint_boundary(zone_raster, method = "zone")

} # }
```
