# Add constraints to the modelled distribution projection using the MigClim approach

This function adds a constraint as defined by the MigClim approach
(Engler et al. 2013) to a
[`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
object to constrain future projections. For a detailed description of
MigClim, please see the respective reference and the UserGuide. **The
default parameters chosen here are suggestions.**

## Usage

``` r
add_constraint_MigClim(
  mod,
  rcThresholdMode = "continuous",
  dispSteps = 1,
  dispKernel = c(1, 0.4, 0.16, 0.06, 0.03),
  barrierType = "strong",
  lddFreq = 0,
  lddRange = c(1000, 10000),
  iniMatAge = 1,
  propaguleProdProb = c(0.2, 0.6, 0.8, 0.95),
  replicateNb = 10,
  dtmp = terra::terraOptions(print = F)$tempdir
)

# S4 method for class 'BiodiversityScenario'
add_constraint_MigClim(
  mod,
  rcThresholdMode = "continuous",
  dispSteps = 1,
  dispKernel = c(1, 0.4, 0.16, 0.06, 0.03),
  barrierType = "strong",
  lddFreq = 0,
  lddRange = c(1000, 10000),
  iniMatAge = 1,
  propaguleProdProb = c(0.2, 0.6, 0.8, 0.95),
  replicateNb = 10,
  dtmp = terra::terraOptions(print = F)$tempdir
)
```

## Arguments

- mod:

  A
  [`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
  object with specified predictors.

- rcThresholdMode:

  A [`character`](https://rdrr.io/r/base/character.html) of either
  **binary** or **continuous** value (Default: **continuous**).

- dispSteps:

  [`numeric`](https://rdrr.io/r/base/numeric.html) parameters on the
  number of dispersal steps. Dispersal steps are executed for each
  timestep (prediction layer). and ideally should be aligned with the
  number of steps for projection. Minimum is `1` (Default) and maximum
  is `99`.

- dispKernel:

  A [`vector`](https://rdrr.io/r/base/vector.html) with the number of
  the dispersal Kernel to be applied. Can be set either to a uniform
  numeric [vector](https://rdrr.io/r/base/vector.html), e.g.
  `c(1,1,1,1)` or to a proportional decline `(1,0.4,0.16,0.06,0.03)`
  (Default). **Depending on the resolution of the raster, this parameter
  needs to be adapted**

- barrierType:

  A [character](https://rdrr.io/r/base/character.html) indicating
  whether any set barrier should be set as `'strong'` or `'weak'`
  barriers. Strong barriers prevent any dispersal across the barrier and
  weak barriers only do so if the whole `"dispKernel"` length is covered
  by the barrier (Default: `'strong'`).

- lddFreq:

  [`numeric`](https://rdrr.io/r/base/numeric.html) parameter indicating
  the frequency of long-distance dispersal (LDD) events. Default is `0`,
  so no long-distance dispersal.

- lddRange:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) value highlighting
  the minimum and maximum distance of LDD events. **Note: The units for
  those distance are in cells, thus the projection units in the
  raster.**

- iniMatAge:

  Initial maturity age. Used together with `propaguleProd` as a proxy of
  population growth. Must be set to the cell age in time units which are
  dispersal steps (Default: `1`).

- propaguleProdProb:

  Probability of a source cell to produce propagules as a function of
  time since colonization. Set as probability vector that defines the
  probability of a cell producing propagules.

- replicateNb:

  Number of replicates to be used for the analysis (Default: `10`).

- dtmp:

  A [`character`](https://rdrr.io/r/base/character.html) to a folder
  where temporary files are to be created.

## Value

Adds a MigClim constraint to a
[`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
object.

## Details

The barrier parameter is defined through `"add_barrier"`.

## References

- Engler R., Hordijk W. and Guisan A. The MIGCLIM R package – seamless
  integration of dispersal constraints into projections of species
  distribution models. Ecography,

- Robin Engler, Wim Hordijk and Loic Pellissier (2013). MigClim:
  Implementing dispersal into species distribution models. R package
  version 1.6.

## See also

Other constraint:
[`add_constraint()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint.md),
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
  add_predictors(env = predictors, transform = 'scale',
                 derivates = "none") |>
  add_constraint_MigClim() |>
  project()
} # }
```
