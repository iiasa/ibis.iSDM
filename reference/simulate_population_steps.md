# Simulate population dynamics following the steps approach

This function adds a flag to a
[`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
object to indicate that species abundances are to be simulated based on
the expected habitat suitability, as well as demography,
density-dependence and dispersal information. The simulation is done
using the *steps* package (Visintin et al. 2020) and conducted after a
habitat suitability projection has been created. *steps* is a spatially
explicit population models coded mostly in R.

For a detailed description of *steps* parameters, please see the
respective reference and help files. Default assumptions underlying this
wrapper are presented in the details

## Usage

``` r
simulate_population_steps(
  mod,
  vital_rates,
  replicates = 1,
  carrying_capacity = NULL,
  initial = NULL,
  dispersal = NULL,
  density_dependence = NULL,
  include_suitability = TRUE
)

# S4 method for class 'BiodiversityScenario,matrix'
simulate_population_steps(
  mod,
  vital_rates,
  replicates = 1,
  carrying_capacity = NULL,
  initial = NULL,
  dispersal = NULL,
  density_dependence = NULL,
  include_suitability = TRUE
)
```

## Arguments

- mod:

  A
  [`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
  object with specified predictors.

- vital_rates:

  A symmetrical demographic matrix. Should have column and row names
  equivalent to the vital stages that are to be estimated.

- replicates:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) vector of the
  number of replicates (Default: `1`).

- carrying_capacity:

  Either
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  or a [`numeric`](https://rdrr.io/r/base/numeric.html) estimate of the
  maximum carrying capacity, e.g. how many adult individual are likely
  to occur per grid cell. If set to
  [`numeric`](https://rdrr.io/r/base/numeric.html), then carrying
  capacity is estimated up to a maximum set (*Note: a more clever way
  would be to use a species-area relationship for scaling. This is not
  yet implemented*).

- initial:

  A
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  giving the initial population size. If not provided, then initial
  populations are guessed (see details) from the projected suitability
  rasters (Default: `NULL`).

- dispersal:

  A dispersal object defined by the `steps` package (Default: `NULL`).

- density_dependence:

  Specification of density dependence defined by the `steps` package
  (Default: `NULL`).

- include_suitability:

  A [`logical`](https://rdrr.io/r/base/logical.html) flag on whether the
  projected suitability estimates should be used (Default: `TRUE`) or
  only the initial conditions set to the first time step.

## Value

Adds flag to a
[`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
object to indicate that further simulations are added during projection.

## Details

In order for this function to work the *steps* package has to be
installed separately. Instructions to do so can be found on
[github](https://github.com/steps-dev/steps).

If initial population lifestages are not provided, then they are
estimated assuming a linear scaling with suitability, a `50:50` split
between sexes and a `1:3` ratio of adults to juveniles. The provision of
different parameters is highly encouraged!

## Note

The *steps* package has multiple options for simulating species
population and not all possible options are represented in this wrapper.

Furthermore, the package still makes use of the `raster` package for
much of its internal data processing. Since *ibis.iSDM* switched to
[terra](https://rspatial.github.io/terra/reference/terra-package.html) a
while ago, there can be efficiency problems as layers need to be
translated between packages.

## References

- Visintin, C., Briscoe, N. J., Woolley, S. N., Lentini, P. E., Tingley,
  R., Wintle, B. A., & Golding, N. (2020). steps: Software for spatially
  and temporally explicit population simulations. Methods in Ecology and
  Evolution, 11(4), 596-603. https://doi.org/10.1111/2041-210X.13354

## See also

Other constraint:
[`add_constraint()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint.md),
[`add_constraint_MigClim()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_MigClim.md),
[`add_constraint_adaptability()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_adaptability.md),
[`add_constraint_boundary()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_boundary.md),
[`add_constraint_connectivity()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_connectivity.md),
[`add_constraint_dispersal()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_dispersal.md),
[`add_constraint_minsize()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_minsize.md),
[`add_constraint_threshold()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_threshold.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Define vital rates
vt <- matrix(c(0.0,0.5,0.75,
               0.5,0.2,0.0,
               0.0,0.5,0.9),
               nrow = 3, ncol = 3, byrow = TRUE)
colnames(vt) <- rownames(vt) <- c('juvenile','subadult','adult')

# Assumes that a trained 'model' object exists
 mod <- scenario(model) |>
  add_predictors(env = predictors, transform = 'scale',
                 derivates = "none") |>
  # Use Vital rates here, but note the other parameters!
  simulate_population_steps(vital_rates = vt) |>
  project()
} # }
```
