# Project a fitted model to a new environment and covariates

Equivalent to
[train](https://iiasa.github.io/ibis.iSDM/reference/train.md), this
function acts as a wrapper to project the model stored in a
[`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
object to newly supplied (future) covariates. Supplied predictors are
usually spatial-temporal predictors which should be prepared via
[`add_predictors()`](https://iiasa.github.io/ibis.iSDM/reference/add_predictors.md)
(e.g. transformations and derivates) in the same way as they have been
during the initial modelling with
[`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md).
Any constraints specified in the scenario object are applied during the
projection.

## Usage

``` r
project.BiodiversityScenario(x, ...)

# S4 method for class 'BiodiversityScenario'
project(
  x,
  date_interpolation = "none",
  stabilize = FALSE,
  stabilize_method = "loess",
  layer = "mean",
  env = NULL,
  verbose = getOption("ibis.setupmessages", default = TRUE),
  ...
)

project.DistributionModel(x, ...)

# S4 method for class 'DistributionModel'
project(
  x,
  env,
  layer = "mean",
  verbose = getOption("ibis.setupmessages", default = TRUE),
  ...
)
```

## Arguments

- x:

  A
  [`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
  object with set predictors. Note that some constrains such as
  `MigClim` can still simulate future change without projections.
  Alternatively a
  [`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
  object can be supplied if other scenario functions are not needed. In
  this case provide an `env` parameter value.

- ...:

  passed on parameters.

- date_interpolation:

  A [`character`](https://rdrr.io/r/base/character.html) on whether
  dates should be interpolated. Options include `"none"` (Default),
  `"annual"`, `"monthly"`, `"daily"`.

- stabilize:

  A [`logical`](https://rdrr.io/r/base/logical.html) value indicating
  whether the suitability projection should be stabilized (Default:
  `FALSE`).

- stabilize_method:

  [`character`](https://rdrr.io/r/base/character.html) stating the
  stabilization method to be applied. Currently supported is
  `` `loess` ``.

- layer:

  A [`character`](https://rdrr.io/r/base/character.html) specifying the
  layer to be projected (Default: `"mean"`).

- env:

  An optional
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  or [`data.frame`](https://rdrr.io/r/base/data.frame.html) object for
  prediction. Ignored unless a
  [`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
  object is supplied (Default: `NULL`).

- verbose:

  Setting this [`logical`](https://rdrr.io/r/base/logical.html) value to
  `TRUE` prints out further information during the model fitting
  (Default: `FALSE`).

## Value

Saves [`stars`](https://rdrr.io/r/graphics/stars.html) objects of the
obtained predictions in mod.

## Details

In the background the function `x$project()` for the respective model
object is called, where `x` is fitted model object. For specifics on the
constraints, see the relevant `constrain` functions, respectively:

- [`add_constraint()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint.md)
  for generic wrapper to add any of the available constrains.

- [`add_constraint_dispersal()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_dispersal.md)
  for specifying dispersal constraint on the temporal projections at
  each step.

- [`add_constraint_MigClim()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_MigClim.md)
  Using the MigClim R-package to simulate dispersal in projections.

- [`add_constraint_connectivity()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_connectivity.md)
  Apply a connectivity constraint at the projection, for instance by
  adding a barrier that prevents migration.

- [`add_constraint_minsize()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_minsize.md)
  Adds a constraint on the minimum area a given thresholded patch should
  have, assuming that smaller areas are in fact not suitable.

- [`add_constraint_adaptability()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_adaptability.md)
  Apply an adaptability constraint to the projection, for instance
  constraining the speed a species is able to adapt to new conditions.

- [`add_constraint_boundary()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_boundary.md)
  To artificially limit the distribution change. Similar as specifying
  projection limits, but can be used to specifically constrain a
  projection within a certain area (e.g. a species range or an island).

Many constrains also requires thresholds to be calculated. Adding
[`threshold()`](https://iiasa.github.io/ibis.iSDM/reference/threshold.md)
to a
[`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
object enables the computation of thresholds at every step based on the
threshold used for the main model (threshold values are taken from
there).

It is also possible to make a complementary simulation with the `steps`
package, which can be provided via
[`simulate_population_steps()`](https://iiasa.github.io/ibis.iSDM/reference/simulate_population_steps.md)
to the
[`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
object. Similar as with thresholds, estimates values will then be added
to the outputs.

Finally this function also allows temporal stabilization across
prediction steps via enabling the parameter `stabilize` and checking the
`stabilize_method` argument. Stabilization can for instance be helpful
in situations where environmental variables are quite dynamic, but
changes in projected suitability are not expected to abruptly increase
or decrease. It is thus a way to smoothen out outliers from the
projection. Options are so far for instance `'loess'` which fits a
[`loess()`](https://rdrr.io/r/stats/loess.html) model per pixel and time
step. This is conducted at the very end of the processing steps and any
thresholds will be recalculated afterwards.

## See also

[`scenario()`](https://iiasa.github.io/ibis.iSDM/reference/scenario.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit a model
fit <- distribution(background) |>
        add_biodiversity_poipa(surveydata) |>
        add_predictors(env = predictors) |>
        engine_breg() |>
        train()

# Fit a scenario
sc <- scenario(fit) |>
        add_predictors(env = future_predictors) |>
        project()
} # }
```
