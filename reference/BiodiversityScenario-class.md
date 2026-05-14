# Class for a biodiversity scenario from a trained model

Base [`R6`](https://r6.r-lib.org/reference/R6Class.html) class for any
biodiversity scenario objects. Serves as container that supplies data
and functions to other
[`R6`](https://r6.r-lib.org/reference/R6Class.html) classes and
functions.

## Note

This sets the threshold method internally to `'fixed'`.

The latent factor is usually obtained from the fitted model object,
unless re-specified and added here to the list.

This requires the `"gganimate"` package.

This requires a set
[`threshold()`](https://iiasa.github.io/ibis.iSDM/reference/threshold.md)
to the scenario object.

This requires set
[`threshold`](https://iiasa.github.io/ibis.iSDM/reference/threshold.md)
prior to projection.

## See also

[`threshold()`](https://iiasa.github.io/ibis.iSDM/reference/threshold.md)

[`add_latent_spatial()`](https://iiasa.github.io/ibis.iSDM/reference/add_latent_spatial.md)

## Public fields

- `modelobject`:

  A name of the model for projection.

- `modelid`:

  An id of the model used for projection.

- `limits`:

  A [`sf`](https://r-spatial.github.io/sf/reference/sf.html) object used
  to constraint the prediction.

- `predictors`:

  A predictor object for projection.

- `constraints`:

  Any constraints set for projection.

- `latentfactors`:

  A [`list`](https://rdrr.io/r/base/list.html) on whether latentfactors
  are used.

- `scenarios`:

  The resulting [`stars`](https://rdrr.io/r/graphics/stars.html)
  objects.

- `log`:

  A [`Log`](https://iiasa.github.io/ibis.iSDM/reference/Log-class.md)
  object for capturing messages.

## Methods

### Public methods

- [`BiodiversityScenario$new()`](#method-BiodiversityScenario-initialize)

- [`BiodiversityScenario$print()`](#method-BiodiversityScenario-print)

- [`BiodiversityScenario$verify()`](#method-BiodiversityScenario-verify)

- [`BiodiversityScenario$show()`](#method-BiodiversityScenario-show)

- [`BiodiversityScenario$get_projection()`](#method-BiodiversityScenario-get_projection)

- [`BiodiversityScenario$get_resolution()`](#method-BiodiversityScenario-get_resolution)

- [`BiodiversityScenario$get_model()`](#method-BiodiversityScenario-get_model)

- [`BiodiversityScenario$get_limits()`](#method-BiodiversityScenario-get_limits)

- [`BiodiversityScenario$rm_limits()`](#method-BiodiversityScenario-rm_limits)

- [`BiodiversityScenario$get_predictor_names()`](#method-BiodiversityScenario-get_predictor_names)

- [`BiodiversityScenario$get_timeperiod()`](#method-BiodiversityScenario-get_timeperiod)

- [`BiodiversityScenario$get_constraints()`](#method-BiodiversityScenario-get_constraints)

- [`BiodiversityScenario$rm_constraints()`](#method-BiodiversityScenario-rm_constraints)

- [`BiodiversityScenario$get_threshold()`](#method-BiodiversityScenario-get_threshold)

- [`BiodiversityScenario$get_thresholdvalue()`](#method-BiodiversityScenario-get_thresholdvalue)

- [`BiodiversityScenario$apply_threshold()`](#method-BiodiversityScenario-apply_threshold)

- [`BiodiversityScenario$set_predictors()`](#method-BiodiversityScenario-set_predictors)

- [`BiodiversityScenario$set_constraints()`](#method-BiodiversityScenario-set_constraints)

- [`BiodiversityScenario$get_log()`](#method-BiodiversityScenario-get_log)

- [`BiodiversityScenario$set_log()`](#method-BiodiversityScenario-set_log)

- [`BiodiversityScenario$get_simulation()`](#method-BiodiversityScenario-get_simulation)

- [`BiodiversityScenario$set_simulation()`](#method-BiodiversityScenario-set_simulation)

- [`BiodiversityScenario$get_predictors()`](#method-BiodiversityScenario-get_predictors)

- [`BiodiversityScenario$rm_predictors()`](#method-BiodiversityScenario-rm_predictors)

- [`BiodiversityScenario$get_data()`](#method-BiodiversityScenario-get_data)

- [`BiodiversityScenario$rm_data()`](#method-BiodiversityScenario-rm_data)

- [`BiodiversityScenario$set_data()`](#method-BiodiversityScenario-set_data)

- [`BiodiversityScenario$set_latent()`](#method-BiodiversityScenario-set_latent)

- [`BiodiversityScenario$get_latent()`](#method-BiodiversityScenario-get_latent)

- [`BiodiversityScenario$rm_latent()`](#method-BiodiversityScenario-rm_latent)

- [`BiodiversityScenario$plot()`](#method-BiodiversityScenario-plot)

- [`BiodiversityScenario$plot_threshold()`](#method-BiodiversityScenario-plot_threshold)

- [`BiodiversityScenario$plot_migclim()`](#method-BiodiversityScenario-plot_migclim)

- [`BiodiversityScenario$plot_animation()`](#method-BiodiversityScenario-plot_animation)

- [`BiodiversityScenario$plot_relative_change()`](#method-BiodiversityScenario-plot_relative_change)

- [`BiodiversityScenario$summary()`](#method-BiodiversityScenario-summary)

- [`BiodiversityScenario$summary_beforeafter()`](#method-BiodiversityScenario-summary_beforeafter)

- [`BiodiversityScenario$plot_scenarios_slope()`](#method-BiodiversityScenario-plot_scenarios_slope)

- [`BiodiversityScenario$calc_scenarios_slope()`](#method-BiodiversityScenario-calc_scenarios_slope)

- [`BiodiversityScenario$mask()`](#method-BiodiversityScenario-mask)

- [`BiodiversityScenario$get_centroid()`](#method-BiodiversityScenario-get_centroid)

- [`BiodiversityScenario$save()`](#method-BiodiversityScenario-save)

- [`BiodiversityScenario$clone()`](#method-BiodiversityScenario-clone)

------------------------------------------------------------------------

### `BiodiversityScenario$new()`

Initializes the object and creates an empty list

#### Usage

    BiodiversityScenario$new()

#### Returns

NULL

------------------------------------------------------------------------

### `BiodiversityScenario$print()`

Print the names and properties of all scenarios.

#### Usage

    BiodiversityScenario$print()

#### Returns

A message on screen

------------------------------------------------------------------------

### `BiodiversityScenario$verify()`

Verify that set Model exist and check self-validity

#### Usage

    BiodiversityScenario$verify()

#### Returns

Invisible

------------------------------------------------------------------------

### `BiodiversityScenario$show()`

Show the name of the Model

#### Usage

    BiodiversityScenario$show()

#### Returns

Model objectname

------------------------------------------------------------------------

### `BiodiversityScenario$get_projection()`

Get projection of the projection.

#### Usage

    BiodiversityScenario$get_projection()

#### Returns

A [`sf`](https://r-spatial.github.io/sf/reference/sf.html) object with
the geographic projection

------------------------------------------------------------------------

### `BiodiversityScenario$get_resolution()`

Get resultion of the projection.

#### Usage

    BiodiversityScenario$get_resolution()

#### Returns

A [`numeric`](https://rdrr.io/r/base/numeric.html) indication of the
resolution.

------------------------------------------------------------------------

### `BiodiversityScenario$get_model()`

Get the actual model used for projection

#### Usage

    BiodiversityScenario$get_model(copy = FALSE)

#### Arguments

- `copy`:

  A [`logical`](https://rdrr.io/r/base/logical.html) flag on whether a
  deep copy should be created.

#### Returns

A DistributionModel object.

------------------------------------------------------------------------

### `BiodiversityScenario$get_limits()`

Get provided projection limits if set.

#### Usage

    BiodiversityScenario$get_limits()

#### Returns

A [`sf`](https://r-spatial.github.io/sf/reference/sf.html) object or
NULL.

------------------------------------------------------------------------

### `BiodiversityScenario$rm_limits()`

Remove current limits.

#### Usage

    BiodiversityScenario$rm_limits()

#### Returns

Invisible

------------------------------------------------------------------------

### `BiodiversityScenario$get_predictor_names()`

Get names of predictors for scenario object.

#### Usage

    BiodiversityScenario$get_predictor_names()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) vector with the
names.

------------------------------------------------------------------------

### `BiodiversityScenario$get_timeperiod()`

Get time period of projection.

#### Usage

    BiodiversityScenario$get_timeperiod(what = "range")

#### Arguments

- `what`:

  [`character`](https://rdrr.io/r/base/character.html) on whether full
  time period or just the range is to be returned.

#### Returns

A time period from start to end.

------------------------------------------------------------------------

### `BiodiversityScenario$get_constraints()`

Get constrains for model

#### Usage

    BiodiversityScenario$get_constraints()

#### Returns

A [`list`](https://rdrr.io/r/base/list.html) with the constraints within
the scenario.

------------------------------------------------------------------------

### `BiodiversityScenario$rm_constraints()`

Remove contraints from model

#### Usage

    BiodiversityScenario$rm_constraints()

#### Returns

Invisible

------------------------------------------------------------------------

### `BiodiversityScenario$get_threshold()`

Get thresholds if specified.

#### Usage

    BiodiversityScenario$get_threshold()

#### Returns

A [`list`](https://rdrr.io/r/base/list.html) with method and value for
the threshold.

------------------------------------------------------------------------

### `BiodiversityScenario$get_thresholdvalue()`

Duplicate function for internal consistency to return threshold

#### Usage

    BiodiversityScenario$get_thresholdvalue()

#### Returns

A [`list`](https://rdrr.io/r/base/list.html) with method and value for
the threshold.

------------------------------------------------------------------------

### `BiodiversityScenario$apply_threshold()`

Apply a new threshold to the projection.

#### Usage

    BiodiversityScenario$apply_threshold(tr = new_waiver())

#### Arguments

- `tr`:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) value with the new
  threshold.

#### Returns

This object.

------------------------------------------------------------------------

### `BiodiversityScenario$set_predictors()`

Set new predictors to this object.

#### Usage

    BiodiversityScenario$set_predictors(x)

#### Arguments

- `x`:

  [`PredictorDataset`](https://iiasa.github.io/ibis.iSDM/reference/PredictorDataset-class.md)
  object to be supplied.

#### Returns

This object.

------------------------------------------------------------------------

### `BiodiversityScenario$set_constraints()`

Set new constrains

#### Usage

    BiodiversityScenario$set_constraints(x)

#### Arguments

- `x`:

  A [`list`](https://rdrr.io/r/base/list.html) object with constraint
  settings.

#### Returns

This object.

------------------------------------------------------------------------

### `BiodiversityScenario$get_log()`

Returns the output filename of the current log object if set.

#### Usage

    BiodiversityScenario$get_log()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) where the output
is returned.

------------------------------------------------------------------------

### `BiodiversityScenario$set_log()`

Set a new log object

#### Usage

    BiodiversityScenario$set_log(x)

#### Arguments

- `x`:

  A [`Log`](https://iiasa.github.io/ibis.iSDM/reference/Log-class.md)
  object.

#### Returns

This object

------------------------------------------------------------------------

### `BiodiversityScenario$get_simulation()`

Get simulation options and parameters if gound

#### Usage

    BiodiversityScenario$get_simulation()

#### Returns

A [`list`](https://rdrr.io/r/base/list.html) with the parameters.

------------------------------------------------------------------------

### `BiodiversityScenario$set_simulation()`

Set simulation objects.

#### Usage

    BiodiversityScenario$set_simulation(x)

#### Arguments

- `x`:

  new simulation entries and options as
  [`list`](https://rdrr.io/r/base/list.html) to be set.

#### Returns

This object.

------------------------------------------------------------------------

### `BiodiversityScenario$get_predictors()`

Get Predictors from the object.

#### Usage

    BiodiversityScenario$get_predictors()

#### Returns

A predictor dataset.

------------------------------------------------------------------------

### `BiodiversityScenario$rm_predictors()`

Remove predictors from the object.

#### Usage

    BiodiversityScenario$rm_predictors(names)

#### Arguments

- `names`:

  A [`character`](https://rdrr.io/r/base/character.html) vector with
  names

#### Returns

This object.

------------------------------------------------------------------------

### `BiodiversityScenario$get_data()`

Get scenario predictions or any other data

#### Usage

    BiodiversityScenario$get_data(what = "scenarios")

#### Arguments

- `what`:

  A [`character`](https://rdrr.io/r/base/character.html) vector with
  names of what

#### Returns

This object.

------------------------------------------------------------------------

### `BiodiversityScenario$rm_data()`

Remove scenario predictions

#### Usage

    BiodiversityScenario$rm_data()

#### Returns

Invisible

------------------------------------------------------------------------

### `BiodiversityScenario$set_data()`

Set new data in object.

#### Usage

    BiodiversityScenario$set_data(x)

#### Arguments

- `x`:

  A new data object measuing scenarios.

#### Returns

This object.

------------------------------------------------------------------------

### `BiodiversityScenario$set_latent()`

Adding latent factors to the object.

#### Usage

    BiodiversityScenario$set_latent(latent)

#### Arguments

- `latent`:

  A [`list`](https://rdrr.io/r/base/list.html) containing the data
  object.

#### Returns

This object.

------------------------------------------------------------------------

### `BiodiversityScenario$get_latent()`

Get latent factors if found in object.

#### Usage

    BiodiversityScenario$get_latent()

#### Returns

A [`list`](https://rdrr.io/r/base/list.html) with the latent settings

------------------------------------------------------------------------

### `BiodiversityScenario$rm_latent()`

Remove latent factors if found in object.

#### Usage

    BiodiversityScenario$rm_latent()

#### Returns

This object.

------------------------------------------------------------------------

### `BiodiversityScenario$plot()`

Plot the predictions made here.

#### Usage

    BiodiversityScenario$plot(what = "suitability", which = NULL, ...)

#### Arguments

- `what`:

  A [`character`](https://rdrr.io/r/base/character.html) describing the
  layers to be plotted.

- `which`:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) subset to any
  specific time steps.

- `...`:

  Any other parameters passed on.

#### Returns

A graphical representation

------------------------------------------------------------------------

### `BiodiversityScenario$plot_threshold()`

Convenience function to plot thresholds if set

#### Usage

    BiodiversityScenario$plot_threshold(which = NULL)

#### Arguments

- `which`:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) subset to any
  specific time steps.

#### Returns

A graphical representation

------------------------------------------------------------------------

### `BiodiversityScenario$plot_migclim()`

Plot Migclim results if existing.

#### Usage

    BiodiversityScenario$plot_migclim()

#### Returns

A graphical representation

------------------------------------------------------------------------

### `BiodiversityScenario$plot_animation()`

Plot animation of scenarios if possible

#### Usage

    BiodiversityScenario$plot_animation(what = "suitability", fname = NULL)

#### Arguments

- `what`:

  A [`character`](https://rdrr.io/r/base/character.html) describing the
  layers to be plotted.

- `fname`:

  An optional filename to write the result.

#### Returns

A graphical representation

------------------------------------------------------------------------

### `BiodiversityScenario$plot_relative_change()`

Plot relative change between baseline and projected thresholds

#### Usage

    BiodiversityScenario$plot_relative_change(
      position = NULL,
      variable = "mean",
      plot = TRUE
    )

#### Arguments

- `position`:

  Which layer to be plotted

- `variable`:

  A [`character`](https://rdrr.io/r/base/character.html) of the variable
  to be plotted

- `plot`:

  [`logical`](https://rdrr.io/r/base/logical.html) flag on whether to
  plot the results or return the object.

#### Returns

A graphical representation or
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html).

------------------------------------------------------------------------

### `BiodiversityScenario$summary()`

Summarize the change in layers between timesteps

#### Usage

    BiodiversityScenario$summary(
      layer = "threshold",
      plot = FALSE,
      relative = FALSE
    )

#### Arguments

- `layer`:

  A [`character`](https://rdrr.io/r/base/character.html) of the variable
  to be plotted

- `plot`:

  [`logical`](https://rdrr.io/r/base/logical.html) flag on whether to
  plot the results or return the coefficients.

- `relative`:

  [`logical`](https://rdrr.io/r/base/logical.html) on coefficients to be
  converted to relative change.

#### Returns

Summarized coefficients as
[`data.frame`](https://rdrr.io/r/base/data.frame.html)

------------------------------------------------------------------------

### `BiodiversityScenario$summary_beforeafter()`

Summarize before-after change of first and last layer.

#### Usage

    BiodiversityScenario$summary_beforeafter()

#### Returns

Summarized coefficients as
[`data.frame`](https://rdrr.io/r/base/data.frame.html)

------------------------------------------------------------------------

### `BiodiversityScenario$plot_scenarios_slope()`

Calculate slopes across the projection

#### Usage

    BiodiversityScenario$plot_scenarios_slope(
      what = "suitability",
      oftype = "stars"
    )

#### Arguments

- `what`:

  A [`character`](https://rdrr.io/r/base/character.html) with layer to
  be plotted (default: `"suitability"`).

- `oftype`:

  [`character`](https://rdrr.io/r/base/character.html) of the output
  type.

#### Returns

A plot of the scenario slopes

------------------------------------------------------------------------

### `BiodiversityScenario$calc_scenarios_slope()`

Calculate slopes across the projection

#### Usage

    BiodiversityScenario$calc_scenarios_slope(
      what = "suitability",
      plot = TRUE,
      oftype = "stars"
    )

#### Arguments

- `what`:

  A [`character`](https://rdrr.io/r/base/character.html) with layer to
  be plotted (default: `"suitability"`).

- `plot`:

  [`logical`](https://rdrr.io/r/base/logical.html) flag on whether to
  plot the results or return the coefficients.

- `oftype`:

  [`character`](https://rdrr.io/r/base/character.html) of the output
  type.

#### Returns

A
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
layer or [`stars`](https://rdrr.io/r/graphics/stars.html) object.

------------------------------------------------------------------------

### `BiodiversityScenario$mask()`

Convenience function to mask all input projections.

#### Usage

    BiodiversityScenario$mask(mask, inverse = FALSE, ...)

#### Arguments

- `mask`:

  A `SpatRaster` or `sf` object.

- `inverse`:

  A `logical` flag if the inverse should be masked instead.

- `...`:

  Any other parameters passed on.

#### Returns

Invisible

------------------------------------------------------------------------

### `BiodiversityScenario$get_centroid()`

Get centroids of projection layers

#### Usage

    BiodiversityScenario$get_centroid(patch = FALSE)

#### Arguments

- `patch`:

  A [`logical`](https://rdrr.io/r/base/logical.html) if centroid should
  be calculated weighted by values.

#### Returns

Returns a [`sf`](https://r-spatial.github.io/sf/reference/sf.html)
object.

------------------------------------------------------------------------

### `BiodiversityScenario$save()`

Save object as output somewhere

#### Usage

    BiodiversityScenario$save(fname, type = "tif", dt = "FLT4S")

#### Arguments

- `fname`:

  An output filename as
  [`character`](https://rdrr.io/r/base/character.html).

- `type`:

  A format as [`character`](https://rdrr.io/r/base/character.html).
  Matched against a list of supported formats.

- `dt`:

  The datatype used, such as float64

#### Returns

Saved spatial prediction on drive.

------------------------------------------------------------------------

### `BiodiversityScenario$clone()`

The objects of this class are cloneable with this method.

#### Usage

    BiodiversityScenario$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
