# Class for the trained Model object

All trained Models inherit the options here plus any additional ones
defined by the engine and inference.

## Note

Could be further pretified and commands outsourced.

## Public fields

- `id`:

  A character id for any trained model

- `name`:

  A description of the model as
  [`character`](https://rdrr.io/r/base/character.html).

- `model`:

  A [`list`](https://rdrr.io/r/base/list.html) containing all input
  datasets and parameters to the model.

- `settings`:

  A
  [`Settings`](https://iiasa.github.io/ibis.iSDM/reference/Settings-class.md)
  object with information on inference.

- `fits`:

  A [`list`](https://rdrr.io/r/base/list.html) containing the prediction
  and fitted model.

- `.internals`:

  A [`list`](https://rdrr.io/r/base/list.html) containing previous
  fitted models.

## Methods

### Public methods

- [`DistributionModel$new()`](#method-DistributionModel-new)

- [`DistributionModel$get_name()`](#method-DistributionModel-get_name)

- [`DistributionModel$print()`](#method-DistributionModel-print)

- [`DistributionModel$show()`](#method-DistributionModel-show)

- [`DistributionModel$plot()`](#method-DistributionModel-plot)

- [`DistributionModel$plot_threshold()`](#method-DistributionModel-plot_threshold)

- [`DistributionModel$show_duration()`](#method-DistributionModel-show_duration)

- [`DistributionModel$summary()`](#method-DistributionModel-summary)

- [`DistributionModel$effects()`](#method-DistributionModel-effects)

- [`DistributionModel$get_equation()`](#method-DistributionModel-get_equation)

- [`DistributionModel$get_data()`](#method-DistributionModel-get_data)

- [`DistributionModel$get_model()`](#method-DistributionModel-get_model)

- [`DistributionModel$set_data()`](#method-DistributionModel-set_data)

- [`DistributionModel$get_thresholdvalue()`](#method-DistributionModel-get_thresholdvalue)

- [`DistributionModel$get_thresholdtype()`](#method-DistributionModel-get_thresholdtype)

- [`DistributionModel$show_rasters()`](#method-DistributionModel-show_rasters)

- [`DistributionModel$get_projection()`](#method-DistributionModel-get_projection)

- [`DistributionModel$get_resolution()`](#method-DistributionModel-get_resolution)

- [`DistributionModel$rm_threshold()`](#method-DistributionModel-rm_threshold)

- [`DistributionModel$calc_suitabilityindex()`](#method-DistributionModel-calc_suitabilityindex)

- [`DistributionModel$get_centroid()`](#method-DistributionModel-get_centroid)

- [`DistributionModel$has_limits()`](#method-DistributionModel-has_limits)

- [`DistributionModel$has_latent()`](#method-DistributionModel-has_latent)

- [`DistributionModel$has_offset()`](#method-DistributionModel-has_offset)

- [`DistributionModel$mask()`](#method-DistributionModel-mask)

- [`DistributionModel$save()`](#method-DistributionModel-save)

- [`DistributionModel$clone()`](#method-DistributionModel-clone)

------------------------------------------------------------------------

### Method `new()`

Initializes the object and creates an empty list

#### Usage

    DistributionModel$new(name)

#### Arguments

- `name`:

  A description of the model as
  [`character`](https://rdrr.io/r/base/character.html).

#### Returns

NULL

------------------------------------------------------------------------

### Method `get_name()`

Return the name of the model

#### Usage

    DistributionModel$get_name()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) with the model
name used.

------------------------------------------------------------------------

### Method [`print()`](https://iiasa.github.io/ibis.iSDM/reference/print.md)

Print the names and summarizes the model within

#### Usage

    DistributionModel$print()

#### Returns

A message on screen

------------------------------------------------------------------------

### Method `show()`

Show the name of the Model.

#### Usage

    DistributionModel$show()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) of the run name.

------------------------------------------------------------------------

### Method [`plot()`](https://iiasa.github.io/ibis.iSDM/reference/plot.md)

Plots the prediction if found.

#### Usage

    DistributionModel$plot(what = "mean")

#### Arguments

- `what`:

  [`character`](https://rdrr.io/r/base/character.html) with the specific
  layer to be plotted.

#### Returns

A graphical representation of the prediction

------------------------------------------------------------------------

### Method `plot_threshold()`

Plots the thresholded prediction if found.

#### Usage

    DistributionModel$plot_threshold(what = 1)

#### Arguments

- `what`:

  [`character`](https://rdrr.io/r/base/character.html) or
  [`numeric`](https://rdrr.io/r/base/numeric.html) for the layer to be
  plotted.

#### Returns

A graphical representation of the thresholded prediction if found.

------------------------------------------------------------------------

### Method `show_duration()`

Show model run time if settings exist

#### Usage

    DistributionModel$show_duration()

#### Returns

A [`numeric`](https://rdrr.io/r/base/numeric.html) estimate of the
duration it took to fit the models.

------------------------------------------------------------------------

### Method [`summary()`](https://iiasa.github.io/ibis.iSDM/reference/summary.md)

Get effects or importance tables from model

#### Usage

    DistributionModel$summary(obj = "fit_best")

#### Arguments

- `obj`:

  A [`character`](https://rdrr.io/r/base/character.html) of which object
  to return.

#### Returns

A [`data.frame`](https://rdrr.io/r/base/data.frame.html) summarizing the
model, usually its coefficient.

------------------------------------------------------------------------

### Method [`effects()`](https://iiasa.github.io/ibis.iSDM/reference/effects.md)

Generic plotting function for effect plots

#### Usage

    DistributionModel$effects(x = "fit_best", what = "fixed", ...)

#### Arguments

- `x`:

  A [`character`](https://rdrr.io/r/base/character.html) for the object
  in question.

- `what`:

  A [`character`](https://rdrr.io/r/base/character.html) for the type of
  coefficients.

- `...`:

  Any other options.

#### Returns

A graphical representation of the coefficents.

------------------------------------------------------------------------

### Method `get_equation()`

Get equation

#### Usage

    DistributionModel$get_equation()

#### Returns

A [`formula`](https://rdrr.io/r/stats/formula.html) of the inferred
model.

------------------------------------------------------------------------

### Method [`get_data()`](https://iiasa.github.io/ibis.iSDM/reference/get_data.md)

Get specific fit from this Model

#### Usage

    DistributionModel$get_data(x = "prediction")

#### Arguments

- `x`:

  A [`character`](https://rdrr.io/r/base/character.html) stating what
  should be returned.

#### Returns

A
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
object with the prediction.

------------------------------------------------------------------------

### Method `get_model()`

Small internal helper function to directly get the model object

#### Usage

    DistributionModel$get_model()

#### Returns

A fitted model if existing.

------------------------------------------------------------------------

### Method `set_data()`

Set new fit for this Model.

#### Usage

    DistributionModel$set_data(x, value)

#### Arguments

- `x`:

  The name of the new fit.

- `value`:

  The
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  layer (or model) to be inserted.

#### Returns

This object.

------------------------------------------------------------------------

### Method `get_thresholdvalue()`

Get the threshold value if calculated

#### Usage

    DistributionModel$get_thresholdvalue()

#### Returns

A [`numeric`](https://rdrr.io/r/base/numeric.html) threshold value.

------------------------------------------------------------------------

### Method `get_thresholdtype()`

Get threshold type and format if calculated.

#### Usage

    DistributionModel$get_thresholdtype()

#### Returns

A vector with a [`character`](https://rdrr.io/r/base/character.html)
method and [`numeric`](https://rdrr.io/r/base/numeric.html) threshold
value.

------------------------------------------------------------------------

### Method `show_rasters()`

List all rasters in object

#### Usage

    DistributionModel$show_rasters()

#### Returns

A [`vector`](https://rdrr.io/r/base/vector.html) with
[`logical`](https://rdrr.io/r/base/logical.html) flags for the various
objects.

------------------------------------------------------------------------

### Method `get_projection()`

Get projection of the background.

#### Usage

    DistributionModel$get_projection()

#### Returns

A geographic projection

------------------------------------------------------------------------

### Method `get_resolution()`

Get the resolution of the projection

#### Usage

    DistributionModel$get_resolution()

#### Returns

[`numeric`](https://rdrr.io/r/base/numeric.html) estimates of the
distribution.

------------------------------------------------------------------------

### Method `rm_threshold()`

Remove calculated thresholds

#### Usage

    DistributionModel$rm_threshold()

#### Returns

Invisible

------------------------------------------------------------------------

### Method `calc_suitabilityindex()`

Calculate a suitability index for a given projection

#### Usage

    DistributionModel$calc_suitabilityindex(method = "normalize")

#### Arguments

- `method`:

  The method used for normalization.

#### Details

Methods can either be normalized by the minimum and maximum. Or the
relative total using the sumof values.

#### Returns

Returns a
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html).

------------------------------------------------------------------------

### Method `get_centroid()`

Get centroids of prediction layers

#### Usage

    DistributionModel$get_centroid(patch = FALSE, layer = "mean")

#### Arguments

- `patch`:

  A [`logical`](https://rdrr.io/r/base/logical.html) if centroid should
  be calculated weighted by values.

- `layer`:

  [`character`](https://rdrr.io/r/base/character.html) of the layer to
  use.

#### Returns

Returns a [`sf::sf`](https://r-spatial.github.io/sf/reference/sf.html)
object.

------------------------------------------------------------------------

### Method `has_limits()`

Logical indication if the prediction was limited.

#### Usage

    DistributionModel$has_limits()

#### Returns

A [`logical`](https://rdrr.io/r/base/logical.html) flag.

------------------------------------------------------------------------

### Method `has_latent()`

Logical indication if the prediction has added latent factors.

#### Usage

    DistributionModel$has_latent()

#### Returns

A [`logical`](https://rdrr.io/r/base/logical.html) flag.

------------------------------------------------------------------------

### Method `has_offset()`

Has a offset been used?

#### Usage

    DistributionModel$has_offset()

#### Returns

A [`logical`](https://rdrr.io/r/base/logical.html) flag.

------------------------------------------------------------------------

### Method [`mask()`](https://iiasa.github.io/ibis.iSDM/reference/mask.md)

Convenience function to mask all input datasets.

#### Usage

    DistributionModel$mask(mask, inverse = FALSE, ...)

#### Arguments

- `mask`:

  A `SpatRaster` or `sf` object.

- `inverse`:

  A `logical` flag if the inverse should be masked instead.

- `...`:

  Any other parameters passed on to mask

#### Returns

Invisible

------------------------------------------------------------------------

### Method [`save()`](https://rdrr.io/r/base/save.html)

Save the prediction as output.

#### Usage

    DistributionModel$save(fname, type = "gtif", dt = "FLT4S")

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

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    DistributionModel$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
