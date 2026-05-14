# PredictorDataset class description

This class describes the PredictorDataset and is used to store
covariates within.

## See also

[`predictor_derivate()`](https://iiasa.github.io/ibis.iSDM/reference/predictor_derivate.md)

[`predictor_transform()`](https://iiasa.github.io/ibis.iSDM/reference/predictor_transform.md)

## Public fields

- `id`:

  The id for this collection as
  [`character`](https://rdrr.io/r/base/character.html).

- `data`:

  A predictor dataset usually as
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html).

- `name`:

  A name for this object.

- `transformed`:

  Saves whether the predictors have been transformed somehow.

- `timeperiod`:

  A timeperiod field

- `is.spatial`:

  A [`logical`](https://rdrr.io/r/base/logical.html) flag on whether the
  predictor is spatial or not.

## Methods

### Public methods

- [`PredictorDataset$new()`](#method-PredictorDataset-initialize)

- [`PredictorDataset$print()`](#method-PredictorDataset-print)

- [`PredictorDataset$get_name()`](#method-PredictorDataset-get_name)

- [`PredictorDataset$get_id()`](#method-PredictorDataset-get_id)

- [`PredictorDataset$get_names()`](#method-PredictorDataset-get_names)

- [`PredictorDataset$get_predictor_names()`](#method-PredictorDataset-get_predictor_names)

- [`PredictorDataset$get_data()`](#method-PredictorDataset-get_data)

- [`PredictorDataset$get_time()`](#method-PredictorDataset-get_time)

- [`PredictorDataset$get_projection()`](#method-PredictorDataset-get_projection)

- [`PredictorDataset$get_resolution()`](#method-PredictorDataset-get_resolution)

- [`PredictorDataset$get_ext()`](#method-PredictorDataset-get_ext)

- [`PredictorDataset$crop_data()`](#method-PredictorDataset-crop_data)

- [`PredictorDataset$mask()`](#method-PredictorDataset-mask)

- [`PredictorDataset$set_data()`](#method-PredictorDataset-set_data)

- [`PredictorDataset$rm_data()`](#method-PredictorDataset-rm_data)

- [`PredictorDataset$show()`](#method-PredictorDataset-show)

- [`PredictorDataset$summary()`](#method-PredictorDataset-summary)

- [`PredictorDataset$has_derivates()`](#method-PredictorDataset-has_derivates)

- [`PredictorDataset$is_transformed()`](#method-PredictorDataset-is_transformed)

- [`PredictorDataset$is_spatial()`](#method-PredictorDataset-is_spatial)

- [`PredictorDataset$get_transformed_params()`](#method-PredictorDataset-get_transformed_params)

- [`PredictorDataset$length()`](#method-PredictorDataset-length)

- [`PredictorDataset$ncell()`](#method-PredictorDataset-ncell)

- [`PredictorDataset$plot()`](#method-PredictorDataset-plot)

- [`PredictorDataset$clone()`](#method-PredictorDataset-clone)

------------------------------------------------------------------------

### `PredictorDataset$new()`

Initializes the object and creates an empty list

#### Usage

    PredictorDataset$new(id, data, transformed = FALSE, ...)

#### Arguments

- `id`:

  The id for this collection as
  [`character`](https://rdrr.io/r/base/character.html).

- `data`:

  A predictor dataset usually as
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html).

- `transformed`:

  A [`logical`](https://rdrr.io/r/base/logical.html) flag if predictors
  have been transformed. Assume not.

- `...`:

  Any other parameters found.

#### Returns

NULL

------------------------------------------------------------------------

### `PredictorDataset$print()`

Print the names and properties of all Biodiversity datasets contained
within

#### Usage

    PredictorDataset$print(format = TRUE)

#### Arguments

- `format`:

  A [`logical`](https://rdrr.io/r/base/logical.html) flag on whether a
  message should be printed.

#### Returns

A message on screen

------------------------------------------------------------------------

### `PredictorDataset$get_name()`

Return name of this object

#### Usage

    PredictorDataset$get_name()

#### Returns

Default [`character`](https://rdrr.io/r/base/character.html) name.

------------------------------------------------------------------------

### `PredictorDataset$get_id()`

Get Id of this object

#### Usage

    PredictorDataset$get_id()

#### Returns

Default [`character`](https://rdrr.io/r/base/character.html) name.

------------------------------------------------------------------------

### `PredictorDataset$get_names()`

Get names of data

#### Usage

    PredictorDataset$get_names()

#### Returns

[`character`](https://rdrr.io/r/base/character.html) names of the data
value.

------------------------------------------------------------------------

### `PredictorDataset$get_predictor_names()`

Alias for get_names

#### Usage

    PredictorDataset$get_predictor_names()

#### Returns

[`character`](https://rdrr.io/r/base/character.html) names of the data
value.

------------------------------------------------------------------------

### `PredictorDataset$get_data()`

Get a specific dataset

#### Usage

    PredictorDataset$get_data(df = FALSE, na.rm = TRUE, ...)

#### Arguments

- `df`:

  [`logical`](https://rdrr.io/r/base/logical.html) on whether data is to
  be returned as [`data.frame`](https://rdrr.io/r/base/data.frame.html).

- `na.rm`:

  [`logical`](https://rdrr.io/r/base/logical.html) if `NA` is to be
  removed from data.frame.

- `...`:

  Any other parameters passed on.

#### Returns

A
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
or [`data.frame`](https://rdrr.io/r/base/data.frame.html).

------------------------------------------------------------------------

### `PredictorDataset$get_time()`

Get time dimension of object.

#### Usage

    PredictorDataset$get_time(...)

#### Arguments

- `...`:

  Any other parameters passed on.

#### Returns

A [`vector`](https://rdrr.io/r/base/vector.html) with the time dimension
of the dataset.

------------------------------------------------------------------------

### `PredictorDataset$get_projection()`

Get Projection

#### Usage

    PredictorDataset$get_projection()

#### Returns

A [`vector`](https://rdrr.io/r/base/vector.html) with the geographical
projection of the object.

------------------------------------------------------------------------

### `PredictorDataset$get_resolution()`

Get Resolution

#### Usage

    PredictorDataset$get_resolution()

#### Returns

A [`numeric`](https://rdrr.io/r/base/numeric.html)
[`vector`](https://rdrr.io/r/base/vector.html) with the spatial
resolution of the data.

------------------------------------------------------------------------

### `PredictorDataset$get_ext()`

Get Extent of predictors

#### Usage

    PredictorDataset$get_ext()

#### Returns

A [`numeric`](https://rdrr.io/r/base/numeric.html)
[`vector`](https://rdrr.io/r/base/vector.html) with the spatial
resolution of the data.

------------------------------------------------------------------------

### `PredictorDataset$crop_data()`

Utility function to clip the predictor dataset by another dataset

#### Usage

    PredictorDataset$crop_data(pol, apply_time = FALSE)

#### Arguments

- `pol`:

  A [`sf`](https://r-spatial.github.io/sf/reference/sf.html) object used
  for cropping the data.

- `apply_time`:

  A [`logical`](https://rdrr.io/r/base/logical.html) flag indicating if
  time should be acknowledged in cropping.

#### Details

This code now also is able to determine the temporally closest layer. In
case a [`data.frame`](https://rdrr.io/r/base/data.frame.html) exists as
predictor, only that is returned.

#### Returns

Invisible TRUE

------------------------------------------------------------------------

### `PredictorDataset$mask()`

Utility function to mask the predictor dataset by another dataset

#### Usage

    PredictorDataset$mask(mask, inverse = FALSE, ...)

#### Arguments

- `mask`:

  A `SpatRaster` or `sf` object.

- `inverse`:

  A `logical` flag if the inverse should be masked instead.

- `...`:

  Any other parameters passed on to masking.

#### Returns

Invisible

------------------------------------------------------------------------

### `PredictorDataset$set_data()`

Add a new Predictor dataset to this collection

#### Usage

    PredictorDataset$set_data(value)

#### Arguments

- `value`:

  A new
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  or [`stars`](https://rdrr.io/r/graphics/stars.html) object.

#### Returns

This object

------------------------------------------------------------------------

### `PredictorDataset$rm_data()`

Remove a specific Predictor by name

#### Usage

    PredictorDataset$rm_data(x)

#### Arguments

- `x`:

  [`character`](https://rdrr.io/r/base/character.html) of the predictor
  name to be removed.

#### Returns

Invisible

------------------------------------------------------------------------

### `PredictorDataset$show()`

Alias for print method

#### Usage

    PredictorDataset$show()

#### Returns

Invisible

------------------------------------------------------------------------

### `PredictorDataset$summary()`

Collect info statistics with optional decimals

#### Usage

    PredictorDataset$summary(digits = 2)

#### Arguments

- `digits`:

  [`numeric`](https://rdrr.io/r/base/numeric.html) Giving the rounding
  precision

#### Returns

A [`data.frame`](https://rdrr.io/r/base/data.frame.html) summarizing the
data.

------------------------------------------------------------------------

### `PredictorDataset$has_derivates()`

Indication if there are any predictors that are derivates of outers

#### Usage

    PredictorDataset$has_derivates()

#### Returns

A [`logical`](https://rdrr.io/r/base/logical.html) flag.

------------------------------------------------------------------------

### `PredictorDataset$is_transformed()`

Predictors have been transformed?

#### Usage

    PredictorDataset$is_transformed()

#### Returns

A [`logical`](https://rdrr.io/r/base/logical.html) flag.

------------------------------------------------------------------------

### `PredictorDataset$is_spatial()`

Is Predictor dataset spatial?

#### Usage

    PredictorDataset$is_spatial()

#### Returns

A [`logical`](https://rdrr.io/r/base/logical.html) flag.

------------------------------------------------------------------------

### `PredictorDataset$get_transformed_params()`

Get transformation params.

#### Usage

    PredictorDataset$get_transformed_params()

#### Returns

A [`matrix`](https://rdrr.io/r/base/matrix.html) flag.

------------------------------------------------------------------------

### `PredictorDataset$length()`

Number of Predictors in object

#### Usage

    PredictorDataset$length()

#### Returns

A [`numeric`](https://rdrr.io/r/base/numeric.html) estimate

------------------------------------------------------------------------

### `PredictorDataset$ncell()`

Number of cells or values in object

#### Usage

    PredictorDataset$ncell()

#### Returns

A [`numeric`](https://rdrr.io/r/base/numeric.html) estimate

------------------------------------------------------------------------

### `PredictorDataset$plot()`

Basic Plotting function

#### Usage

    PredictorDataset$plot()

#### Returns

A graphical interpretation of the predictors in this object.

------------------------------------------------------------------------

### `PredictorDataset$clone()`

The objects of this class are cloneable with this method.

#### Usage

    PredictorDataset$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
