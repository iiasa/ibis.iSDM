# BiodiversityDataset prototype description

BiodiversityDataset prototype description

## Public fields

- `name`:

  The default name of this dataset as
  [`character`](https://rdrr.io/r/base/character.html).

- `id`:

  A [`character`](https://rdrr.io/r/base/character.html) with the unique
  id for this dataset.

- `equation`:

  A [`formula`](https://rdrr.io/r/stats/formula.html) object containing
  the equation of how this dataset is modelled.

- `family`:

  The family used for this dataset as
  [`character`](https://rdrr.io/r/base/character.html).

- `link`:

  The link function used for this data as
  [`character`](https://rdrr.io/r/base/character.html).

- `type`:

  A [`character`](https://rdrr.io/r/base/character.html) with the type
  as [`character`](https://rdrr.io/r/base/character.html).

- `weight`:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) containing custom
  weights per observation for this dataset.

- `field_occurrence`:

  A [`character`](https://rdrr.io/r/base/character.html) with the name
  of the column name containing observations.

- `data`:

  Contains the observational data in
  [`sf`](https://r-spatial.github.io/sf/reference/sf.html) format.

- `use_intercept`:

  A [`logical`](https://rdrr.io/r/base/logical.html) flag on whether
  intercepts are included for this dataset.

- `pseudoabsence_settings`:

  Optionally provided pseudoabsence settings.

## Methods

### Public methods

- [`BiodiversityDataset$new()`](#method-BiodiversityDataset-initialize)

- [`BiodiversityDataset$print()`](#method-BiodiversityDataset-print)

- [`BiodiversityDataset$set_equation()`](#method-BiodiversityDataset-set_equation)

- [`BiodiversityDataset$get_equation()`](#method-BiodiversityDataset-get_equation)

- [`BiodiversityDataset$show_equation()`](#method-BiodiversityDataset-show_equation)

- [`BiodiversityDataset$get_id()`](#method-BiodiversityDataset-get_id)

- [`BiodiversityDataset$get_type()`](#method-BiodiversityDataset-get_type)

- [`BiodiversityDataset$get_column_occ()`](#method-BiodiversityDataset-get_column_occ)

- [`BiodiversityDataset$get_family()`](#method-BiodiversityDataset-get_family)

- [`BiodiversityDataset$get_link()`](#method-BiodiversityDataset-get_link)

- [`BiodiversityDataset$get_data()`](#method-BiodiversityDataset-get_data)

- [`BiodiversityDataset$get_weight()`](#method-BiodiversityDataset-get_weight)

- [`BiodiversityDataset$show()`](#method-BiodiversityDataset-show)

- [`BiodiversityDataset$get_observations()`](#method-BiodiversityDataset-get_observations)

- [`BiodiversityDataset$mask()`](#method-BiodiversityDataset-mask)

- [`BiodiversityDataset$clone()`](#method-BiodiversityDataset-clone)

------------------------------------------------------------------------

### `BiodiversityDataset$new()`

Initializes the object and creates an empty list

#### Usage

    BiodiversityDataset$new(
      name,
      id,
      equation,
      family,
      link,
      type,
      weight,
      field_occurrence,
      data,
      use_intercept,
      pseudoabsence_settings
    )

#### Arguments

- `name`:

  The default name of this dataset as
  [`character`](https://rdrr.io/r/base/character.html).

- `id`:

  A [`character`](https://rdrr.io/r/base/character.html) with the unique
  id for this dataset.

- `equation`:

  A [`formula`](https://rdrr.io/r/stats/formula.html) object containing
  the equation of how this dataset is modelled.

- `family`:

  The family used for this dataset as
  [`character`](https://rdrr.io/r/base/character.html).

- `link`:

  The link function used for this data as
  [`character`](https://rdrr.io/r/base/character.html).

- `type`:

  A [`character`](https://rdrr.io/r/base/character.html) with the type
  as [`character`](https://rdrr.io/r/base/character.html).

- `weight`:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) containing custom
  weights per observation for this dataset.

- `field_occurrence`:

  A [`character`](https://rdrr.io/r/base/character.html) with the name
  of the column name containing observations.

- `data`:

  Contains the observational data in
  [`sf`](https://r-spatial.github.io/sf/reference/sf.html) format.

- `use_intercept`:

  A [`logical`](https://rdrr.io/r/base/logical.html) flag on whether
  intercepts are included for this dataset.

- `pseudoabsence_settings`:

  Optionally provided pseudoabsence settings.

#### Returns

NULL

------------------------------------------------------------------------

### `BiodiversityDataset$print()`

Print the names and properties of all Biodiversity datasets contained
within

#### Usage

    BiodiversityDataset$print()

#### Returns

A message on screen

------------------------------------------------------------------------

### `BiodiversityDataset$set_equation()`

Set new equation and writes it into `formula`

#### Usage

    BiodiversityDataset$set_equation(x)

#### Arguments

- `x`:

  A new [`formula`](https://rdrr.io/r/stats/formula.html) object.

#### Returns

Invisible

------------------------------------------------------------------------

### `BiodiversityDataset$get_equation()`

Get equation

#### Usage

    BiodiversityDataset$get_equation()

#### Returns

A placeholder or [`formula`](https://rdrr.io/r/stats/formula.html)
object.

------------------------------------------------------------------------

### `BiodiversityDataset$show_equation()`

Function to print the equation

#### Usage

    BiodiversityDataset$show_equation()

#### Returns

A message on screen.

------------------------------------------------------------------------

### `BiodiversityDataset$get_id()`

Get Id within the dataset

#### Usage

    BiodiversityDataset$get_id()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) with the id.

------------------------------------------------------------------------

### `BiodiversityDataset$get_type()`

Get type of the dataset.

#### Usage

    BiodiversityDataset$get_type(short = FALSE)

#### Arguments

- `short`:

  A [`logical`](https://rdrr.io/r/base/logical.html) flag if this should
  be formatted in shortform.

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) with the type

------------------------------------------------------------------------

### `BiodiversityDataset$get_column_occ()`

Get field with occurrence information

#### Usage

    BiodiversityDataset$get_column_occ()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) with the
occurence field

------------------------------------------------------------------------

### `BiodiversityDataset$get_family()`

Get family

#### Usage

    BiodiversityDataset$get_family()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) with the family
for the dataset

------------------------------------------------------------------------

### `BiodiversityDataset$get_link()`

Get custom link function

#### Usage

    BiodiversityDataset$get_link()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) with the family
for the dataset

------------------------------------------------------------------------

### `BiodiversityDataset$get_data()`

Get data from the object

#### Usage

    BiodiversityDataset$get_data()

#### Returns

A [`sf`](https://r-spatial.github.io/sf/reference/sf.html) object with
the data

------------------------------------------------------------------------

### `BiodiversityDataset$get_weight()`

Get weight

#### Usage

    BiodiversityDataset$get_weight()

#### Returns

A [`numeric`](https://rdrr.io/r/base/numeric.html) with the weights
within the dataset.

------------------------------------------------------------------------

### `BiodiversityDataset$show()`

Print input messages

#### Usage

    BiodiversityDataset$show()

#### Returns

A message on screen.

------------------------------------------------------------------------

### `BiodiversityDataset$get_observations()`

Collect info statistics about number of observations

#### Usage

    BiodiversityDataset$get_observations()

#### Returns

A [`numeric`](https://rdrr.io/r/base/numeric.html) with the number of
observations.

------------------------------------------------------------------------

### `BiodiversityDataset$mask()`

Convenience function to mask all input datasets.

#### Usage

    BiodiversityDataset$mask(mask, inverse = FALSE, ...)

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

### `BiodiversityDataset$clone()`

The objects of this class are cloneable with this method.

#### Usage

    BiodiversityDataset$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
