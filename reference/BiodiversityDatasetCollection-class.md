# BiodiversityDatasetCollection super class description

Acts a container for a specified set of BiodiversityDataset contained
within. Functions are provided to summarize across the
BiodiversityDataset-class objects.

## Note

This can likely be beautified further.

## Public fields

- `data`:

  A [`list`](https://rdrr.io/r/base/list.html) of
  [`BiodiversityDataset`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDataset-class.md)
  objects.

- `name`:

  The default name of this collection as
  [`character`](https://rdrr.io/r/base/character.html).

## Methods

### Public methods

- [`BiodiversityDatasetCollection$new()`](#method-BiodiversityDatasetCollection-new)

- [`BiodiversityDatasetCollection$print()`](#method-BiodiversityDatasetCollection-print)

- [`BiodiversityDatasetCollection$show()`](#method-BiodiversityDatasetCollection-show)

- [`BiodiversityDatasetCollection$get_types()`](#method-BiodiversityDatasetCollection-get_types)

- [`BiodiversityDatasetCollection$get_names()`](#method-BiodiversityDatasetCollection-get_names)

- [`BiodiversityDatasetCollection$set_data()`](#method-BiodiversityDatasetCollection-set_data)

- [`BiodiversityDatasetCollection$get_data_object()`](#method-BiodiversityDatasetCollection-get_data_object)

- [`BiodiversityDatasetCollection$get_data()`](#method-BiodiversityDatasetCollection-get_data)

- [`BiodiversityDatasetCollection$get_coordinates()`](#method-BiodiversityDatasetCollection-get_coordinates)

- [`BiodiversityDatasetCollection$mask()`](#method-BiodiversityDatasetCollection-mask)

- [`BiodiversityDatasetCollection$rm_data()`](#method-BiodiversityDatasetCollection-rm_data)

- [`BiodiversityDatasetCollection$length()`](#method-BiodiversityDatasetCollection-length)

- [`BiodiversityDatasetCollection$get_observations()`](#method-BiodiversityDatasetCollection-get_observations)

- [`BiodiversityDatasetCollection$get_equations()`](#method-BiodiversityDatasetCollection-get_equations)

- [`BiodiversityDatasetCollection$get_families()`](#method-BiodiversityDatasetCollection-get_families)

- [`BiodiversityDatasetCollection$get_links()`](#method-BiodiversityDatasetCollection-get_links)

- [`BiodiversityDatasetCollection$get_columns_occ()`](#method-BiodiversityDatasetCollection-get_columns_occ)

- [`BiodiversityDatasetCollection$get_weights()`](#method-BiodiversityDatasetCollection-get_weights)

- [`BiodiversityDatasetCollection$get_ids()`](#method-BiodiversityDatasetCollection-get_ids)

- [`BiodiversityDatasetCollection$get_id_byType()`](#method-BiodiversityDatasetCollection-get_id_byType)

- [`BiodiversityDatasetCollection$get_id_byName()`](#method-BiodiversityDatasetCollection-get_id_byName)

- [`BiodiversityDatasetCollection$show_equations()`](#method-BiodiversityDatasetCollection-show_equations)

- [`BiodiversityDatasetCollection$plot()`](#method-BiodiversityDatasetCollection-plot)

- [`BiodiversityDatasetCollection$clone()`](#method-BiodiversityDatasetCollection-clone)

------------------------------------------------------------------------

### Method `new()`

Initializes the object and creates an empty list

#### Usage

    BiodiversityDatasetCollection$new()

#### Returns

NULL

------------------------------------------------------------------------

### Method [`print()`](https://iiasa.github.io/ibis.iSDM/reference/print.md)

Print the names and properties of all Biodiversity datasets contained
within

#### Usage

    BiodiversityDatasetCollection$print(format = TRUE)

#### Arguments

- `format`:

  A [`logical`](https://rdrr.io/r/base/logical.html) flag on whether a
  message should be printed.

#### Returns

A message on screen

------------------------------------------------------------------------

### Method `show()`

Aliases that calls print.

#### Usage

    BiodiversityDatasetCollection$show()

#### Returns

A message on screen

------------------------------------------------------------------------

### Method `get_types()`

Types of all biodiversity datasets included in this

#### Usage

    BiodiversityDatasetCollection$get_types(short = FALSE)

#### Arguments

- `short`:

  A [`logical`](https://rdrr.io/r/base/logical.html) flag whether types
  should be in short format.

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) vector.

------------------------------------------------------------------------

### Method `get_names()`

Get names and format them if necessary

#### Usage

    BiodiversityDatasetCollection$get_names(format = FALSE)

#### Arguments

- `format`:

  A [`logical`](https://rdrr.io/r/base/logical.html) flag whether names
  are to be formatted

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) vector.

------------------------------------------------------------------------

### Method `set_data()`

Add a new Biodiversity dataset to this collection.

#### Usage

    BiodiversityDatasetCollection$set_data(x, value)

#### Arguments

- `x`:

  A [`character`](https://rdrr.io/r/base/character.html) with the name
  or id of this dataset.

- `value`:

  A BiodiversityDataset

#### Returns

Invisible

------------------------------------------------------------------------

### Method `get_data_object()`

Get a specific Biodiversity dataset by id

#### Usage

    BiodiversityDatasetCollection$get_data_object(id)

#### Arguments

- `id`:

  A [`character`](https://rdrr.io/r/base/character.html) with a given id
  for the dataset.

#### Returns

Returns a BiodiversityDataset.

------------------------------------------------------------------------

### Method [`get_data()`](https://iiasa.github.io/ibis.iSDM/reference/get_data.md)

Get all biodiversity observations from a given dataset.

#### Usage

    BiodiversityDatasetCollection$get_data(id)

#### Arguments

- `id`:

  A [`character`](https://rdrr.io/r/base/character.html) with a given id
  for the dataset.

#### Returns

Returns all data from a set BiodiversityDataset.

------------------------------------------------------------------------

### Method `get_coordinates()`

Get coordinates for a given biodiversity dataset. Else return a wkt
object

#### Usage

    BiodiversityDatasetCollection$get_coordinates(id)

#### Arguments

- `id`:

  A [`character`](https://rdrr.io/r/base/character.html) with a given id
  for the dataset.

#### Returns

All coordinates from a given object in `data.frame`.

------------------------------------------------------------------------

### Method [`mask()`](https://iiasa.github.io/ibis.iSDM/reference/mask.md)

Convenience function to mask all input datasets.

#### Usage

    BiodiversityDatasetCollection$mask(mask, inverse = FALSE)

#### Arguments

- `mask`:

  A `SpatRaster` or `sf` object.

- `inverse`:

  A `logical` flag if the inverse should be masked instead.

#### Returns

Invisible

------------------------------------------------------------------------

### Method `rm_data()`

Remove a specific biodiversity dataset by id

#### Usage

    BiodiversityDatasetCollection$rm_data(id)

#### Arguments

- `id`:

  A [`character`](https://rdrr.io/r/base/character.html) with a given id
  for the dataset.

#### Returns

Invisible

------------------------------------------------------------------------

### Method [`length()`](https://rdrr.io/r/base/length.html)

Number of Biodiversity Datasets in connection

#### Usage

    BiodiversityDatasetCollection$length()

#### Returns

A [`numeric`](https://rdrr.io/r/base/numeric.html) with the number of
datasets.

------------------------------------------------------------------------

### Method `get_observations()`

Get number of observations of all datasets

#### Usage

    BiodiversityDatasetCollection$get_observations()

#### Returns

A [`numeric`](https://rdrr.io/r/base/numeric.html) with the number of
observations across datasets.

------------------------------------------------------------------------

### Method `get_equations()`

Get equations from all datasets

#### Usage

    BiodiversityDatasetCollection$get_equations()

#### Returns

A [`list`](https://rdrr.io/r/base/list.html) vector with all equations
across datasets.

------------------------------------------------------------------------

### Method `get_families()`

Get families from datasets.

#### Usage

    BiodiversityDatasetCollection$get_families()

#### Returns

A [`list`](https://rdrr.io/r/base/list.html) vector with all families
across datasets.

------------------------------------------------------------------------

### Method `get_links()`

Get custom link functions

#### Usage

    BiodiversityDatasetCollection$get_links()

#### Returns

A [`list`](https://rdrr.io/r/base/list.html) vector with all link
functions across datasets.

------------------------------------------------------------------------

### Method `get_columns_occ()`

Get fields with observation columns

#### Usage

    BiodiversityDatasetCollection$get_columns_occ()

#### Returns

A [`list`](https://rdrr.io/r/base/list.html) vector with the names of
observation columns.

------------------------------------------------------------------------

### Method `get_weights()`

Get the weights across datasets.

#### Usage

    BiodiversityDatasetCollection$get_weights()

#### Returns

A [`list`](https://rdrr.io/r/base/list.html) vector with the weights if
set per dataset.

------------------------------------------------------------------------

### Method `get_ids()`

Get ids of all assets in the collection.

#### Usage

    BiodiversityDatasetCollection$get_ids()

#### Returns

A [`list`](https://rdrr.io/r/base/list.html) vector with the ids of all
datasets.

------------------------------------------------------------------------

### Method `get_id_byType()`

Search for a specific biodiversity dataset with type

#### Usage

    BiodiversityDatasetCollection$get_id_byType(type)

#### Arguments

- `type`:

  A [`character`](https://rdrr.io/r/base/character.html) for a given
  data type.

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) with the id(s) of
datasets with the given type.

------------------------------------------------------------------------

### Method `get_id_byName()`

Get id by name

#### Usage

    BiodiversityDatasetCollection$get_id_byName(name)

#### Arguments

- `name`:

  A [`character`](https://rdrr.io/r/base/character.html) for a given
  name.

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) with the id(s) of
datasets with the given name.

------------------------------------------------------------------------

### Method `show_equations()`

Show equations of all datasets

#### Usage

    BiodiversityDatasetCollection$show_equations(msg = TRUE)

#### Arguments

- `msg`:

  A [`logical`](https://rdrr.io/r/base/logical.html) on whether to use
  print a message instead.

#### Returns

Shows equations on screen or as
[`character`](https://rdrr.io/r/base/character.html).

------------------------------------------------------------------------

### Method [`plot()`](https://iiasa.github.io/ibis.iSDM/reference/plot.md)

Plot the whole collection

#### Usage

    BiodiversityDatasetCollection$plot()

#### Returns

Invisible

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    BiodiversityDatasetCollection$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
