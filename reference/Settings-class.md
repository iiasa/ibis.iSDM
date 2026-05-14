# Prototype for model settings object

Basic [`R6`](https://r6.r-lib.org/reference/R6Class.html) object for
Settings object, a List that stores settings used related to model
training.

## Public fields

- `name`:

  The default name of this settings as
  [`character`](https://rdrr.io/r/base/character.html).

- `modelid`:

  A [`character`](https://rdrr.io/r/base/character.html) of the model id
  this belongs to.

- `data`:

  A [`list`](https://rdrr.io/r/base/list.html) of contained settings.

## Methods

### Public methods

- [`Settings$new()`](#method-Settings-initialize)

- [`Settings$print()`](#method-Settings-print)

- [`Settings$show()`](#method-Settings-show)

- [`Settings$length()`](#method-Settings-length)

- [`Settings$duration()`](#method-Settings-duration)

- [`Settings$summary()`](#method-Settings-summary)

- [`Settings$get()`](#method-Settings-get)

- [`Settings$set()`](#method-Settings-set)

- [`Settings$clone()`](#method-Settings-clone)

------------------------------------------------------------------------

### `Settings$new()`

Initializes the object and creates an empty list

#### Usage

    Settings$new()

#### Returns

NULL

------------------------------------------------------------------------

### `Settings$print()`

Print the names and properties of all Biodiversity datasets contained
within

#### Usage

    Settings$print()

#### Returns

A message on screen

------------------------------------------------------------------------

### `Settings$show()`

Shows the name and the settings

#### Usage

    Settings$show()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) of the name and
settings.

------------------------------------------------------------------------

### `Settings$length()`

Number of options

#### Usage

    Settings$length()

#### Returns

A [`numeric`](https://rdrr.io/r/base/numeric.html) with the number of
options.

------------------------------------------------------------------------

### `Settings$duration()`

Computation duration convenience function

#### Usage

    Settings$duration()

#### Returns

The amount of time passed for model fitting if found.

------------------------------------------------------------------------

### `Settings$summary()`

Summary call of the contained parameters

#### Usage

    Settings$summary()

#### Returns

A [`list`](https://rdrr.io/r/base/list.html) with the parameters in this
object.

------------------------------------------------------------------------

### `Settings$get()`

Get a specific setting

#### Usage

    Settings$get(what)

#### Arguments

- `what`:

  A [`character`](https://rdrr.io/r/base/character.html) with the
  respective setting.

#### Returns

The setting if found in the object.

------------------------------------------------------------------------

### `Settings$set()`

Set new settings

#### Usage

    Settings$set(what, x, copy = FALSE)

#### Arguments

- `what`:

  A [`character`](https://rdrr.io/r/base/character.html) with the name
  for the new settings.

- `x`:

  The new setting to be stored. Can be any object.

- `copy`:

  [`logical`](https://rdrr.io/r/base/logical.html) on whether a new
  settings object is to be created.

#### Returns

The setting if found in the object.

------------------------------------------------------------------------

### `Settings$clone()`

The objects of this class are cloneable with this method.

#### Usage

    Settings$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
