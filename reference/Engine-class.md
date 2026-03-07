# Engine class description

Basic object for engine, all other engines inherit from here.

## Public fields

- `engine`:

  The class name of the engine.

- `name`:

  The name of the engine

- `data`:

  Any data or parameters necessary to make this engine work.

## Methods

### Public methods

- [`Engine$new()`](#method-Engine-new)

- [`Engine$print()`](#method-Engine-print)

- [`Engine$show()`](#method-Engine-show)

- [`Engine$get_class()`](#method-Engine-get_class)

- [`Engine$get_data()`](#method-Engine-get_data)

- [`Engine$list_data()`](#method-Engine-list_data)

- [`Engine$set_data()`](#method-Engine-set_data)

- [`Engine$get_self()`](#method-Engine-get_self)

- [`Engine$clone()`](#method-Engine-clone)

------------------------------------------------------------------------

### Method `new()`

Initializes the object and creates an empty list

#### Usage

    Engine$new(engine, name)

#### Arguments

- `engine`:

  The class name of the engine.

- `name`:

  The name of the engine

#### Returns

NULL

------------------------------------------------------------------------

### Method [`print()`](https://iiasa.github.io/ibis.iSDM/reference/print.md)

Print the Engine name

#### Usage

    Engine$print()

#### Returns

A message on screen

------------------------------------------------------------------------

### Method `show()`

Aliases that calls print.

#### Usage

    Engine$show()

#### Returns

A message on screen

------------------------------------------------------------------------

### Method `get_class()`

Get class description

#### Usage

    Engine$get_class()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) with the class as
saved in engine

------------------------------------------------------------------------

### Method [`get_data()`](https://iiasa.github.io/ibis.iSDM/reference/get_data.md)

Get specific data from this engine

#### Usage

    Engine$get_data(x)

#### Arguments

- `x`:

  A respecified data to be added to the engine.

#### Returns

A [`list`](https://rdrr.io/r/base/list.html) with the data.

------------------------------------------------------------------------

### Method `list_data()`

List all data

#### Usage

    Engine$list_data()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) vector of the
data entries.

------------------------------------------------------------------------

### Method `set_data()`

Set data for this engine

#### Usage

    Engine$set_data(x, value)

#### Arguments

- `x`:

  A [`character`](https://rdrr.io/r/base/character.html) with the name
  or id of this dataset.

- `value`:

  A new [`list`](https://rdrr.io/r/base/list.html) of parameters.

#### Returns

Invisible

------------------------------------------------------------------------

### Method `get_self()`

Dummy function to get self object

#### Usage

    Engine$get_self()

#### Returns

This object

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    Engine$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
