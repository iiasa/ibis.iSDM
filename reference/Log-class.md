# Log prototype.

Basic [`R6::R6`](https://r6.r-lib.org/reference/R6Class.html) object for
Log, any Log inherit from here

## Public fields

- `filename`:

  A [`character`](https://rdrr.io/r/base/character.html) of where the
  log is to be stored.

- `output`:

  The log content.

## Methods

### Public methods

- [`Log$new()`](#method-Log-new)

- [`Log$print()`](#method-Log-print)

- [`Log$open()`](#method-Log-open)

- [`Log$close()`](#method-Log-close)

- [`Log$get_filename()`](#method-Log-get_filename)

- [`Log$set_filename()`](#method-Log-set_filename)

- [`Log$delete()`](#method-Log-delete)

- [`Log$open_system()`](#method-Log-open_system)

- [`Log$clone()`](#method-Log-clone)

------------------------------------------------------------------------

### Method `new()`

Initializes the object and specifies some default parameters.

#### Usage

    Log$new(filename, output)

#### Arguments

- `filename`:

  A [`character`](https://rdrr.io/r/base/character.html) of where the
  log is to be stored.

- `output`:

  The log content.

#### Returns

NULL

------------------------------------------------------------------------

### Method [`print()`](https://iiasa.github.io/ibis.iSDM/reference/print.md)

Print message with filename

#### Usage

    Log$print()

#### Returns

A message on screen

------------------------------------------------------------------------

### Method [`open()`](https://rdrr.io/r/base/connections.html)

Opens the connection to the output filename.

#### Usage

    Log$open(type = c("output", "message"))

#### Arguments

- `type`:

  A [`character`](https://rdrr.io/r/base/character.html) vector of the
  output types.

#### Returns

Invisible TRUE

------------------------------------------------------------------------

### Method [`close()`](https://rdrr.io/r/base/connections.html)

Closes the connection to the output file

#### Usage

    Log$close()

#### Returns

Invisible TRUE

------------------------------------------------------------------------

### Method `get_filename()`

Get output filename

#### Usage

    Log$get_filename()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) with the filename

------------------------------------------------------------------------

### Method `set_filename()`

Set a new output filename

#### Usage

    Log$set_filename(value)

#### Arguments

- `value`:

  A [`character`](https://rdrr.io/r/base/character.html) with the new
  filename.

#### Returns

Invisible TRUE

------------------------------------------------------------------------

### Method `delete()`

Delete log file

#### Usage

    Log$delete()

#### Returns

Invisible TRUE

------------------------------------------------------------------------

### Method `open_system()`

Open log with system viewer

#### Usage

    Log$open_system()

#### Returns

Invisible TRUE

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    Log$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
