# Log prototype.

Basic [`R6`](https://r6.r-lib.org/reference/R6Class.html) object for
Log, any Log inherit from here

## Value

An [`R6::R6Class`](https://r6.r-lib.org/reference/R6Class.html)
generator object.

## Public fields

- `filename`:

  A [`character`](https://rdrr.io/r/base/character.html) of where the
  log is to be stored.

- `output`:

  The log content.

## Methods

### Public methods

- [`Log$new()`](#method-Log-initialize)

- [`Log$print()`](#method-Log-print)

- [`Log$open()`](#method-Log-open)

- [`Log$close()`](#method-Log-close)

- [`Log$get_filename()`](#method-Log-get_filename)

- [`Log$set_filename()`](#method-Log-set_filename)

- [`Log$delete()`](#method-Log-delete)

- [`Log$open_system()`](#method-Log-open_system)

- [`Log$clone()`](#method-Log-clone)

------------------------------------------------------------------------

### `Log$new()`

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

### `Log$print()`

Print message with filename

#### Usage

    Log$print()

#### Returns

A message on screen

------------------------------------------------------------------------

### `Log$open()`

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

### `Log$close()`

Closes the connection to the output file

#### Usage

    Log$close()

#### Returns

Invisible TRUE

------------------------------------------------------------------------

### `Log$get_filename()`

Get output filename

#### Usage

    Log$get_filename()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) with the filename

------------------------------------------------------------------------

### `Log$set_filename()`

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

### `Log$delete()`

Delete log file

#### Usage

    Log$delete()

#### Returns

Invisible TRUE

------------------------------------------------------------------------

### `Log$open_system()`

Open log with system viewer

#### Usage

    Log$open_system()

#### Returns

Invisible TRUE

------------------------------------------------------------------------

### `Log$clone()`

The objects of this class are cloneable with this method.

#### Usage

    Log$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
log <- Log$new(tempfile(fileext = ".txt"), new_waiver())
log$get_filename()
#> [1] "file1f5d1db17f5.txt"
```
