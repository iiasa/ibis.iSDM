# Base Prior class

This class sets up the base class for priors which will be inherited by
all priors.

## Value

Defines a Prior object.

## Note

This functionality likely is deprecated or checks have been superseeded.

## Public fields

- `id`:

  A [`character`](https://rdrr.io/r/base/character.html) with the id of
  the prior.

- `name`:

  A [`character`](https://rdrr.io/r/base/character.html) with the name
  of the prior.

- `type`:

  A [`character`](https://rdrr.io/r/base/character.html) with the type
  of the prior.

- `variable`:

  A [`character`](https://rdrr.io/r/base/character.html) with the
  variable name for the prior.

- `distribution`:

  A [`character`](https://rdrr.io/r/base/character.html) with the
  distribution of the prior if relevant.

- `value`:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) or
  [`character`](https://rdrr.io/r/base/character.html) with the prior
  value, e.g. the hyper-parameters.

- `prob`:

  Another [`numeric`](https://rdrr.io/r/base/numeric.html) entry on the
  prior field. The inclusion probability.

- `lims`:

  A limitation on the lower and upper bounds of a numeric value.

## Methods

### Public methods

- [`Prior$new()`](#method-Prior-initialize)

- [`Prior$print()`](#method-Prior-print)

- [`Prior$validate()`](#method-Prior-validate)

- [`Prior$get()`](#method-Prior-get)

- [`Prior$set()`](#method-Prior-set)

- [`Prior$get_id()`](#method-Prior-get_id)

- [`Prior$get_name()`](#method-Prior-get_name)

- [`Prior$clone()`](#method-Prior-clone)

------------------------------------------------------------------------

### `Prior$new()`

Initializes the object and prepared the various prior variables

#### Usage

    Prior$new(
      id,
      name,
      variable,
      value,
      type = NULL,
      distribution = NULL,
      prob = NULL,
      lims = NULL
    )

#### Arguments

- `id`:

  A [`character`](https://rdrr.io/r/base/character.html) with the id of
  the prior.

- `name`:

  A [`character`](https://rdrr.io/r/base/character.html) with the name
  of the prior.

- `variable`:

  A [`character`](https://rdrr.io/r/base/character.html) with the
  variable name for the prior.

- `value`:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) or
  [`character`](https://rdrr.io/r/base/character.html) with the prior
  value, e.g. the hyper-parameters.

- `type`:

  A [`character`](https://rdrr.io/r/base/character.html) with the type
  of the prior.

- `distribution`:

  A [`character`](https://rdrr.io/r/base/character.html) with the
  distribution of the prior if relevant.

- `prob`:

  Another [`numeric`](https://rdrr.io/r/base/numeric.html) entry on the
  prior field. The inclusion probability.

- `lims`:

  A limitation on the lower and upper bounds of a numeric value.

#### Returns

NULL

------------------------------------------------------------------------

### `Prior$print()`

Print out the prior type and variable.

#### Usage

    Prior$print()

#### Returns

A message on screen

------------------------------------------------------------------------

### `Prior$validate()`

Generic validation function for a provided value.

#### Usage

    Prior$validate(x)

#### Arguments

- `x`:

  A new prior value.

#### Returns

Invisible TRUE

------------------------------------------------------------------------

### `Prior$get()`

Get prior values

#### Usage

    Prior$get(what = "value")

#### Arguments

- `what`:

  A [`character`](https://rdrr.io/r/base/character.html) with the entry
  to be returned (Default: `value`).

#### Returns

Invisible TRUE

------------------------------------------------------------------------

### `Prior$set()`

Set prior

#### Usage

    Prior$set(x)

#### Arguments

- `x`:

  A new prior value as [`numeric`](https://rdrr.io/r/base/numeric.html)
  or [`character`](https://rdrr.io/r/base/character.html).

#### Returns

Invisible TRUE

------------------------------------------------------------------------

### `Prior$get_id()`

Get a specific ID from a prior.

#### Usage

    Prior$get_id()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) id.

------------------------------------------------------------------------

### `Prior$get_name()`

Get Name of object

#### Usage

    Prior$get_name()

#### Returns

Returns a [`character`](https://rdrr.io/r/base/character.html) with the
class name.

------------------------------------------------------------------------

### `Prior$clone()`

The objects of this class are cloneable with this method.

#### Usage

    Prior$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
prior <- Prior$new(name = "ExamplePrior", variable = "temperature", value = 1)
prior$get("value")
#> [1] 1
```
