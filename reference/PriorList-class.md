# List of Priors supplied to an class

This class represents a collection of
[`Prior`](https://iiasa.github.io/ibis.iSDM/reference/Prior-class.md)
objects. It provides methods for accessing, adding and removing priors
from the list

## Value

A PriorList object.

## Public fields

- `priors`:

  A list of
  [`Prior`](https://iiasa.github.io/ibis.iSDM/reference/Prior-class.md)
  object.

## Methods

### Public methods

- [`PriorList$new()`](#method-PriorList-initialize)

- [`PriorList$print()`](#method-PriorList-print)

- [`PriorList$show()`](#method-PriorList-show)

- [`PriorList$length()`](#method-PriorList-length)

- [`PriorList$ids()`](#method-PriorList-ids)

- [`PriorList$varnames()`](#method-PriorList-varnames)

- [`PriorList$classes()`](#method-PriorList-classes)

- [`PriorList$types()`](#method-PriorList-types)

- [`PriorList$exists()`](#method-PriorList-exists)

- [`PriorList$add()`](#method-PriorList-add)

- [`PriorList$get()`](#method-PriorList-get)

- [`PriorList$collect()`](#method-PriorList-collect)

- [`PriorList$rm()`](#method-PriorList-rm)

- [`PriorList$summary()`](#method-PriorList-summary)

- [`PriorList$combine()`](#method-PriorList-combine)

- [`PriorList$clone()`](#method-PriorList-clone)

------------------------------------------------------------------------

### `PriorList$new()`

Initializes the object

#### Usage

    PriorList$new(priors)

#### Arguments

- `priors`:

  A list of
  [`Prior`](https://iiasa.github.io/ibis.iSDM/reference/Prior-class.md)
  object.

#### Returns

NULL

------------------------------------------------------------------------

### `PriorList$print()`

Print out summary statistics

#### Usage

    PriorList$print()

#### Returns

A message on screen

------------------------------------------------------------------------

### `PriorList$show()`

Aliases that calls print.

#### Usage

    PriorList$show()

#### Returns

A message on screen

------------------------------------------------------------------------

### `PriorList$length()`

Number of priors in object

#### Usage

    PriorList$length()

#### Returns

A [`numeric`](https://rdrr.io/r/base/numeric.html) with the number of
priors set

------------------------------------------------------------------------

### `PriorList$ids()`

Ids of prior objects

#### Usage

    PriorList$ids()

#### Returns

A list with ids of the priors objects for query

------------------------------------------------------------------------

### `PriorList$varnames()`

Variable names of priors in object

#### Usage

    PriorList$varnames()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) list with the
variable names of the priors.

------------------------------------------------------------------------

### `PriorList$classes()`

Function to return the classes of all contained priors

#### Usage

    PriorList$classes()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) list with the
class names of the priors.

------------------------------------------------------------------------

### `PriorList$types()`

Get types of all contained priors

#### Usage

    PriorList$types()

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) list with the
type names of the priors.

------------------------------------------------------------------------

### `PriorList$exists()`

Does a certain variable or type combination exist as prior ?

#### Usage

    PriorList$exists(variable, type = NULL)

#### Arguments

- `variable`:

  A [`character`](https://rdrr.io/r/base/character.html) with the
  variable name.

- `type`:

  A [`character`](https://rdrr.io/r/base/character.html) with the type.

#### Returns

A [`character`](https://rdrr.io/r/base/character.html) id.

------------------------------------------------------------------------

### `PriorList$add()`

Add a new prior to the object.

#### Usage

    PriorList$add(p)

#### Arguments

- `p`:

  A
  [`Prior`](https://iiasa.github.io/ibis.iSDM/reference/Prior-class.md)
  object.

#### Returns

Invisible TRUE

------------------------------------------------------------------------

### `PriorList$get()`

Get specific prior values from the list if set

#### Usage

    PriorList$get(variable, type = NULL, what = "value")

#### Arguments

- `variable`:

  A [`character`](https://rdrr.io/r/base/character.html) with the
  variable name.

- `type`:

  A [`character`](https://rdrr.io/r/base/character.html) with the type
  name

- `what`:

  A [`character`](https://rdrr.io/r/base/character.html) on the specific
  entry to return (Default: `prior value`).

#### Returns

The prior object.

------------------------------------------------------------------------

### `PriorList$collect()`

Collect priors for a given id or multiple.

#### Usage

    PriorList$collect(id)

#### Arguments

- `id`:

  A [`character`](https://rdrr.io/r/base/character.html) with the prior
  id.

#### Returns

A `PriorList` object.

------------------------------------------------------------------------

### `PriorList$rm()`

Remove a set prior by id

#### Usage

    PriorList$rm(id)

#### Arguments

- `id`:

  A [`character`](https://rdrr.io/r/base/character.html) with the prior
  id.

#### Returns

Invisible TRUE

------------------------------------------------------------------------

### `PriorList$summary()`

Summary function that lists all priors

#### Usage

    PriorList$summary()

#### Returns

A [`data.frame`](https://rdrr.io/r/base/data.frame.html) with the
summarized priors.

------------------------------------------------------------------------

### `PriorList$combine()`

Combining function to combine this PriorList with another new one

#### Usage

    PriorList$combine(x)

#### Arguments

- `x`:

  A new `PriorList` object.

#### Returns

Invisible TRUE

------------------------------------------------------------------------

### `PriorList$clone()`

The objects of this class are cloneable with this method.

#### Usage

    PriorList$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
if (FALSE) { # \dontrun{
priors(
    INLAPrior('var1','normal',c(0,0.1)),
    INLAPrior('var2','normal',c(0,0.1))
   )
} # }
```
