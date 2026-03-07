# Adds a log file to distribution object

This function allows to specify a file as
[Log](https://iiasa.github.io/ibis.iSDM/reference/Log-class.md) file,
which is used to save all console outputs, prints and messages.

## Usage

``` r
add_log(x, filename)

# S4 method for class 'BiodiversityDistribution,character'
add_log(x, filename)
```

## Arguments

- x:

  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

- filename:

  A [`character`](https://rdrr.io/r/base/character.html) object. The
  destination must be writeable and filename ends with `'txt'`.

## Value

Adds a log file to a
[`distribution`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
object.

## Examples

``` r
if (FALSE) { # \dontrun{
 x <- distribution(background) |>
    add_log()
 x
} # }
```
