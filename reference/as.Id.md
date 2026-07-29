# As Id

As Id

## Usage

``` r
as.Id(x, ...)

# S3 method for class 'character'
as.Id(x, ...)
```

## Arguments

- x:

  A [`character`](https://rdrr.io/r/base/character.html) to be converted
  as id.

- ...:

  Other arguements

## Value

An object of class `"Id"`.

## Examples

``` r
id <- as.Id("example-id")
is.Id(id)
#> [1] TRUE
```
