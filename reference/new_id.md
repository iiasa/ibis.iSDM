# Identifier

Generate a new unique identifier.

## Usage

``` r
new_id()
```

## Value

`"Id"` object.

## Details

Identifiers are made using the
[`uuid::UUIDgenerate()`](https://rdrr.io/pkg/uuid/man/UUIDgenerate.html).

## See also

[`uuid::UUIDgenerate()`](https://rdrr.io/pkg/uuid/man/UUIDgenerate.html).

## Examples

``` r
# create new id
i <- new_id()

# print id
print(i)
#> id: 516350c6-b595-4eeb-8d74-1c7e4fcc1eea

# convert to character
as.character(i)
#> [1] "516350c6-b595-4eeb-8d74-1c7e4fcc1eea"

# check if it is an Id object
is.Id(i)
#> [1] TRUE
```
