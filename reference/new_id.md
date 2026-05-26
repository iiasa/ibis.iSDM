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
#> id: 36102951-7a36-43f0-9b0a-4d21f97c91a6

# convert to character
as.character(i)
#> [1] "36102951-7a36-43f0-9b0a-4d21f97c91a6"

# check if it is an Id object
is.Id(i)
#> [1] TRUE
```
