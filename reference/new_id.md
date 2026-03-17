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
#> id: 3f37cf67-8cbb-4c49-9691-88d981bbd11e

# convert to character
as.character(i)
#> [1] "3f37cf67-8cbb-4c49-9691-88d981bbd11e"

# check if it is an Id object
is.Id(i)
#> [1] TRUE
```
