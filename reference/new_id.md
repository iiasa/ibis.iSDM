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
#> id: fa3fc517-7846-4d72-a4e5-afdf37f8236b

# convert to character
as.character(i)
#> [1] "fa3fc517-7846-4d72-a4e5-afdf37f8236b"

# check if it is an Id object
is.Id(i)
#> [1] TRUE
```
