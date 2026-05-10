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
#> id: 9d54b2f1-e024-4ef4-ba49-17fed21e5243

# convert to character
as.character(i)
#> [1] "9d54b2f1-e024-4ef4-ba49-17fed21e5243"

# check if it is an Id object
is.Id(i)
#> [1] TRUE
```
