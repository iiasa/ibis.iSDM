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
#> id: c24f1d59-ec51-4e07-b1cf-607547ca8612

# convert to character
as.character(i)
#> [1] "c24f1d59-ec51-4e07-b1cf-607547ca8612"

# check if it is an Id object
is.Id(i)
#> [1] TRUE
```
