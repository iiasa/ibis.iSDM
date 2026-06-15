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
#> id: 1d538bcc-1246-4c80-a02c-d1970393c149

# convert to character
as.character(i)
#> [1] "1d538bcc-1246-4c80-a02c-d1970393c149"

# check if it is an Id object
is.Id(i)
#> [1] TRUE
```
