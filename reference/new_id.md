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
#> id: 83c02868-2f94-4cdd-a9c0-21d62490491f

# convert to character
as.character(i)
#> [1] "83c02868-2f94-4cdd-a9c0-21d62490491f"

# check if it is an Id object
is.Id(i)
#> [1] TRUE
```
