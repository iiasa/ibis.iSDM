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
#> id: c763d8a8-a1b8-4611-8265-f468c31af038

# convert to character
as.character(i)
#> [1] "c763d8a8-a1b8-4611-8265-f468c31af038"

# check if it is an Id object
is.Id(i)
#> [1] TRUE
```
