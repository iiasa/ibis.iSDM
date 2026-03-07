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
#> id: 9bf64093-f62c-43ac-838a-8c26f1022a7a

# convert to character
as.character(i)
#> [1] "9bf64093-f62c-43ac-838a-8c26f1022a7a"

# check if it is an Id object
is.Id(i)
#> [1] TRUE
```
