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
#> id: f9ef7eca-450f-4c46-9712-0cc304ab51e0

# convert to character
as.character(i)
#> [1] "f9ef7eca-450f-4c46-9712-0cc304ab51e0"

# check if it is an Id object
is.Id(i)
#> [1] TRUE
```
