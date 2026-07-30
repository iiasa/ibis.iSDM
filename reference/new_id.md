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
#> id: d73c2f5c-1858-4f31-a890-4ea78ed3fb1a

# convert to character
as.character(i)
#> [1] "d73c2f5c-1858-4f31-a890-4ea78ed3fb1a"

# check if it is an Id object
is.Id(i)
#> [1] TRUE
```
