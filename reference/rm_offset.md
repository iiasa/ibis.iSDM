# Function to remove an offset

This is just a wrapper function for removing specified offsets from a
[`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
object.

## Usage

``` r
rm_offset(x, layer = NULL)

# S4 method for class 'BiodiversityDistribution'
rm_offset(x, layer = NULL)
```

## Arguments

- x:

  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

- layer:

  A `character` pointing to the specific layer to be removed. If set to
  `NULL`, then all offsets are removed from the object.

## Value

Removes an offset from a
[`distribution`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
object.

## See also

Other offset:
[`add_offset()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset.md),
[`add_offset_bias()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset_bias.md),
[`add_offset_elevation()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset_elevation.md),
[`add_offset_range()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset_range.md)

## Examples

``` r
if (FALSE) { # \dontrun{
 rm_offset(model) -> model
} # }
```
