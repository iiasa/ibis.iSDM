# Function to remove a latent effect

This is just a wrapper function for removing specified offsets from a
[`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
object.

## Usage

``` r
rm_latent(x)

# S4 method for class 'BiodiversityDistribution'
rm_latent(x)

# S4 method for class 'BiodiversityScenario'
rm_latent(x)
```

## Arguments

- x:

  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

## Value

Removes a latent spatial effect from a
[`distribution`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
object.

## See also

add_latent_spatial

## Examples

``` r
if (FALSE) { # \dontrun{
 rm_latent(model) -> model
} # }
```
