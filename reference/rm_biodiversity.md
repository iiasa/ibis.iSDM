# Remove specific BiodiversityDataset from a [distribution](https://iiasa.github.io/ibis.iSDM/reference/distribution.md) object

Remove a particular dataset (or all) from an
[distribution](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
object with a
[`BiodiversityDatasetCollection`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDatasetCollection-class.md).

## Usage

``` r
rm_biodiversity(x, name, id)

# S4 method for class 'BiodiversityDistribution'
rm_biodiversity(x, name, id)
```

## Arguments

- x:

  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

- name:

  A [`character`](https://rdrr.io/r/base/character.html) with the name
  of the biodiversity dataset.

- id:

  A [`character`](https://rdrr.io/r/base/character.html) with the id of
  the biodiversity dataset.

## Value

A
[`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md)
object with matching biodiversity data removed.

## Examples

``` r
if (FALSE) { # \dontrun{
distribution(background) |>
 add_biodiversity_poipa(species, "Duckus communus")
 rm_biodiversity(names = "Duckus communus")
} # }
```
