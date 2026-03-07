# Identify local limiting factor

Calculates a
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
of locally limiting factors from a given projected model. To calculate
this first the
[`spartial`](https://iiasa.github.io/ibis.iSDM/reference/spartial.md)
effect of each individual covariate in the model is calculated.

The effect is estimated as that variable most responsible for decreasing
suitability at that cell. The decrease in suitability is calculated, for
each predictor in turn, relative to the suitability that would be
achieved if that predictor took the value equal to the mean. The
predictor associated with the largest decrease in suitability is the
most limiting factor.

## Usage

``` r
limiting(mod, plot = TRUE)

# S4 method for class 'ANY'
limiting(mod, plot = TRUE)
```

## Arguments

- mod:

  A fitted `'DistributionModel'` object from which limited factors are
  to be identified.

- plot:

  Should the result be plotted? (Default: `TRUE`).

## Value

A `terra` object of the most important variable for a given grid cell.

## References

- Elith, J., Kearney, M. and Phillips, S. (2010), The art of modelling
  range-shifting species. Methods in Ecology and Evolution, 1: 330-342.
  doi: 10.1111/j.2041-210X.2010.00036.x

## Examples

``` r
if (FALSE) { # \dontrun{
o <- limiting(fit)
plot(o)
} # }
```
