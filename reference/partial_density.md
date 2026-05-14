# Visualize the density of the data over the environmental data

Based on a fitted model, plot the density of observations over the
estimated variable and environmental space. Opposed to the
[partial](https://iiasa.github.io/ibis.iSDM/reference/partial.md) and
[spartial](https://iiasa.github.io/ibis.iSDM/reference/spartial.md)
functions, which are rather low-level interfaces, this function provides
more detail in the light of the data. It is also able to contrast
different variables against each other and show the used data.

## Usage

``` r
partial_density(mod, x.var, df = FALSE, ...)

# S4 method for class 'ANY,character'
partial_density(mod, x.var, df = FALSE, ...)
```

## Arguments

- mod:

  A trained `DistributionModel` object. Requires a fitted model and
  inferred prediction.

- x.var:

  A [character](https://rdrr.io/r/base/character.html) indicating the
  variable to be investigated. Can be a
  [`vector`](https://rdrr.io/r/base/vector.html) of length `1` or `2`.

- df:

  [`logical`](https://rdrr.io/r/base/logical.html) if plotting data
  should be returned instead (Default: `FALSE`).

- ...:

  Other engine specific parameters.

## Value

A
[`ggplot2`](https://ggplot2.tidyverse.org/reference/ggplot2-package.html)
object showing the marginal response in light of the data.

## Details

This functions calculates the observed density of presence and absence
points over the whole surface of a specific variable. It can be used to
visually inspect the fit of the model to data.

## Note

By default all variables that are not `x.var` are hold constant at the
mean.

## References

- Warren, D.L., Matzke, N.J., Cardillo, M., Baumgartner, J.B., Beaumont,
  L.J., Turelli, M., Glor, R.E., Huron, N.A., Simões, M., Iglesias, T.L.
  Piquet, J.C., and Dinnage, R. 2021. ENMTools 1.0: an R package for
  comparative ecological biogeography. Ecography, 44(4), pp.504-511.

## See also

[partial](https://iiasa.github.io/ibis.iSDM/reference/partial.md)

## Examples

``` r
if (FALSE) { # \dontrun{
 # Do a partial calculation of a trained model
 partial_density(fit, x.var = "Forest.cover")
 # Or with two variables
 partial_density(fit, x.var = c("Forest.cover", "bio01"))
} # }
```
