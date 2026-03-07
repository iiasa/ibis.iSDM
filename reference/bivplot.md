# Bivariate prediction plot for distribution objects

Often there is an intention to display not only the predictions made
with a SDM, but also the uncertainty of the prediction. Uncertainty be
estimated either directly by the model or by calculating the variation
in prediction values among a set of models.

In particular Bayesian engines can produce not only mean estimates of
fitted responses, but also pixel-based estimates of uncertainty from the
posterior such as the standard deviation (SD) or the coefficient of
variation of a given prediction.

This function makes use of the `"biscale"` R-package to create bivariate
plots of the fitted distribution object, allowing to visualize two
variables at once. It is mostly thought of as a convenience function to
create such bivariate plots for quick visualization.

Supported Inputs are either single trained Bayesian
[`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
with uncertainty or the output of an
[`ensemble()`](https://iiasa.github.io/ibis.iSDM/reference/ensemble.md)
call. In both cases, users have to make sure that `"xvar"` and `"yvar"`
are set accordingly.

## Usage

``` r
bivplot(
  mod,
  xvar = "mean",
  yvar = "sd",
  plot = TRUE,
  fname = NULL,
  title = NULL,
  col = "BlueGold",
  ...
)

# S4 method for class 'ANY'
bivplot(
  mod,
  xvar = "mean",
  yvar = "sd",
  plot = TRUE,
  fname = NULL,
  title = NULL,
  col = "BlueGold",
  ...
)
```

## Arguments

- mod:

  A trained
  [`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
  or alternatively a
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object with `prediction` model within.

- xvar:

  A [`character`](https://rdrr.io/r/base/character.html) denoting the
  value on the x-axis (Default: `'mean'`).

- yvar:

  A [`character`](https://rdrr.io/r/base/character.html) denoting the
  value on the y-axis (Default: `'sd'`).

- plot:

  A [`logical`](https://rdrr.io/r/base/logical.html) indication of
  whether the result is to be plotted (Default: `TRUE`)?

- fname:

  A [`character`](https://rdrr.io/r/base/character.html) specifying the
  output filename a created figure should be written to.

- title:

  Allows to respecify the title through a
  [`character`](https://rdrr.io/r/base/character.html) (Default:`NULL`).

- col:

  A [`character`](https://rdrr.io/r/base/character.html) stating the
  colour palette to use. Has to be either a predefined value or a vector
  of colours. See `"biscale::bi_pal_manual"`. Default: `"BlueGold"`.

- ...:

  Other engine specific parameters.

## Value

Saved bivariate plot in `'fname'` if specified, otherwise plot.

## Note

**This function requires the biscale package to be installed.** Although
a work around without the package could be developed, it was not deemed
necessary at this point. See also this
[gist](https://gist.github.com/scbrown86/2779137a9378df7b60afd23e0c45c188).

## See also

[partial](https://iiasa.github.io/ibis.iSDM/reference/partial.md),
[plot.DistributionModel](https://iiasa.github.io/ibis.iSDM/reference/plot.md)
