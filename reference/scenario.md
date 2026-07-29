# Create a new scenario based on trained model parameters

This function creates a new
[BiodiversityScenario](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
object that contains the projections of a model.

## Usage

``` r
scenario(fit, limits = NULL, reuse_limits = FALSE, copy_model = FALSE)

# S4 method for class 'ANY'
scenario(fit, limits = NULL, reuse_limits = FALSE, copy_model = FALSE)
```

## Arguments

- fit:

  A
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md)
  object containing a trained model.

- limits:

  A
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  or [`sf`](https://r-spatial.github.io/sf/reference/sf.html) object
  that limits the projection surface when intersected with the
  prediction data (Default: `NULL`). This can for instance be set as an
  expert-delineated constrain to limit spatial projections.

- reuse_limits:

  A [`logical`](https://rdrr.io/r/base/logical.html) on whether to reuse
  limits if found in the trained
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md)
  object (Default: `FALSE`). See also notes!

- copy_model:

  A [`logical`](https://rdrr.io/r/base/logical.html) of whether the
  model object is to be copied to the scenario object. Note that setting
  this option to `TRUE` can increase the required amount of memory
  (Default: `FALSE`).

## Value

A
[`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
object.

## Note

If a limit has been defined already during
[`train()`](https://iiasa.github.io/ibis.iSDM/reference/train.md), for
example by adding an extrapolation limit
[`add_limits_extrapolation()`](https://iiasa.github.io/ibis.iSDM/reference/add_limits_extrapolation.md),
this zonal layer can be reused for the projections. **Note: This
effectively fixes the projections to certain areas.**

## Examples

``` r
if (FALSE) { # \dontrun{
  scenario(fit, limits = island_area)
} # }
```
