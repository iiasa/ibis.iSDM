# Engine for process models using scampr

Similar to others, this engine enables the fitting and prediction of
log-Gaussian Cox process (LGCP) and Inhomogeneous Poisson process (IPP)
processes. It uses the `scampr` package, which uses maximum likelihood
estimation fitted via `TMB` (Template Model Builder).

It also support the addition of spatial latent effects which can be
added via Gaussian fields and approximated by 'FRK' (Fixed Rank Kriging)
and are integrated out using either variational or Laplace
approximation.

The main use case for this engine is as an alternative to
[`engine_inlabru()`](https://iiasa.github.io/ibis.iSDM/reference/engine_inlabru.md)
and
[`engine_inla()`](https://iiasa.github.io/ibis.iSDM/reference/engine_inla.md)
for fitting iSDMs, e.g. those combining both presence-only and
presence-absence point occurrence data.

## Usage

``` r
engine_scampr(x, type = "response", dens = "posterior", maxit = 500, ...)
```

## Arguments

- x:

  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

- type:

  The mode used for creating (posterior or prior) predictions. Either
  setting `"link"` or `"response"` (Default: `"response"`).

- dens:

  A [`character`](https://rdrr.io/r/base/character.html) on how
  predictions are made, either from the `"posterior"` (Default) or
  `"prior"`.

- maxit:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) on the number of
  iterations for the optimizer (Default: `500`).

- ...:

  Other parameters passed on.

## Value

An
[Engine](https://iiasa.github.io/ibis.iSDM/reference/Engine-class.md).

## Details

This engine may only be used to predict for one or two datasets at most.
It supports only presence-only PPMs and presence/absence Binary GLMs, or
'IDM' (for an integrated data model).

## Note

- The package can currently be installed from github directly only
  `"ElliotDovers/scampr"`

- Presence-absence models in SCAMPR currently only support cloglog link
  functions!

## References

- Dovers, E., Popovic, G. C., & Warton, D. I. (2024). A fast method for
  fitting integrated species distribution models. Methods in Ecology and
  Evolution, 15(1), 191-203.

- Dovers, E., Stoklosa, D., and Warton D. I. (2024). Fitting
  log-Gaussian Cox processes using generalized additive model software.
  The American Statistician, 1-17.

## See also

Other engine:
[`engine_bart()`](https://iiasa.github.io/ibis.iSDM/reference/engine_bart.md),
[`engine_breg()`](https://iiasa.github.io/ibis.iSDM/reference/engine_breg.md),
[`engine_gdb()`](https://iiasa.github.io/ibis.iSDM/reference/engine_gdb.md),
[`engine_glm()`](https://iiasa.github.io/ibis.iSDM/reference/engine_glm.md),
[`engine_glmnet()`](https://iiasa.github.io/ibis.iSDM/reference/engine_glmnet.md),
[`engine_inla()`](https://iiasa.github.io/ibis.iSDM/reference/engine_inla.md),
[`engine_inlabru()`](https://iiasa.github.io/ibis.iSDM/reference/engine_inlabru.md),
[`engine_stan()`](https://iiasa.github.io/ibis.iSDM/reference/engine_stan.md),
[`engine_xgboost()`](https://iiasa.github.io/ibis.iSDM/reference/engine_xgboost.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Load background
background <- terra::rast(system.file('extdata/europegrid_50km.tif',
package='ibis.iSDM',mustWork = TRUE))

# Add GLM as an engine
x <- distribution(background) |> engine_scampr()
} # }
```
