# Engine for Generalized linear models (GLM)

This engine implements a basic generalized linear model (GLM) for
creating species distribution models. The main purpose of this engine is
to support a basic, dependency-free method for inference and projection
that can be used within the package for examples and vignettes. That
being said, the engine is fully functional as any other engine.

The basic implementation of GLMs here is part of a general class of
linear models and has - with exception of offsets - only minimal options
to integrate other sources of information such as priors or joint
integration. The general recommendation is to use
[`engine_glmnet()`](https://iiasa.github.io/ibis.iSDM/reference/engine_glmnet.md)
instead for regularization support. However basic GLMs can in some cases
be useful for quick projections or for
[`ensemble()`](https://iiasa.github.io/ibis.iSDM/reference/ensemble.md)
of small models (a practice common for rare species).

## Usage

``` r
engine_glm(x, control = NULL, type = "response", ...)
```

## Arguments

- x:

  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

- control:

  A [`list`](https://rdrr.io/r/base/list.html) containing parameters for
  controlling the fitting process (Default: `NULL`).

- type:

  The mode used for creating posterior predictions. Either making
  `"link"` or `"response"` (Default: `"response"`).

- ...:

  Other parameters passed on to
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html).

## Value

An
[Engine](https://iiasa.github.io/ibis.iSDM/reference/Engine-class.md).

## Details

This engine is essentially a wrapper for
[`stats::glm.fit()`](https://rdrr.io/r/stats/glm.html), however with
customized settings to support offsets and weights.

If `"optim_hyperparam"` is set to `TRUE` in
[`train()`](https://iiasa.github.io/ibis.iSDM/reference/train.md), then
a AIC based step-wise (backwards) model selection is performed.
Generally however
[`engine_glmnet`](https://iiasa.github.io/ibis.iSDM/reference/engine_glmnet.md)
should be the preferred package for models with more than `>3`
covariates.

## References

- Hastie, T. J. and Pregibon, D. (1992) Generalized linear models.
  Chapter 6 of Statistical Models in S eds J. M. Chambers and T. J.
  Hastie, Wadsworth & Brooks/Cole.

## See also

Other engine:
[`engine_bart()`](https://iiasa.github.io/ibis.iSDM/reference/engine_bart.md),
[`engine_breg()`](https://iiasa.github.io/ibis.iSDM/reference/engine_breg.md),
[`engine_gdb()`](https://iiasa.github.io/ibis.iSDM/reference/engine_gdb.md),
[`engine_glmnet()`](https://iiasa.github.io/ibis.iSDM/reference/engine_glmnet.md),
[`engine_inla()`](https://iiasa.github.io/ibis.iSDM/reference/engine_inla.md),
[`engine_inlabru()`](https://iiasa.github.io/ibis.iSDM/reference/engine_inlabru.md),
[`engine_scampr()`](https://iiasa.github.io/ibis.iSDM/reference/engine_scampr.md),
[`engine_stan()`](https://iiasa.github.io/ibis.iSDM/reference/engine_stan.md),
[`engine_xgboost()`](https://iiasa.github.io/ibis.iSDM/reference/engine_xgboost.md)

## Examples

``` r
# Load background
background <- terra::rast(system.file('extdata/europegrid_50km.tif',
package='ibis.iSDM',mustWork = TRUE))

# Add GLM as an engine
x <- distribution(background) |> engine_glm()
#> [Setup] 2026-05-26 20:27:20.219293 | Creating distribution object...
print(x)
#> <Biodiversity distribution model>
#> Background extent: 
#>      xmin: -16.064, xmax: 34.95,
#>      ymin: 36.322, ymax: 71.535
#>    projection: +proj=longlat +datum=WGS84 +no_defs
#>  --------- 
#> Biodiversity data:
#>    None
#>  --------- 
#>   predictors:     None
#>   priors:         <Default>
#>   latent:         None
#>   log:            <Console>
#>   engine:         <GLM>
```
