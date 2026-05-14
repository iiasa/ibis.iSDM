# Add latent spatial effect to the model equation

In general we understand under latent spatial effects the occurrence of
spatial dependency in the observations, which might either be caused by
spatial biases, similarities in the underlying sampling processes or
unmeasured latent covariates, e.g. those that have not been quantified.

This package supports a range of different spatial effects, however they
differ from another by their impact on the estimated prediction. Some
effects simply add the spatial dependence as covariate, others make use
of spatial random effects to account for spatial dependence in the
predictions. By default these effects are added to each dataset as
covariate or shared spatial field (e.g. SPDE). See details for an
explanation of the available options.

## Usage

``` r
add_latent_spatial(
  x,
  method = "spde",
  priors = NULL,
  separate_spde = FALSE,
  ...
)

# S4 method for class 'BiodiversityDistribution'
add_latent_spatial(
  x,
  method = "spde",
  priors = NULL,
  separate_spde = FALSE,
  ...
)

# S4 method for class 'BiodiversityScenario'
add_latent_spatial(x, layer = NULL, reuse_latent = TRUE, ...)
```

## Arguments

- x:

  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

- method:

  A [`character`](https://rdrr.io/r/base/character.html) describing what
  kind of spatial effect is to be added to the model. See details.

- priors:

  A `"Prior-List"` object supplied to the latent effect. Relevant only
  for
  [`engine_inla`](https://iiasa.github.io/ibis.iSDM/reference/engine_inla.md)
  and `NULL` equates the use of default priors.

- separate_spde:

  A [`logical`](https://rdrr.io/r/base/logical.html) parameter
  indicating whether, in the case of SPDE effects, separate effects for
  each likelihood are being fitted. Default (`FALSE`) uses a copy of the
  first added likelihood.

- ...:

  Other parameters passed down

- layer:

  A
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  layer describing alternative latent effects to be used instead if
  `"reuse_latent"` is set to `FALSE`.

- reuse_latent:

  A [`logical`](https://rdrr.io/r/base/logical.html) flag on whether any
  latent effects found in the fitted model should be reused (Default
  `TRUE`).

## Value

Adds latent spatial effect to a
[`distribution`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
object.

## Details

There are several different options some of which depend on the engine
used. In case a unsupported method for an engine is chosen this is
modified to the next similar method.

**Available are:**

- `"spde"` - stochastic partial differential equation (SPDE) for
  [`engine_inla`](https://iiasa.github.io/ibis.iSDM/reference/engine_inla.md)
  and
  [`engine_inlabru`](https://iiasa.github.io/ibis.iSDM/reference/engine_inlabru.md).
  SPDE effects aim at capturing the variation of the response variable
  in space, once all of the covariates are accounted for. Examining the
  spatial distribution of the spatial error can reveal which covariates
  might be missing. For example, if elevation is positively correlated
  with the response variable, but is not included in the model, we could
  see a higher posterior mean in areas with higher elevation. Note that
  calculations of SPDE's can be computationally costly.

- `"car"` - conditional autocorrelative errors (CAR) for
  [`engine_inla`](https://iiasa.github.io/ibis.iSDM/reference/engine_inla.md).
  Not yet implemented in full.

- `"kde"` - additional covariate of the kernel density of input point
  observations.

- `"poly"` - spatial trend correction by adding coordinates as
  polynominal transformation. This method assumed that a transformation
  of spatial coordinates can if - included as additional predictor -
  explain some of the variance in the distribution. This method does not
  interact with species occurrences.

- `"nnd"` - nearest neighbour distance. This function calculates the
  euclidean distance from each point to the nearest other grid cell with
  known species occurrence. Originally proposed by Allouche et
  al. (2008) and can be applied across all datasets in the
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

## References

- Allouche, O.; Steinitz, O.; Rotem, D.; Rosenfeld, A.; Kadmon, R.
  (2008). Incorporating distance constraints into species distribution
  models. Journal of Applied Ecology, 45(2), 599-609.
  doi:10.1111/j.1365-2664.2007.01445.x

- Mendes, P., Velazco, S. J. E., de Andrade, A. F. A., &
  Júnior, P. D. M. (2020). Dealing with overprediction in species
  distribution models: How adding distance constraints can improve model
  accuracy. Ecological Modelling, 431, 109180.

## Examples

``` r
if (FALSE) { # \dontrun{
 distribution(background) |> add_latent_spatial(method = "poly")
} # }
```
