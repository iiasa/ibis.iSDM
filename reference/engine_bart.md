# Engine for use of Bayesian Additive Regression Trees (BART)

The Bayesian regression approach to a sum of complementary trees is to
shrink the said fit of each tree through a regularization prior. BART
models provide non-linear highly flexible estimation and have been shown
to compare favourable among machine learning algorithms (Dorie et al.
2019). Default prior preference is for trees to be small (few terminal
nodes) and shrinkage towards `0`.

This package requires the `"dbarts"` R-package to be installed. Many of
the functionalities of this engine have been inspired by the
`"embarcadero"` R-package. Users are therefore advised to cite if they
make heavy use of BART.

## Usage

``` r
engine_bart(x, iter = 1000, nburn = 250, chains = 4, type = "response", ...)
```

## Arguments

- x:

  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

- iter:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) estimate of the
  number of trees to be used in the sum-of-trees formulation (Default:
  `1000`).

- nburn:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) estimate of the
  burn in samples (Default: `250`).

- chains:

  A number of the number of chains to be used (Default: `4`).

- type:

  The type used for creating posterior predictions. Either `"link"` or
  `"response"` (Default: `"response"`).

- ...:

  Other options.

## Value

An
[Engine](https://iiasa.github.io/ibis.iSDM/reference/Engine-class.md).

## Details

Prior distributions can furthermore be set for:

- probability that a tree stops at a node of a given depth (Not yet
  implemented)

- probability that a given variable is chosen for a splitting rule

- probability of splitting that variable at a particular value (Not yet
  implemented)

## References

- Carlson, CJ. embarcadero: Species distribution modelling with Bayesian
  additive regression trees in r. Methods Ecol Evol. 2020; 11: 850– 858.
  https://doi.org/10.1111/2041-210X.13389

- Dorie, V., Hill, J., Shalit, U., Scott, M., & Cervone, D. (2019).
  Automated versus do-it-yourself methods for causal inference: Lessons
  learned from a data analysis competition. Statistical Science, 34(1),
  43-68.

- Vincent Dorie (2020). dbarts: Discrete Bayesian Additive Regression
  Trees Sampler. R package version 0.9-19.
  https://CRAN.R-project.org/package=dbarts

## See also

Other engine:
[`engine_breg()`](https://iiasa.github.io/ibis.iSDM/reference/engine_breg.md),
[`engine_gdb()`](https://iiasa.github.io/ibis.iSDM/reference/engine_gdb.md),
[`engine_glm()`](https://iiasa.github.io/ibis.iSDM/reference/engine_glm.md),
[`engine_glmnet()`](https://iiasa.github.io/ibis.iSDM/reference/engine_glmnet.md),
[`engine_inla()`](https://iiasa.github.io/ibis.iSDM/reference/engine_inla.md),
[`engine_inlabru()`](https://iiasa.github.io/ibis.iSDM/reference/engine_inlabru.md),
[`engine_scampr()`](https://iiasa.github.io/ibis.iSDM/reference/engine_scampr.md),
[`engine_stan()`](https://iiasa.github.io/ibis.iSDM/reference/engine_stan.md),
[`engine_xgboost()`](https://iiasa.github.io/ibis.iSDM/reference/engine_xgboost.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Add BART as an engine
x <- distribution(background) |> engine_bart(iter = 100)
} # }
```
