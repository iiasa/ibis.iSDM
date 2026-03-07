# Use inlabru as engine

Model components are specified with general inputs and mapping methods
to the latent variables, and the predictors are specified via general R
expressions, with separate expressions for each observation likelihood
model in multi-likelihood models. The inlabru engine - similar as the
[`engine_inla`](https://iiasa.github.io/ibis.iSDM/reference/engine_inla.md)
function acts a wrapper for INLA, albeit `"inlabru"` has a number of
convenience functions implemented that make in particular predictions
with new data much more straight forward (e.g. via posterior simulation
instead of fitting). Since more recent versions `"inlabru"` also
supports the addition of multiple likelihoods, therefore allowing full
integrated inference.

## Usage

``` r
engine_inlabru(
  x,
  optional_mesh = NULL,
  max.edge = NULL,
  offset = NULL,
  cutoff = NULL,
  proj_stepsize = NULL,
  strategy = "auto",
  int.strategy = "eb",
  area = "gpc2",
  timeout = NULL,
  type = "response",
  ...
)
```

## Source

<https://inlabru-org.github.io/inlabru/articles/>

## Arguments

- x:

  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

- optional_mesh:

  A directly supplied `"INLA"` mesh (Default: `NULL`)

- max.edge:

  The largest allowed triangle edge length, must be in the same scale
  units as the coordinates. Default is an educated guess (Default:
  `NULL`).

- offset:

  interpreted as a numeric factor relative to the approximate data
  diameter. Default is an educated guess (Default: `NULL`).

- cutoff:

  The minimum allowed distance between points on the mesh. Default is an
  educated guess (Default: `NULL`).

- proj_stepsize:

  The stepsize in coordinate units between cells of the projection grid
  (Default: `NULL`)

- strategy:

  Which approximation to use for the joint posterior. Options are
  `"auto"` ("default"), `"adaptative"`, `"gaussian"`,
  `"simplified.laplace"` & `"laplace"`.

- int.strategy:

  Integration strategy. Options are `"auto"`, `"grid"`, `"eb"`
  ("default") & `"ccd"`.

- area:

  Accepts a [`character`](https://rdrr.io/r/base/character.html)
  denoting the type of area calculation to be done on the mesh (Default:
  `'gpc2'`).

- timeout:

  Specify a timeout for INLA models in sec. Afterwards it passed.

- type:

  The mode used for creating posterior predictions. Either summarizing
  the linear `"predictor"` or `"response"` (Default:`"response"`).

- ...:

  Other variables

## Value

An
[Engine](https://iiasa.github.io/ibis.iSDM/reference/Engine-class.md).

## Details

All `INLA` engines require the specification of a mesh that needs to be
provided to the `"optional_mesh"` parameter. Otherwise the mesh will be
created based on best guesses of the data spread. A good mesh needs to
have triangles as regular as possible in size and shape: equilateral.

- `"max.edge"`: The largest allowed triangle edge length, must be in the
  same scale units as the coordinates. Lower bounds affect the density
  of triangles.

- `"offset"`: The automatic extension distance of the mesh. If positive:
  same scale units. If negative, interpreted as a factor relative to the
  approximate data diameter, i.e., a value of -0.10 will add a 10% of
  the data diameter as outer extension.

- `"cutoff"`: The minimum allowed distance between points, it means that
  points at a closer distance than the supplied value are replaced by a
  single vertex. It is critical when there are some points very close to
  each other, either for point locations or in the domain boundary.

- `"proj_stepsize"`: The stepsize for spatial predictions, which affects
  the spatial grain of any outputs created.

Priors can be set via
[INLAPrior](https://iiasa.github.io/ibis.iSDM/reference/INLAPrior.md).

## Note

**How INLA Meshes are generated, substantially influences prediction
outcomes. See Dambly et al. (2023).**

## References

- Bachl, F. E., Lindgren, F., Borchers, D. L., & Illian, J. B. (2019).
  inlabru: an R package for Bayesian spatial modelling from ecological
  survey data. Methods in Ecology and Evolution, 10(6), 760-766.

- Simpson, Daniel, Janine B. Illian, S. H. Sørbye, and Håvard Rue. 2016.
  “Going Off Grid: Computationally Efficient Inference for Log-Gaussian
  Cox Processes.” Biometrika 1 (103): 49–70.

- Dambly, L. I., Isaac, N. J., Jones, K. E., Boughey, K. L., &
  O'Hara, R. B. (2023). Integrated species distribution models fitted in
  INLA are sensitive to mesh parameterisation. Ecography, e06391.

## See also

Other engine:
[`engine_bart()`](https://iiasa.github.io/ibis.iSDM/reference/engine_bart.md),
[`engine_breg()`](https://iiasa.github.io/ibis.iSDM/reference/engine_breg.md),
[`engine_gdb()`](https://iiasa.github.io/ibis.iSDM/reference/engine_gdb.md),
[`engine_glm()`](https://iiasa.github.io/ibis.iSDM/reference/engine_glm.md),
[`engine_glmnet()`](https://iiasa.github.io/ibis.iSDM/reference/engine_glmnet.md),
[`engine_inla()`](https://iiasa.github.io/ibis.iSDM/reference/engine_inla.md),
[`engine_scampr()`](https://iiasa.github.io/ibis.iSDM/reference/engine_scampr.md),
[`engine_stan()`](https://iiasa.github.io/ibis.iSDM/reference/engine_stan.md),
[`engine_xgboost()`](https://iiasa.github.io/ibis.iSDM/reference/engine_xgboost.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Add inlabru as an engine
x <- distribution(background) |> engine_inlabru()
} # }
```
