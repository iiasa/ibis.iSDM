# Use INLA as engine

Allows a full Bayesian analysis of linear and additive models using
Integrated Nested Laplace Approximation (INLA). This engine has been
largely superseded by the
[`engine_inlabru`](https://iiasa.github.io/ibis.iSDM/reference/engine_inlabru.md)
engine and users are advised to use that one instead, unless specific
options are required.

## Usage

``` r
engine_inla(
  x,
  optional_mesh = NULL,
  optional_projstk = NULL,
  max.edge = NULL,
  offset = NULL,
  cutoff = NULL,
  proj_stepsize = NULL,
  timeout = NULL,
  strategy = "auto",
  int.strategy = "eb",
  barrier = FALSE,
  type = "response",
  area = "gpc2",
  nonconvex.bdry = FALSE,
  nonconvex.convex = -0.15,
  nonconvex.concave = -0.05,
  nonconvex.res = 40,
  ...
)
```

## Arguments

- x:

  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

- optional_mesh:

  A directly supplied `"INLA"` mesh (Default: `NULL`)

- optional_projstk:

  A directly supplied projection stack. Useful if projection stack is
  identical for multiple species (Default: `NULL`)

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
  (Default: `NULL`).

- timeout:

  Specify a timeout for INLA models in sec. Afterwards it passed.

- strategy:

  Which approximation to use for the joint posterior. Options are
  `"auto"` ("default"), `"adaptative"`, `"gaussian"`,
  `"simplified.laplace"` & `"laplace"`.

- int.strategy:

  Integration strategy. Options are `"auto"`,`"grid"`, `"eb"`
  ("default") & `"ccd"`. See also
  https://groups.google.com/g/r-inla-discussion-group/c/hDboQsJ1Mls

- barrier:

  Should a barrier model be added to the model?

- type:

  The mode used for creating posterior predictions. Either summarizing
  the linear `"predictor"` or `"response"` (Default: `"response"`).

- area:

  Accepts a [`character`](https://rdrr.io/r/base/character.html)
  denoting the type of area calculation to be done on the mesh (Default:
  `'gpc2'`).

- nonconvex.bdry:

  Create a non-convex boundary hulls instead (Default: `FALSE`) **Not
  yet implemented**

- nonconvex.convex:

  Non-convex minimal extension radius for convex curvature **Not yet
  implemented**

- nonconvex.concave:

  Non-convex minimal extension radius for concave curvature **Not yet
  implemented**

- nonconvex.res:

  Computation resolution for nonconvex.hulls **Not yet implemented**

- ...:

  Other options.

## Value

An engine.

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

- Havard Rue, Sara Martino, and Nicholas Chopin (2009), Approximate
  Bayesian Inference for Latent Gaussian Models Using Integrated Nested
  Laplace Approximations (with discussion), Journal of the Royal
  Statistical Society B, 71, 319-392.

- Finn Lindgren, Havard Rue, and Johan Lindstrom (2011). An Explicit
  Link Between Gaussian Fields and Gaussian Markov Random Fields: The
  Stochastic Partial Differential Equation Approach (with discussion),
  Journal of the Royal Statistical Society B, 73(4), 423-498.

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
[`engine_inlabru()`](https://iiasa.github.io/ibis.iSDM/reference/engine_inlabru.md),
[`engine_scampr()`](https://iiasa.github.io/ibis.iSDM/reference/engine_scampr.md),
[`engine_stan()`](https://iiasa.github.io/ibis.iSDM/reference/engine_stan.md),
[`engine_xgboost()`](https://iiasa.github.io/ibis.iSDM/reference/engine_xgboost.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Add INLA as an engine (with a custom mesh)
x <- distribution(background) |> engine_inla(mesh = my_mesh)
} # }
```
