# Calculate environmental similarity of reference datasets to predictors.

Calculate the environmental similarity of the provided covariates with
respect to a reference dataset. Currently supported is Multivariate
Environmental Similarity index and the multivariate combination novelty
index (NT2) based on the Mahalanobis divergence (see references).

## Usage

``` r
similarity(
  obj,
  ref,
  ref_type = "poipo",
  method = "mess",
  predictor_names = NULL,
  full = FALSE,
  plot = TRUE,
  ...
)

# S4 method for class 'BiodiversityDistribution'
similarity(
  obj,
  ref,
  ref_type = "poipo",
  method = "mess",
  predictor_names = NULL,
  full = FALSE,
  plot = TRUE,
  ...
)

# S4 method for class 'SpatRaster'
similarity(
  obj,
  ref,
  ref_type = "poipo",
  method = "mess",
  predictor_names = NULL,
  full = FALSE,
  plot = TRUE,
  ...
)
```

## Arguments

- obj:

  A
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md),
  [`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
  or alternatively a
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object.

- ref:

  A
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md),
  [`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
  or alternatively a
  [`data.frame`](https://rdrr.io/r/base/data.frame.html) with extracted
  values (corresponding to those given in `obj`).

- ref_type:

  A [`character`](https://rdrr.io/r/base/character.html) specifying the
  type of biodiversity to use when obj is a
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md).

- method:

  A specifc method for similarity calculation. Currently supported:
  `'mess'`, `'nt'`.

- predictor_names:

  An optional [`character`](https://rdrr.io/r/base/character.html)
  specifying the covariates to be used (Default: `NULL`).

- full:

  should similarity values be returned for all variables
  (Default:`FALSE`)?

- plot:

  Should the result be plotted? Otherwise return the output list
  (Default: `TRUE`).

- ...:

  other options (Non specified).

## Value

This function returns a list containing:

- `similarity`: A `SpatRaster` object with multiple layers giving the
  environmental similarities for each variable in `x` (only included
  when `"full=TRUE"`);

- `mis`: a `SpatRaster` layer giving the minimum similarity value across
  all variables for each location (i.e. the MESS);

- `exip`: a `SpatRaster` layer indicating whether any model would
  interpolate or extrapolate to this location based on environmental
  surface;

- `mod`: a factor `SpatRaster` layer indicating which variable was most
  dissimilar to its reference range (i.e. the MoD map, Elith et al.
  2010); and

- `mos`: a factor `SpatRaster` layer indicating which variable was most
  similar to its reference range.

## Details

`similarity` implements the MESS algorithm described in Appendix S3 of
Elith et al. (2010) as well as the Mahalanobis dissimilarity described
in Mesgaran et al. (2014).

## References

- Elith, J., Kearney, M., and Phillips, S. (2010) "The art of modelling
  range-shifting species". *Methods in Ecology and Evolution*, 1:
  330-342. https://doi.org/10.1111/j.2041-210X.2010.00036.x

- Mesgaran, M.B., Cousens, R.D. and Webber, B.L. (2014) "Here be
  dragons: a tool for quantifying novelty due to covariate range and
  correlation change when projecting species distribution models".
  *Diversity and Distributions*, 20: 1147-1159.
  https://doi.org/10.1111/ddi.12209

## See also

dismo R-package.

## Examples

``` r
 if (FALSE) { # \dontrun{
plot(
  similarity(x) # Where x is a distribution or Raster object
)
 } # }
```
