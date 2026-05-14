# Threshold a continuous prediction to a categorical layer

It is common in many applications of species distribution modelling that
estimated continuous suitability surfaces are converted into discrete
representations of where suitable habitat might or might not exist. This
so called *threshold'ing* can be done in various ways which are further
described in the details.

In case a
[`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
is provided as input in this function for `obj`, it is furthermore
necessary to provide a
[`sf`](https://r-spatial.github.io/sf/reference/sf.html) object for
validation as there is no
[`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
to read this information from.

**Note:** This of course also allows to estimate the threshold based on
withheld data, for instance those created from an a-priori
cross-validation procedure.

For
[`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
objects, adding this function to the processing pipeline stores a
threshold attribute in the created
[scenario](https://iiasa.github.io/ibis.iSDM/reference/scenario.md)
object.

For
[BiodiversityScenario](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
objects a set `threshold()` simply indicates that the projection should
create and use thresholds as part of the results. The threshold values
for this are either taken from the provided model or through an optional
provide parameter `value`.

If instead the aim is to apply thresholds to each step of the
suitability projection, see
[`add_constraint_threshold()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_threshold.md).

## Usage

``` r
threshold(
  obj,
  method = "mtp",
  value = NULL,
  point = NULL,
  field_occurrence = "observed",
  format = "binary",
  return_threshold = FALSE,
  ...
)

# S4 method for class 'ANY'
threshold(
  obj,
  method = "mtp",
  value = NULL,
  point = NULL,
  field_occurrence = "observed",
  format = "binary",
  return_threshold = FALSE,
  ...
)

# S4 method for class 'SpatRaster'
threshold(
  obj,
  method = "fixed",
  value = NULL,
  point = NULL,
  field_occurrence = "observed",
  format = "binary",
  return_threshold = FALSE
)

# S4 method for class 'BiodiversityScenario'
threshold(
  obj,
  method = "mtp",
  value = NULL,
  point = NULL,
  field_occurrence = "observed",
  format = "binary",
  return_threshold = FALSE,
  ...
)
```

## Arguments

- obj:

  A
  [BiodiversityScenario](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
  object to which an existing threshold is to be added.

- method:

  A specifc method for thresholding. See details for available options.

- value:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) value specifying
  the specific threshold for scenarios (Default: `NULL` grabs the value
  from `obj`).

- point:

  A [`sf`](https://r-spatial.github.io/sf/reference/sf.html) object
  containing observational data used for model training.

- field_occurrence:

  A [`character`](https://rdrr.io/r/base/character.html) location of
  biodiversity point records.

- format:

  [`character`](https://rdrr.io/r/base/character.html) indication of
  whether `"binary"`, `"normalize"` or `"percentile"` formatted
  thresholds are to be created (Default: `"binary"`). Also see
  Muscatello et al. (2021).

- return_threshold:

  Should threshold value be returned instead (Default: `FALSE`)

- ...:

  Any other parameter. Used to fetch value if set somehow.

## Value

A
[SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
if a
[SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
object as input. Otherwise the threshold is added to the respective
[`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
or
[`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
object.

## Details

The following options are currently implemented:

- `'fixed'` = applies a single pre-determined threshold. Requires
  `value` to be set.

- `'mtp'` = minimum training presence is used to find and set the lowest
  predicted suitability for any occurrence point.

- `'percentile'` = For a percentile threshold. A `value` as parameter
  has to be set here.

- `'min.cv'` = Threshold the raster so to minimize the coefficient of
  variation (cv) of the posterior. Uses the lowest tercile of the cv in
  space. Only feasible with Bayesian engines.

- `'TSS'` = Determines the optimal TSS (True Skill Statistic). Requires
  the `"modEvA"` package to be installed.

- `'kappa'` = Determines the optimal kappa value (Kappa). Requires the
  `"modEvA"` package to be installed.

- `'F1score'` = Determines the optimal F1score (also known as Sorensen
  similarity). Requires the `"modEvA"` package to be installed.

- `'Sensitivity'` = Determines the optimal sensitivity of presence
  records. Requires the `"modEvA"` package to be installed.

- `'Specificity'` = Determines the optimal specificity of presence
  records. Requires the `"modEvA"` package to be installed.

- `'AUC'` = Determines the optimal AUC of presence records. Requires the
  `"modEvA"` package to be installed.

- `'kmeans'` = Determines a threshold based on a 2 cluster k-means
  clustering. The presence class is assumed to be the cluster with the
  larger mean.

## References

- Lawson, C.R., Hodgson, J.A., Wilson, R.J., Richards, S.A., 2014.
  Prevalence, thresholds and the performance of presence-absence models.
  Methods Ecol. Evol. 5, 54–64. https://doi.org/10.1111/2041-210X.12123

- Liu, C., White, M., Newell, G., 2013. Selecting thresholds for the
  prediction of species occurrence with presence-only data. J. Biogeogr.
  40, 778–789. https://doi.org/10.1111/jbi.12058

- Muscatello, A., Elith, J., Kujala, H., 2021. How decisions about
  fitting species distribution models affect conservation outcomes.
  Conserv. Biol. 35, 1309–1320. https://doi.org/10.1111/cobi.13669

## See also

`"modEvA"`

## Examples

``` r
if (FALSE) { # \dontrun{
 # Where mod is an estimated DistributionModel
 tr <- threshold(mod)
 tr$plot_threshold()
} # }
```
