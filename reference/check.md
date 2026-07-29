# Check objects in the package for common errors or issues

Not always is there enough data or sufficient information to robustly
infer the suitable habitat or niche of a species. As many SDM algorithms
are essentially regression models, similar assumptions about model
convergence, homogeneity of residuals and inference usually apply
(although often ignored). This function simply checks the respective
input object for common issues or mistakes.

## Usage

``` r
check(obj, stoponwarning = FALSE)

# S4 method for class 'ANY'
check(obj, stoponwarning = FALSE)
```

## Arguments

- obj:

  A
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md),
  [`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
  or
  [`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
  object.

- stoponwarning:

  [`logical`](https://rdrr.io/r/base/logical.html) Should check return a
  stop if warning is raised? (Default: `FALSE`).

## Value

Message outputs

## Details

Different checks are implemented depending on the supplied object

- [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md)

&nbsp;

- Checks if there are less than 200 observations

- TODO: Add rm_insufficient_covs link

&nbsp;

- [`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)

&nbsp;

- Check model convergence

- Check if model is found

- Check if coefficients exist

- Check if there are unusual outliers in prediction (using 10 median
  absolute deviation)

- Check if threshold is larger than layer

&nbsp;

- [`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)

&nbsp;

- Check if scenario layers are valid

## Note

This function will likely be expanded with additional checks in the
future. If you have ideas, please let them know per issue.

## Examples

``` r
if (FALSE) { # \dontrun{
 # Where mod is an estimated DistributionModel
 check(mod)
} # }
```
