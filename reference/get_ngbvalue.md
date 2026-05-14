# Function to extract nearest neighbour predictor values of provided points

This function performs nearest neighbour matching between biodiversity
observations and independent predictors, and operates directly on
provided data.frames. **Note that despite being parallized this function
can be rather slow for large data volumes of data!**

## Usage

``` r
get_ngbvalue(
  coords,
  env,
  longlat = TRUE,
  field_space = c("x", "y"),
  cheap = FALSE,
  ...
)
```

## Arguments

- coords:

  A [`matrix`](https://rdrr.io/r/base/matrix.html),
  [`data.frame`](https://rdrr.io/r/base/data.frame.html) or
  [`sf`](https://r-spatial.github.io/sf/reference/sf.html) object.

- env:

  A [`data.frame`](https://rdrr.io/r/base/data.frame.html) object with
  the predictors.

- longlat:

  A [`logical`](https://rdrr.io/r/base/logical.html) variable indicating
  whether the projection is long-lat.

- field_space:

  A [`vector`](https://rdrr.io/r/base/vector.html) highlight the columns
  from which coordinates are to be extracted (Default: `c('x','y')`).

- cheap:

  A [`logical`](https://rdrr.io/r/base/logical.html) variable whether
  the dataset is considered to be large and faster computation could
  help.

- ...:

  other options.

## Value

A [`data.frame`](https://rdrr.io/r/base/data.frame.html) with the
extracted covariate data from each provided data point.

## Details

Nearest neighbour matching is done via the
[geodist](https://hypertidy.github.io/geodist/reference/geodist.html)
R-package
([`geodist::geodist`](https://hypertidy.github.io/geodist/reference/geodist.html)).

## Note

If multiple values are of equal distance during the nearest neighbour
check, then the results is by default averaged.

## References

- Mark Padgham and Michael D. Sumner (2021). geodist: Fast,
  Dependency-Free Geodesic Distance Calculations. R package version
  0.0.7. https://CRAN.R-project.org/package=geodist

## Examples

``` r
if (FALSE) { # \dontrun{
 # Create matchup table
tab <- get_ngbvalue( coords = coords, # Coordinates
                     env = env # Data.frame with covariates and coordinates
                  )
} # }
```
