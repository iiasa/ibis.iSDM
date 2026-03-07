# Specify elevational preferences as offset

This function implements the elevation preferences offset defined in
Ellis‐Soto et al. (2021). The code here was adapted from the Supporting
materials script.

## Usage

``` r
add_offset_elevation(x, elev, pref, rate = 0.0089, add = TRUE)

# S4 method for class 'BiodiversityDistribution,SpatRaster,numeric'
add_offset_elevation(x, elev, pref, rate = 0.0089, add = TRUE)
```

## Arguments

- x:

  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object.

- elev:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  with the elevation for a given background.

- pref:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) vector of length
  `2` giving the lower and upper bound of known elevational preferences.
  Can be set to `Inf` if unknown.

- rate:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) for the rate used
  in the offset (Default: `.0089`). This parameter specifies the decay
  to near zero probability at elevation above and below the expert
  limits.

- add:

  [`logical`](https://rdrr.io/r/base/logical.html) specifying whether
  new offset is to be added. Setting this parameter to `FALSE` replaces
  the current offsets with the new one (Default: `TRUE`).

## Value

Adds a elevational offset to a
[`distribution`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
object.

## Details

Specifically this functions calculates a continuous decay and decreasing
probability of a species to occur from elevation limits. It requires a
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
with elevation information. A generalized logistic transform (aka
Richard's curve) is used to calculate decay from the suitable
elevational areas, with the `"rate"` parameter allowing to vary the
steepness of decline.

Note that all offsets created by this function are by default
log-transformed before export. In addition this function also
mean-centers the output as recommended by Ellis-Soto et al.

## References

- Ellis‐Soto, D., Merow, C., Amatulli, G., Parra, J.L., Jetz, W., 2021.
  Continental‐scale 1 km hummingbird diversity derived from fusing point
  records with lateral and elevational expert information. Ecography
  (Cop.). 44, 640–652. https://doi.org/10.1111/ecog.05119

- Merow, C., Allen, J.M., Aiello-Lammens, M., Silander, J.A., 2016.
  Improving niche and range estimates with Maxent and point process
  models by integrating spatially explicit information. Glob. Ecol.
  Biogeogr. 25, 1022–1036. https://doi.org/10.1111/geb.12453

## See also

Other offset:
[`add_offset()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset.md),
[`add_offset_bias()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset_bias.md),
[`add_offset_range()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset_range.md),
[`rm_offset()`](https://iiasa.github.io/ibis.iSDM/reference/rm_offset.md)

## Examples

``` r
if (FALSE) { # \dontrun{
 # Adds the offset to a distribution object
 distribution(background) |> add_offset_elevation(dem, pref = c(400, 1200))
} # }
```
