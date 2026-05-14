# Add pseudo-absence points to a point data set

For most engines, background or pseudo-absence points are necessary. The
distinction lies in how the absence data are handled. For
[`poisson`](https://rdrr.io/r/stats/family.html) distributed responses,
absence points are considered background points over which the intensity
of sampling (`lambda`) is integrated (in a classical Poisson
point-process model).

In contrast in [`binomial`](https://rdrr.io/r/stats/family.html)
distributed responses, the absence information is assumed to be an
adequate representation of the true absences and treated by the model as
such... Here it is advised to specify absence points in a way that they
represent potential true absence, such as for example through targeted
background sampling or by sampling them within/outside a given range.

## Usage

``` r
add_pseudoabsence(
  df,
  field_occurrence = "observed",
  template = NULL,
  settings = getOption("ibis.pseudoabsence")
)
```

## Arguments

- df:

  A [`sf`](https://r-spatial.github.io/sf/reference/sf.html),
  [`data.frame`](https://rdrr.io/r/base/data.frame.html) or
  [`tibble::tibble`](https://tibble.tidyverse.org/reference/tibble.html)
  object containing point data.

- field_occurrence:

  A [`character`](https://rdrr.io/r/base/character.html) name of the
  column containing the presence information (Default: `observed`).

- template:

  A
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object that is aligned with the predictors (Default: `NULL`). If set
  to `NULL`, then `background` in the
  [`pseudoabs_settings()`](https://iiasa.github.io/ibis.iSDM/reference/pseudoabs_settings.md)
  has to be a
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object.

- settings:

  A
  [`pseudoabs_settings()`](https://iiasa.github.io/ibis.iSDM/reference/pseudoabs_settings.md)
  objects. Absence settings are taken from
  [ibis_options](https://iiasa.github.io/ibis.iSDM/reference/ibis_options.md)
  otherwise (Default).

## Value

A [`data.frame`](https://rdrr.io/r/base/data.frame.html) containing the
newly created pseudo absence points.

## Details

A
[`pseudoabs_settings()`](https://iiasa.github.io/ibis.iSDM/reference/pseudoabs_settings.md)
object can be added to setup how absence points should be sampled. A
`bias` parameter can be set to specify a bias layer to sample from, for
instance a layer of accessibility. Note that when modelling several
datasets, it might make sense to check across all datasets whether
certain areas are truly absent. By default, the pseudo-absence points
are not sampled in areas in which there are already presence points.

## Note

This method removes all columns from the input `df` object other than
the `field_occurrence` column and the coordinate columns (which will be
created if not already present).

## References

- Stolar, J., & Nielsen, S. E. (2015). Accounting for spatially biased
  sampling effort in presence‐only species distribution modelling.
  Diversity and Distributions, 21(5), 595-608.

- Bird, T.J., Bates, A.E., Lefcheck, J.S., Hill, N.A., Thomson, R.J.,
  Edgar, G.J., Stuart-Smith, R.D., Wotherspoon, S., Krkosek, M.,
  Stuart-Smith, J.F. and Pecl, G.T., 2014. Statistical solutions for
  error and bias in global citizen science datasets. Biological
  Conservation, 173, pp.144-154.
