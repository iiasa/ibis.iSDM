# Sanitize variable names

Prepared covariates often have special characters in their variable
names which can or can not be used in formulas or cause errors for
certain engines. This function converts special characters of variable
names into a format

## Usage

``` r
sanitize_names(names)
```

## Arguments

- names:

  A [`vector`](https://rdrr.io/r/base/vector.html) of
  [`character`](https://rdrr.io/r/base/character.html) vectors to be
  sanitized.

## Value

A [`vector`](https://rdrr.io/r/base/vector.html) of sanitized
[`character`](https://rdrr.io/r/base/character.html).

## Examples

``` r
# Correct variable names
vars <- c("Climate-temperature2015", "Elevation__sealevel", "Landuse.forest..meanshare")
sanitize_names(vars)
#> [1] "Climate_temperature2015"   "Elevation_sealevel"       
#> [3] "Landuse.forest..meanshare"
```
