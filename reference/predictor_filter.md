# Filter a set of correlated predictors to fewer ones

This function helps to remove highly correlated variables from a set of
predictors. It supports multiple options some of which require both
environmental predictors and observations, others only predictors.

Some of the options require different packages to be pre-installed, such
as `ranger` or `Boruta`.

## Usage

``` r
predictor_filter(env, keep = NULL, method = "pearson", ...)
```

## Arguments

- env:

  A
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  or alternatively
  [`data.frame`](https://rdrr.io/r/base/data.frame.html) or
  [`matrix`](https://rdrr.io/r/base/matrix.html) with extracted
  environmental covariates for a given species.

- keep:

  A [`vector`](https://rdrr.io/r/base/vector.html) with variables to
  keep regardless. These are usually variables for which prior
  information is known.

- method:

  Which method to use for constructing the correlation matrix (Options:
  `'pearson'` (Default), `'spearman'`\| `'kendal'`), `"abess"`, or
  `"boruta"`.

- ...:

  Other options for a specific method

## Value

A [`character`](https://rdrr.io/r/base/character.html)
[`vector`](https://rdrr.io/r/base/vector.html) of variable names to be
excluded. If the function fails due to some reason return `NULL`.

## Details

Available options are:

- `"none"` No prior variable removal is performed (Default).

- `"pearson"`, `"spearman"` or `"kendall"` Makes use of pairwise
  comparisons to identify and remove highly collinear predictors
  (Pearson's `r >= 0.7`).

- `"abess"` A-priori adaptive best subset selection of covariates via
  the `abess` package (see References). Note that this effectively fits
  a separate generalized linear model to reduce the number of
  covariates.

- `"boruta"` Uses the `Boruta` package to identify non-informative
  features.

## Note

Using this function on predictors effectively means that a separate
model is fitted on the data with all the assumptions that come with in
(e.g. linearity, appropriateness of response, normality, etc).

## Examples

``` r
if (FALSE) { # \dontrun{
 # Remove highly correlated predictors
 env <- predictor_filter(env, option = "pearson")
} # }
```
