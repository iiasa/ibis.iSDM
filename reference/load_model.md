# Load a pre-computed model

The `load_model` function (opposed to the `write_model`) loads previous
saved
[`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md).
It is essentially a wrapper to
[`readRDS`](https://rspatial.github.io/terra/reference/serialize.html).

When models are loaded, they are briefly checked for their validity and
presence of necessary components.

## Usage

``` r
load_model(fname, verbose = getOption("ibis.setupmessages", default = TRUE))

# S4 method for class 'character'
load_model(fname, verbose = getOption("ibis.setupmessages", default = TRUE))
```

## Arguments

- fname:

  A [`character`](https://rdrr.io/r/base/character.html) depicting an
  output filename.

- verbose:

  [`logical`](https://rdrr.io/r/base/logical.html) indicating whether
  messages should be shown. Overwrites `getOption("ibis.setupmessages")`
  (Default: `TRUE`).

## Value

A
[`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
object.

## See also

write_model

## Examples

``` r
if (FALSE) { # \dontrun{
# Load model
mod <- load_model("testmodel.rds")

summary(mod)
} # }
```
