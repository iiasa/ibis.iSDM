# Install ibis dependencies

Some of the dependencies (R-Packages) that ibis.iSDM relies on are by
intention not added to the Description of the file to keep the number of
mandatory dependencies small and enable the package to run even on
systems that might not have all libraries pre-installed.

This function provides a convenience wrapper to install those missing
dependencies as needed. It furthermore checks which packages require
updating and updates them as needed.

## Usage

``` r
ibis_dependencies(deps = getOption("ibis.dependencies"), update = TRUE)
```

## Arguments

- deps:

  A [`vector`](https://rdrr.io/r/base/vector.html) with the names of the
  packages to be installed (Default: `"ibis.dependencies"` in
  [`ibis_options`](https://iiasa.github.io/ibis.iSDM/reference/ibis_options.md)).

- update:

  A [`logical`](https://rdrr.io/r/base/logical.html) flag of whether all
  (installed) packages should also be checked for updates (Default:
  `TRUE`).

## Value

Nothing. Packages will be installed.

## Note

INLA is handled in a special way as it is not available via cran.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Install and update all dependencies
  ibis_dependencies()
} # }
```
