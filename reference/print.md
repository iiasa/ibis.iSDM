# Print

Display information about any object created through the ibis.iSDM
R-package.

## Usage

``` r
# S3 method for class 'distribution'
print(x, ...)

# S3 method for class 'BiodiversityDistribution'
print(x, ...)

# S3 method for class 'BiodiversityDatasetCollection'
print(x, ...)

# S3 method for class 'BiodiversityDataset'
print(x, ...)

# S3 method for class 'PredictorDataset'
print(x, ...)

# S3 method for class 'DistributionModel'
print(x, ...)

# S3 method for class 'BiodiversityScenario'
print(x, ...)

# S3 method for class 'Prior'
print(x, ...)

# S3 method for class 'PriorList'
print(x, ...)

# S3 method for class 'Engine'
print(x, ...)

# S3 method for class 'Settings'
print(x, ...)

# S3 method for class 'Log'
print(x, ...)

# S3 method for class 'Id'
print(x, ...)

# S4 method for class 'Id'
print(x, ...)
```

## Arguments

- x:

  Any object created through the package.

- ...:

  not used.

## Value

Object specific.

## See also

[`base::print()`](https://rdrr.io/r/base/print.html).

## Examples

``` r
if (FALSE) { # \dontrun{
# Where mod is fitted object
mod
print(mod)
} # }
```
