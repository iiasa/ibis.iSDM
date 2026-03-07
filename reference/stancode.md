# Show the stan code from a trained model

This helper function shows the code from a trained
[DistributionModel](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
using the
[`engine_stan`](https://iiasa.github.io/ibis.iSDM/reference/engine_stan.md).
This function is emulated after a similar functionality in the brms
R-package. **It only works with models inferred with stan!**

## Usage

``` r
stancode(obj, ...)

stancode.DistributionModel(obj, ...)
```

## Arguments

- obj:

  Any prepared object.

- ...:

  not used.

## Value

None.

## See also

rstan, cmdstanr, brms
