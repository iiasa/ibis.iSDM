# Set the parallel processing flag to TRUE

Small helper function to enable parallel processing. If set to `TRUE`,
then parallel inference (if supported by engines) and projection is
enabled across the package. For enabling prediction support beyond
sequential prediction see the
[`ibis_future`](https://iiasa.github.io/ibis.iSDM/reference/ibis_future.md)
function.

## Usage

``` r
ibis_enable_parallel()
```

## Value

Invisible

## See also

[future](https://future.futureverse.org/reference/future.html),
[ibis_future](https://iiasa.github.io/ibis.iSDM/reference/ibis_future.md)
