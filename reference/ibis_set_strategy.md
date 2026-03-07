# Set the number of threads for parallel processing.

Small helper function to respecify the strategy for parallel processing
(Default: `'sequential'`).

## Usage

``` r
ibis_set_strategy(strategy = "sequential")
```

## Arguments

- strategy:

  A [`character`](https://rdrr.io/r/base/character.html) with the
  strategy.

## Value

Invisible

## Details

Currently supported strategies are:

- `"sequential"` = Resolves futures sequentially in the current R
  process (Package default).

- `"multisession"` = Resolves futures asynchronously across `'cores'`
  sessions.

- `"multicore"` = Resolves futures asynchronously across on forked
  processes. Only works on UNIX systems!

- `"cluster"` = Resolves futures asynchronously in sessions on this or
  more machines.

- `"slurm"` = To be implemented: Slurm linkage via batchtools.

## See also

[future::future](https://future.futureverse.org/reference/future.html),
[ibis_future](https://iiasa.github.io/ibis.iSDM/reference/ibis_future.md)
