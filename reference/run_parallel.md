# Parallel computation of function

Some computations take considerable amount of time to execute. This
function provides a helper wrapper for running functions of the
[`apply`](https://rdrr.io/r/base/apply.html) family to specified
outputs.

## Usage

``` r
run_parallel(
  X,
  FUN,
  cores = 1,
  approach = "future",
  export_packages = NULL,
  ...
)
```

## Arguments

- X:

  A [`list`](https://rdrr.io/r/base/list.html),
  [`data.frame`](https://rdrr.io/r/base/data.frame.html) or
  [`matrix`](https://rdrr.io/r/base/matrix.html) object to be fed to a
  single core or parallel [apply](https://rdrr.io/r/base/apply.html)
  call.

- FUN:

  A [`function`](https://rdrr.io/r/base/function.html) passed on for
  computation.

- cores:

  A [numeric](https://rdrr.io/r/base/numeric.html) of the number of
  cores to use (Default: `1`).

- approach:

  [`character`](https://rdrr.io/r/base/character.html) for the
  parallelization approach taken (Options: `"parallel"` or `"future"`).

- export_packages:

  A [`vector`](https://rdrr.io/r/base/vector.html) with packages to
  export for use on parallel nodes (Default: `NULL`).

- ...:

  Any other parameter passed on.

## Value

A [`list`](https://rdrr.io/r/base/list.html) with one result per element
or row in `X`.

## Details

By default, the parallel package is used for parallel computation,
however an option exists to use the
[future](https://future.futureverse.org/reference/future.html) package
instead.

## Examples

``` r
run_parallel(list(1:3, 4:6), mean, cores = 1)
#> The run_parallel function is likely deprecated and is only kept for reference...
#> [[1]]
#> [1] 2
#> 
#> [[2]]
#> [1] 5
#> 
```
