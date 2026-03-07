# Internal function to enable (a)synchronous parallel processing

This function checks if parallel processing can be set up and enables
it. **Ideally this is done by the user for more control!** In the
package parallelization is usually only used for predictions and
projections, but not for inference in which case parallel inference
should be handled by the engine.

## Usage

``` r
ibis_future(
  plan_exists = FALSE,
  cores = getOption("ibis.nthread", default = 2),
  strategy = getOption("ibis.futurestrategy"),
  workers = NULL
)
```

## Arguments

- plan_exists:

  A [`logical`](https://rdrr.io/r/base/logical.html) check on whether an
  existing
  [`future::future`](https://future.futureverse.org/reference/future.html)
  plan exists (Default: `FALSE`).

- cores:

  A [`numeric`](https://rdrr.io/r/base/numeric.html) number stating the
  number of cores to use.

- strategy:

  A [`character`](https://rdrr.io/r/base/character.html) denoting the
  strategy to be used for future. See help of
  [`future::future`](https://future.futureverse.org/reference/future.html)
  for options. (Default: `"multisession"`).

- workers:

  An optional list of remote machines or workers, e.g.
  `"c(remote.server.org)"`. Alternatively a `"cluster"` object can be
  provided.

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

## Note

The `'plan'` set by
[future::future](https://future.futureverse.org/reference/future.html)
exists after the function has been executed.

If the aim is to parallize across many species, this is better done in a
scripted solution. Make sure not to parallize predictions within
existing clusters to avoid out-of-memory issues.

## See also

[future::future](https://future.futureverse.org/reference/future.html)

## Examples

``` r
if (FALSE) { # \dontrun{
# Starts future job. F in this case is a prediction function.
ibis_future(cores = 4, strategy = "multisession")
} # }
```
