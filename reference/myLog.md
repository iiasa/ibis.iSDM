# Custom messaging function for scripts

This functions prints a message with a custom header and colour.

## Usage

``` r
myLog(title = "[Processing]", col = "green", ...)
```

## Arguments

- title:

  The title in the log output

- col:

  A [`character`](https://rdrr.io/r/base/character.html) indicating the
  text colour to be used. Supported are `'green'` / `'yellow'` / `'red'`

- ...:

  Any additional outputs or words for display

## Examples

``` r
if (FALSE) { # \dontrun{
myLog("[Setup]", "red", "Some error occurred during data preparation.")
} # }
```
