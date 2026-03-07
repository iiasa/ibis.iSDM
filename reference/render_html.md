# render_html

Renders
[DistributionModel](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
to HTML

## Usage

``` r
render_html(mod, file, title = NULL, author = NULL, notes = "-", ...)

# S4 method for class 'ANY'
render_html(mod, file, title = NULL, author = NULL, notes = "-", ...)
```

## Arguments

- mod:

  Any object belonging to
  [DistributionModel](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)

- file:

  `Character` with path to file.

- title:

  `Character` with title of document.

- author:

  `Character` with name of author.

- notes:

  `Character` with notes added at the beginning of the document.

- ...:

  Currently not used

## Value

Writes HTML file

## Details

Renders a HTML file with several summaries of a trained
[DistributionModel](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md).
The file paths must be an HTML file ending. The function creates a
temporary Rmd file that gets rendered as HTML using the `file` argument.

## Examples

``` r
if (FALSE) { # \dontrun{
mod <- distribution(background) |>
  add_biodiversity_poipo(species) |>
  add_predictors(predictors) |>
  engine_glmnet() |>
  train()

render_html(mod, file = "Test.html")
} # }
```
