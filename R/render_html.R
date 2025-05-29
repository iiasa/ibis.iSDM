#' render_html
#'
#' @description Renders [DistributionModel] to HTML
#'
#' @param mod Any object belonging to [DistributionModel]
#' @param file \code{Character} with path to file.
#' @param title \code{Character} with title of document.
#' @param author \code{Character} with name of author.
#' @param notes \code{Character} with notes added at the beginning of the document.
#' @param ... Currently not used
#'
#' @details Renders a HTML file with several summaries of a trained [DistributionModel].
#' The file paths must be an HTML file ending. The functions creates a temporary
#' Rmd file that gets renders as a HTML using the \code{file} argument.
#'
#' @return Writes HTML file
#'
#' @keywords misc
#'
#' @examples
#' \dontrun{
#' mod <- distribution(background) |>
#'   add_biodiversity_poipo(species) |>
#'   add_predictors(predictors) |>
#'   engine_glmnet() |>
#'   train()
#'
#' render_html(mod, file = "Test.html")
#' }
#' @export
#' @name render_html
NULL

#' @rdname render_html
#' @export
methods::setGeneric(
  "render_html",
  signature = methods::signature("mod"),
  function(mod, file, title = NULL, author = NULL, notes = "-", ...) standardGeneric("render_html"))

#' @rdname render_html
methods::setMethod(
  "render_html",
  methods::signature(mod = "ANY"),
  function(mod, file, title = NULL, author = NULL, notes = "-",...) {

    # create path for rmd file
    file_rmd <- paste0(tempfile(), ".rmd")

    assertthat::assert_that(inherits(mod, "DistributionModel"),
                            msg = "Please provide DistributionModel")

    if(!assertthat::has_extension(path = file, ext = "html")) cli::cli_abort("Please use .html file.",
                                                                   call. = FALSE)

    # get some meta information document
    if (is.null(author)) author <- Sys.getenv("USERNAME")
    if (is.null(title)) title <- "Overview SDM"
    date <- Sys.Date()

    # get meta information of the model
    engine <- mod$get_name()
    res <- mod$get_resolution()
    thres <- mod$get_thresholdtype()
    thres <- ifelse(test = inherits(x = thres, what = "Waiver"), yes = "none", no = thres)
    thres_value <- ifelse(test = thres == "none", yes = "none", no = mod$get_thresholdvalue())
    offsets <- mod$has_offset()

    # get fit
    coef <- mod$summary()
    eq <- mod$fits$fit_best_equation

    # calculate partial response
    part_res <- partial(mod = mod, plot = FALSE)

    # setup some objects
    col_sdm <- ibis_colours[['sdm_colour']]
    col_thres <- c("grey", "black")

    # The indentation has to be exactly as it is inside the cat
    cat("---
title: \'", title, "\'
author: \'", author, "\'
date:  \'", as.character(date), "\'
output:
  html_document:
    theme: simplex
---

<!-- set defaults -->
\`\`\`{r setup, include = FALSE}
knitr::opts_chunk$set(echo = TRUE)
\`\`\`

*This report was created automatically by the `ibis.iSDM` package.*

<br>

---

General notes: ", notes, "

---

<br>

## Model information

Used engine: ", engine, "

Resolution: x=", res[1], " y=", res[2], "

<!-- Threshold: ", thres, ", Threshold value: ", thres_value, " -->

Offsets?: ", offsets, "

## Prediction map
\`\`\`{r plot_pred, echo = FALSE, out.width = \'75%\', fig.align = \'center\'}
if (thres == \'none\') {
plot(mod$fits$prediction, main = \'Prediction\', col = col_sdm)} else {
par(mfrow = c(1,2))
plot(mod$fits$prediction, main = \'Prediction\', col = col_sdm)
plot(mod$fits[[4]], main = \'Threshold\', col = col_thres)
par(mfrow = c(1,1))}
\`\`\`

## Coefficients

The best fit used the following equation and fitted coefficients:
\`\`\`{r eq, echo = FALSE}
eq
\`\`\`

\`\`\`{r table_coef, echo = FALSE}
knitr::kable(coef)
\`\`\`

## Partial responses

\`\`\`{r plot_res, echo = FALSE, out.width = \'65%\', fig.align = \'center\'}
ggplot2::ggplot(data = part_res, ggplot2::aes(x = partial_effect)) +
ggplot2::theme_classic() +
ggplot2::geom_line(ggplot2::aes(y = mean)) +
ggplot2::facet_wrap(. ~ variable, scales = \'free\') +
ggplot2::labs(x = \'Variable\', y = \'Partial effect\')
\`\`\`", file = file_rmd, sep = "")

    # render_html file
    rmarkdown::render(input = file_rmd, output_file = file)

    # make sure the temp file gets deleted
    unlink(x = file_rmd)

  }
)
