#' @include utils.R
NULL

#' Plot wrappers
#'
#' @description
#' Plots information from a given object where a plotting object is available.
#'
#' @param x Any object belonging to [DistributionModel], [BiodiversityDatasetCollection], [PredictorDataset] or [BiodiversityScenario].
#' @param what In case a [RasterLayer] is supplied, this parameter specifies the layer to be shown (Default: \code{"mean"}).
#' @param ... Further arguments passed on to \code{x$plot}.
#'
#' @details
#' The plotted outputs vary depending on what object is being plotted.
#' For example for a fitted [DistributionModel] the output is usually the fitted spatial prediction (Default: \code{'mean'}).
#' @return Graphical output
#' @keywords misc
#' @name plot
NULL

#' @rdname plot
#' @method plot DistributionModel
#' @keywords misc
#' @export
plot.DistributionModel <- function(x, ...) x$plot(...)

#' @rdname plot
#' @method plot BiodiversityDatasetCollection
#' @keywords misc
#' @export
plot.BiodiversityDatasetCollection <- function(x, ...) x$plot(...)

#' @rdname plot
#' @method plot PredictorDataset
#' @keywords misc
#' @export
plot.PredictorDataset <- function(x, ...) x$plot(...)

#' @rdname plot
#' @method plot Engine
#' @keywords misc
#' @export
plot.Engine <- function(x,...) x$plot(...)

#' @rdname plot
#' @method plot BiodiversityScenario
#' @keywords misc
#' @export
plot.BiodiversityScenario <- function(x,...) x$plot(...)


# ------------ '
#' Bivariate plot wrapper for distribution objects
#'
#' @description
#' In particular Bayesian engines can produce not only mean estimates of a fitted response,
#' but also estimates of uncertainty such as the standard deviation of the posterio or the coefficient
#' of variation of a given prediction.
#'
#' This function makes use of the [`biscale`] R-package to create bivariate plots of the fitted distribution object.
#' It is mostly thought of as a convenience function to create such plots for quick visualization.
#'
#' @param mod A trained `DistributionModel` object with \code{fit_best} model within.
#' @param xvar A [`character`] denoting the value on the x-axis (Default: \code{'mean'}).
#' @param yvar A [`character`] denoting the value on the y-axis (Default: \code{'sd'}).
#' @param plot A [`logical`] indication of whether the result is to be plotted?
#' @param fname A [`character`] specifying the output filename a created figure should be written to.
#' @param title Allows to respecify the title through a [`character`] (Default: \code{NULL}).
#' @param col A [`character`] stating the colour palette to use. Has to be either a predefined value or a
#' vector of colours. See \code{"biscale::bi_pal_manual"}. Default: \code{"BlueGold"}.
#' @param ... Other engine specific parameters
#' @seealso [partial], [plot.DistributionModel]
#' @note
#' **This function requires the biscale package to be installed.**
#' Although a work around without the package could be developed, it was not deemed necessary at this point.
#' See also this [gist](https://gist.github.com/scbrown86/2779137a9378df7b60afd23e0c45c188).
#' @return A [RasterLayer] with the created partial response.
#' @keywords misc
#' @export
#' @name bivplot
methods::setGeneric(
  "bivplot",
  signature = methods::signature("mod"),
  function(mod, xvar = "mean", yvar = "sd", plot = TRUE, fname = NULL, col = "BlueGold",...) standardGeneric("bivplot"))

#' @name bivplot
#' @rdname bivplot
#' @usage \S4method{bivplot}{ANY}(mod)
methods::setMethod(
  "bivplot",
  methods::signature(mod = "ANY"),
  function(mod, xvar = "mean", yvar = "sd", plot = TRUE, fname = NULL, title = NULL, col = "BlueGold",...) {
    assertthat::assert_that(inherits(mod, "DistributionModel"),
                            msg = "The bivplot function currently only works with fitted distribution objects!")
    # Generic checks
    assertthat::assert_that(is.logical(plot),
                            is.character(xvar),
                            is.character(yvar),
                            is.character(col) || is.vector(col),
                            is.null(title) || is.character(title),
                            is.null(fname) || is.character(fname)
    )
    # Check that distribution object has a prediction
    assertthat::assert_that("prediction" %in% mod$show_rasters(),
                            is.Raster(mod$get_data()),
                            msg = "No prediction found in the provided object.")
    obj <- mod$get_data()
    # Check that at least mean and standard deviation is available
    assertthat::assert_that(xvar %in% names(obj),
                            yvar %in% names(obj),
                            msg = "Specified (default?) variables for xvar/yvar not found in model!")

    # Check whether biscale package is available
    check_package('biscale')
    check_package('ggplot2')
    check_package("cowplot")
    if(!("biscale" %in% loadedNamespaces()) || ('biscale' %notin% sessionInfo()$otherPkgs) ) {
      try({requireNamespace('biscale');attachNamespace("biscale")},silent = TRUE)
    }

    # Check provided colours
    if(is.character(col)){
      choices <- c("Bluegill", "BlueGold", "BlueOr", "BlueYl", "Brown",
                   "Brown2", "DkBlue","DkBlue2", "DkCyan","DkCyan2", "DkViolet",
                   "DkViolet2", "GrPink","GrPink2", "PinkGrn", "PurpleGrn", "PurpleOr")
      col <- match.arg(col, choices, several.ok = FALSE)
    }

    # Define default title
    if(is.null(title)){
      title <- paste("Bivariate plot of prediction\n (",mod$model$runname,')')
    }

    # Create dimensions
    legend <- biscale::bi_legend(pal = col,
                        dim = 3,
                        xlab = paste0("Larger ", xvar),
                        ylab = paste0("Larger ", yvar),
                        size = 16)

    # Create data for plotting
    df <- obj[[c(xvar,yvar)]] |> predictor_transform(option = "norm") |>
      raster::as.data.frame(xy = TRUE)
    names(df)[3:4] <- c("var1", "var2")
    # df <- subset(df, complete.cases(df))
    df <- biscale::bi_class(df, x = var1, y = var2, dim = 3, style = "quantile")

    map <- ggplot2::ggplot() +
      ggplot2::geom_raster(data = df , ggplot2::aes(x = x, y = y, fill = bi_class)) +
      biscale::bi_theme(base_size = 16) +
      biscale::bi_scale_fill(pal = col, na.value = "transparent") +
      # coord_quickmap() +
      ggplot2::labs(
        title = title,
        x = "",
        y = ""
      ) +
      ggplot2::theme(legend.position = "none")

    # Add legend with cowplot
    finalPlot <- cowplot::ggdraw() +
      cowplot::draw_plot(map, 0, 0, 1, 1) +
      cowplot::draw_plot(legend, 0.2, .65, 0.2, 0.2)

    # Print the plot
    if(plot){
      print(finalPlot)
    }
    if(is.character(fname)){
      cowplot::ggsave2(filename = fname, plot = finalPlot)
    }

    return(finalPlot)
  }
)

