#' Plot wrappers
#'
#' @description Plots information from a given object where a plotting object is
#' available.
#'
#' @param x Any object belonging to [DistributionModel], [BiodiversityDatasetCollection],
#' [PredictorDataset] or [BiodiversityScenario].
#' @param what In case a [SpatRaster] is supplied, this parameter specifies the layer
#' to be shown (Default: \code{"mean"}).
#' @param ... Further arguments passed on to \code{x$plot}.
#'
#' @details The plotted outputs vary depending on what object is being plotted.
#' For example for a fitted [DistributionModel] the output is usually the fitted
#' spatial prediction (Default: \code{'mean'}).
#'
#' @return Graphical output
#'
#' @keywords misc
#'
#' @examples
#' \dontrun{
#' # Build and train a model
#' mod <- distribution(background) |>
#'   add_biodiversity_poipo(species) |>
#'   add_predictors(predictors) |>
#'   engine_glmnet() |>
#'   train()
#' # Plot the resulting model
#' plot(mod)
#' }
#'
#' @name plot
NULL

#' @rdname plot
#' @method plot DistributionModel
#' @export
plot.DistributionModel <- function(x, what = "mean", ...) x$plot(what,...)

#' @rdname plot
#' @method plot BiodiversityDatasetCollection
#' @export
plot.BiodiversityDatasetCollection <- function(x, ...) x$plot(...)

#' @rdname plot
#' @method plot PredictorDataset
#' @export
plot.PredictorDataset <- function(x, ...) x$plot(...)

#' @rdname plot
#' @method plot Engine
#' @export
plot.Engine <- function(x,...) x$plot(...)

#' @rdname plot
#' @method plot BiodiversityScenario
#' @export
plot.BiodiversityScenario <- function(x,...) x$plot(...)

#' Bivariate plot wrapper for distribution objects
#'
#' @description Often there is an intention to display not only the predictions
#' made with a SDM, but also the uncertainty of the prediction. Uncertainty be
#' estimated either directly by the model or by calculating the variation in
#' prediction values among a set of models.
#'
#' In particular Bayesian engines can produce not only mean estimates of fitted
#' responses, but also pixel-based estimates of uncertainty from the posterior
#' such as the standard deviation (SD) or the coefficient of variation of a
#' given prediction.
#'
#' This function makes use of the \code{"biscale"} R-package to create bivariate
#' plots of the fitted distribution object, allowing to visualize two variables
#' at once. It is mostly thought of as a convenience function to create such
#' bivariate plots for quick visualization.
#'
#' Supported Inputs are either single trained Bayesian [`DistributionModel`]
#' with uncertainty or the output of an [`ensemble()`] call. In both cases,
#' users have to make sure that \code{"xvar"} and \code{"yvar"} are set
#' accordingly.
#'
#' @param mod A trained [`DistributionModel`] or alternatively a [`SpatRaster`]
#' object with \code{prediction} model within.
#' @param xvar A [`character`] denoting the value on the x-axis (Default: \code{'mean'}).
#' @param yvar A [`character`] denoting the value on the y-axis (Default: \code{'sd'}).
#' @param plot A [`logical`] indication of whether the result is to be plotted
#' (Default: \code{TRUE})?
#' @param fname A [`character`] specifying the output filename a created figure
#' should be written to.
#' @param title Allows to respecify the title through a [`character`] (Default:\code{NULL}).
#' @param col A [`character`] stating the colour palette to use. Has to be either
#' a predefined value or a vector of colours. See \code{"biscale::bi_pal_manual"}.
#' Default: \code{"BlueGold"}.
#' @param ... Other engine specific parameters.
#'
#' @note
#' **This function requires the biscale package to be installed.**
#' Although a work around without the package could be developed, it was not
#' deemed necessary at this point. See also this
#' [gist](https://gist.github.com/scbrown86/2779137a9378df7b60afd23e0c45c188).
#'
#' @return Saved bivariate plot in \code{'fname'} if specified, otherwise plot.
#'
#' @seealso [partial], [plot.DistributionModel]
#' @keywords misc
#'
#' @export
#' @name bivplot
NULL

#' @rdname bivplot
#' @export
methods::setGeneric(
  "bivplot",
  signature = methods::signature("mod"),
  function(mod, xvar = "mean", yvar = "sd", plot = TRUE, fname = NULL, title = NULL, col = "BlueGold",...) standardGeneric("bivplot"))

#' @rdname bivplot
methods::setMethod(
  "bivplot",
  methods::signature(mod = "ANY"),
  function(mod, xvar = "mean", yvar = "sd", plot = TRUE, fname = NULL, title = NULL, col = "BlueGold",...) {
    # Generic checks
    assertthat::assert_that(is.logical(plot),
                            is.character(xvar),
                            is.character(yvar),
                            is.character(col) || is.vector(col),
                            is.null(title) || is.character(title),
                            is.null(fname) || is.character(fname),
                            isTRUE(plot) || is.character(fname)
    )
    # Check whether object is a raster, otherwise extract object
    if(is.Raster(mod)){
      assertthat::assert_that(terra::nlyr(mod)>1)
      obj <- mod
      # If number of layers equal to 2 (output from ensemble?), change xvar and yvar
      if(terra::nlyr(mod)==2 && !(xvar %in% names(obj))){
        if(getOption('ibis.setupmessages')) myLog('[Parameter]','yellow','Variable not found. Changing to layer names...')
        xvar <- names(obj)[1]; yvar <- names(obj)[2]
      }
    } else {
      assertthat::assert_that(inherits(mod, "DistributionModel"),
                              msg = "The bivplot function currently only works with fitted distribution objects!")
      # Check that distribution object has a prediction
      assertthat::assert_that("prediction" %in% mod$show_rasters(),
                              is.Raster(mod$get_data()),
                              msg = "No prediction found in the provided object.")
      obj <- mod$get_data()
    }

    # Check that at least mean and standard deviation is available
    assertthat::assert_that(xvar %in% names(obj),
                            yvar %in% names(obj),
                            msg = "Specified (default?) variables for xvar/yvar not found in model!")

    # Check whether biscale package is available
    check_package('biscale')
    check_package("cowplot")
    if(!("biscale" %in% loadedNamespaces()) || ('biscale' %notin% utils::sessionInfo()$otherPkgs) ) {
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
      terra::as.data.frame(xy = TRUE)
    names(df)[3:4] <- c("var1", "var2")
    suppressWarnings(
      df <- biscale::bi_class(df, x = var1, y = var2, dim = 3, style = "quantile")
    )

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
