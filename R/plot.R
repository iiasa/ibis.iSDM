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

#' Bivariate prediction plot for distribution objects
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
                            isTRUE(plot)
    )
    # Check whether object is a raster, otherwise extract object
    if(is.Raster(mod)){
      assertthat::assert_that(terra::nlyr(mod)>1,msg = "SpatRaster object has less than 2 layers. Use plot().")
      obj <- mod
      # If number of layers equal to 2 (output from ensemble?), change xvar and yvar
      if(terra::nlyr(mod)==2 && !(xvar %in% names(obj))){
        if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Parameter]','yellow','Variable not found. Changing to layer names...')
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
      if(is.Raster(mod)) tt <- "" else tt <- paste0("\n (",mod$model$runname,")")
      title <- paste("Bivariate plot of prediction",tt)
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

#' Niche plot for distribution objects
#'
#' @description
#' The suitability of any given area for a biodiversity feature can in
#' many instances be complex and non-linear. Visualizing obtained suitability
#' predictions (e.g. from [`train()`]) against underlying predictors might help
#' to explain the underlying gradients of the niche.
#'
#' Supported Inputs for this function are either single trained \code{ibis.iSDM}
#' [`DistributionModel`] objects or alternatively a set of three [`SpatRaster`] objects.
#' In both cases, users can specify \code{"xvar"} and \code{"yvar"} explicitly
#' or leave them empty. In the latter case a principal component analysis (PCA)
#' is conducted on the full environmental stack (loaded from [`DistributionModel`]
#' or supplied separately).
#'
#' @param mod A trained [`DistributionModel`] or alternatively a [`SpatRaster`]
#' object with \code{prediction} model within.
#' @param xvar A [`character`] denoting the predictor on the x-axis. Alternatively a [`SpatRaster`]
#' object can be provided.
#' @param yvar A [`character`] denoting the predictor on the y-axis. Alternatively a [`SpatRaster`]
#' object can be provided.
#' @param envvars A [`SpatRaster`] object containing all environmental variables. Only
#' used if \code{xvar} and \code{yvar} is empty (Default: \code{NULL}).
#' @param overlay_data A [`logical`] on whether training data should be overlaid
#' on the plot. Only used for [`DistributionModel`] objects (Default: \code{FALSE}).
#' @param plot A [`logical`] indication of whether the result is to be plotted
#' (Default: \code{TRUE})?
#' @param fname A [`character`] specifying the output file name a created figure
#' should be written to.
#' @param title Allows to respecify the title through a [`character`] (Default: \code{NULL}).
#' @param pal An optional [`vector`] with continuous custom colours (Default: \code{NULL}).
#' @param ... Other engine specific parameters.
#'
#' @return Saved niche plot in \code{'fname'} if specified, otherwise plot.
#'
#' @seealso [partial], [plot.DistributionModel]
#' @keywords misc
#' @examples
#' # Make quick prediction
#' background <- terra::rast(system.file('extdata/europegrid_50km.tif',
#' package='ibis.iSDM',mustWork = TRUE))
#' virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM'), 'points',quiet = TRUE)
#' ll <- list.files(system.file('extdata/predictors/',package = 'ibis.iSDM',mustWork = TRUE),full.names = TRUE)
#'
#' # Load them as rasters
#' predictors <- terra::rast(ll);names(predictors) <- tools::file_path_sans_ext(basename(ll))
#'
#' # Add GLM as an engine and predict
#' fit <- distribution(background) |>
#' add_biodiversity_poipo(virtual_points, field_occurrence = 'Observed',
#' name = 'Virtual points',docheck = FALSE) |>
#' add_predictors(predictors, transform = 'none',derivates = 'none') |>
#' engine_glm() |>
#' train()
#'
#' # Plot niche for prediction for temperature and forest cover
#' nicheplot(fit, xvar = "bio01_mean_50km", yvar = "CLC3_312_mean_50km" )
#' @export
#' @name nicheplot
NULL

#' @rdname nicheplot
#' @export
methods::setGeneric(
  "nicheplot",
  signature = methods::signature("mod"),
  function(mod, xvar = NULL, yvar = NULL, envvars = NULL, overlay_data = FALSE,
           plot = TRUE, fname = NULL, title = NULL, pal = NULL, ...) standardGeneric("nicheplot"))

#' @rdname nicheplot
methods::setMethod(
  "nicheplot",
  methods::signature(mod = "ANY"),
  function(mod, xvar = NULL, yvar = NULL, envvars = NULL, overlay_data = FALSE,
           plot = TRUE, fname = NULL, title = NULL, pal = NULL, ...) {
    # Generic checks
    assertthat::assert_that(
                            is.null(xvar) || (is.character(xvar) || is.Raster(xvar)),
                            is.null(yvar) || (is.character(yvar) || is.Raster(yvar)),
                            is.null(envvars) || is.Raster(envvars),
                            is.logical(overlay_data),
                            is.null(title) || is.character(title),
                            is.null(fname) || is.character(fname),
                            is.logical(plot),
                            is.null(pal) || is.vector(pal),
                            msg = "Provide correct parameters (see help file)."
    )
    check_package("ggplot2") # Should be loaded by default

    # Check specific on x/y variables
    assertthat::assert_that((is.null(xvar) && is.null(yvar)) ||
                              (is.character(xvar) && is.character(yvar)),
                            msg = "Both x and y variable names must be specified!")

    # Check whether object is a raster, otherwise extract object
    if(is.Raster(mod)){
      obj <- mod
      # Check x and y variables are correct
      assertthat::assert_that(
        is.Raster(xvar) && is.Raster(yvar),
        msg = "SpatRaster objects need to be supplied as xvar and yvar!"
      )

      # Check if set
      if(!is.null(xvar) && !is.null(yvar)){
        # Align if mismatching
        if(!is_comparable_raster(obj, xvar)){
          warning('xvariable not aligned with prediction. Aligning them now...')
          xvar <- alignRasters(xvar, obj, method = 'bilinear', func = mean, cl = FALSE)
        }
        if(!is_comparable_raster(obj, yvar)){
          warning('yvariable not aligned with prediction. Aligning them now...')
          yvar <- alignRasters(yvar, obj, method = 'bilinear', func = mean, cl = FALSE)
        }
      } else {
        assertthat::assert_that(is.Raster(envvars),
                                terra::nlyr(envvars)>1,
                                msg = "A multi layer environmental stack has to be supplied directly!")

        if(!is_comparable_raster(obj, envvars)){
          warning('Predictorstack not aligned with prediction. Aligning them now...')
          envvars <- alignRasters(envvars, obj, method = 'bilinear', func = mean, cl = FALSE)
        }
      }

    } else {
      # Distribution model objects #
      assertthat::assert_that(inherits(mod, "DistributionModel"),
                              is.Raster(mod$get_data()),
                              msg = "The nicheplot function currently only works with fitted distribution objects!")
      # Check that distribution object has a prediction
      assertthat::assert_that("prediction" %in% mod$show_rasters(),
                              is.Raster(mod$get_data()),
                              msg = "No prediction found in the provided object.")
      obj <- mod$get_data()[[1]] # Get the first layer

      # Check if set
      if(!is.null(xvar) && !is.null(yvar)){
        assertthat::assert_that(is.character(xvar), is.character(yvar),
                                msg = "Specify predictor names for both xvar and yvar.")
        # Check that variables are present
        assertthat::assert_that(
          xvar %in% mod$model$predictors_names, yvar %in% mod$model$predictors_names,
          msg = "Variables not used in underlying model?"
        )
        # Also get the xvar/yvar
        if(is.character(xvar)) xvar <- mod$model$predictors_object$get_data()[[xvar]]
        if(is.character(yvar)) yvar <- mod$model$predictors_object$get_data()[[yvar]]
      } else {
        if(is.null(envvars)){
          envvars <- mod$model$predictors_object$get_data()
        } else {
          assertthat::assert_that(is.Raster(envvars),
                                  terra::nlyr(envvars)>1,
                                  msg = "A multi layer environmental stack has to be supplied directly!")        }
      }
    }

    # Get training data if set
    if(overlay_data){
      assertthat::assert_that(inherits(mod, "DistributionModel"),
                              msg = "Data overlay currently only works with a fitted model object!")
      # Collect Biodiversity occurrence data
      occ <- collect_occurrencepoints(mod$model,
                                      include_absences = FALSE,
                                      addName = TRUE,tosf = TRUE)
    }

    # Define default title
    if(is.null(title)){
      if(is.Raster(mod)) tt <- names(obj) else tt <- paste0("\n (",mod$model$runname,")")
      title <- paste("Niche plot for prediction ",tt)
    }

    # Define colour palette if not set
    if(is.null(pal)){
      pal <- ibis_colours$sdm_colour
    } else {
      assertthat::assert_that(length(pal)>=1)
    }

    # ----------- #
    if(is.null(xvar) && is.null(yvar)){
      # No specific variables found. Conduct a principle component analysis
      # and predict for the first two axes.
      assertthat::assert_that(is.Raster(envvars))
      pca <- terra::prcomp(envvars, maxcell = 10000)
      rpca <- terra::predict(envvars, pca, index=1:2)

      # Now check number of cells and extract. If too large, sample at random
      o <- c(obj, rpca)
      names(o) <- c("mean", "PC1", "PC2")
      if(terra::ncell(o)>10000){
        # Messenger
        if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Visualization]','green','Sampling at random grid cells for extraction.')
        ex <- terra::spatSample(o, size = 10000, method = "random",
                                as.df = TRUE, na.rm = TRUE)
      } else {
        # Extract
        ex <- terra::as.data.frame(o, xy = FALSE, na.rm = TRUE, time = FALSE)
      }
      assertthat::assert_that(nrow(ex)>0)

      # Define variable names
      xvar_lab <- "PC1"
      yvar_lab <- "PC2"
      col_lab <- names(obj)

      # Now plot
      viz <- ggplot2::ggplot() +
        ggplot2::theme_bw(base_size = 20) +
        ggplot2::geom_point(data = ex, ggplot2::aes(x = PC1, y = PC2, colour = mean, alpha=mean)) +
        ggplot2::scale_colour_gradientn(colours = pal) +
        ggplot2::guides(colour = ggplot2::guide_colorbar(title = col_lab), alpha = "none") +
        ggplot2::theme(legend.position = "bottom",
                       legend.title = ggplot2::element_text(vjust = 1),
                       legend.key.size = ggplot2::unit(1, "cm")) +
        ggplot2::labs(
          title = title,
          x = xvar_lab,
          y = yvar_lab
        )

      # Should the training data be overlaid?
      if(overlay_data){
        pp <- terra::extract(o, occ, ID = FALSE)
        pp$name <- as.factor( occ$name )

        viz <- viz +
          ggplot2::geom_point(data = pp,
                              ggplot2::aes(x = PC1, y = PC2),
                              col = "black",
                              size = 1.5,show.legend = TRUE) +
          ggplot2::labs(subtitle = "Training data (black)")
      }

    } else {
      # Make plot for two variables which should have been collected above.
      # Check that all Raster objects are there
      assertthat::assert_that(
        is.Raster(xvar), is.Raster(yvar), is.Raster(obj),
        terra::hasValues(obj),
        msg = "Layers are not in spatial format?"
      )

      # Define variable names
      xvar_lab <- names(xvar)
      yvar_lab <- names(yvar)
      col_lab <- names(obj)

      # Now check number of cells and extract. If too large, sample at random
      o <- c(obj, xvar, yvar)
      names(o) <- c("mean", "xvar", "yvar")
      if(terra::ncell(o)>10000){
        # Messenger
        if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Visualization]','green','Sampling at random grid cells for extraction.')
        ex <- terra::spatSample(o, size = 10000, method = "random",
                                as.df = TRUE, na.rm = TRUE)
      } else {
        # Extract
        ex <- terra::as.data.frame(o, xy = FALSE, na.rm = TRUE, time = FALSE)
      }
      assertthat::assert_that(nrow(ex)>0)


      # Now plot
      viz <- ggplot2::ggplot() +
        ggplot2::theme_bw(base_size = 20) +
        ggplot2::geom_point(data = ex, ggplot2::aes(x = xvar, y = yvar, colour = mean, alpha = mean)) +
        ggplot2::scale_colour_gradientn(colours = pal) +
        ggplot2::guides(colour = ggplot2::guide_colorbar(title = col_lab), alpha = "none") +
        ggplot2::theme(legend.position = "bottom",
                       legend.title = ggplot2::element_text(vjust = 1),
                       legend.key.size = ggplot2::unit(1, "cm")) +
        ggplot2::labs(
          title = title,
          x = xvar_lab,
          y = yvar_lab
        )

      # Should the training data be overlaid?
      if(overlay_data){
        pp <- terra::extract(o, occ, ID = FALSE)
        pp$name <- as.factor( occ$name )

        viz <- viz +
          ggplot2::geom_point(data = pp,
                              ggplot2::aes(x = xvar, y = yvar),
                              col = "black",
                              size = 1.5,show.legend = TRUE) +
          ggplot2::labs(subtitle = "Training data (black)")
      }
    }

    # Print the plot
    if(plot){
      return(viz)
    }
    if(is.character(fname)){
      cowplot::ggsave2(filename = fname, plot = viz)
    }
  }
)
