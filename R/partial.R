#' @include utils.R
NULL

#' Obtain partial effects of trained model
#'
#' @description Create a partial response or effect plot of a trained model.
#'
#' @param mod A trained `DistributionModel` object with \code{fit_best} model
#'   within.
#' @param x.var A [character] indicating the variable for which a partial effect
#'   is to be calculated.
#' @param constant A [numeric] constant to be inserted for all other variables.
#'   Default calculates a mean per variable.
#' @param variable_length [numeric] The interpolation depth (nr. of points) to
#'   be used (Default: \code{100}).
#' @param values [numeric] Directly specified values to compute partial effects
#'   for. If this parameter is set to anything other than \code{NULL}, the
#'   parameter \code{"variable_length"} is ignored (Default: \code{NULL}).
#' @param newdata An optional [data.frame] with provided data for partial
#'   estimation (Default: \code{NULL}).
#' @param plot A [`logical`] indication of whether the result is to be plotted?
#' @param type A specified type, either \code{'response'} or \code{'predictor'}.
#'   Can be missing.
#' @param ... Other engine specific parameters.
#' @seealso [partial]
#' @details By default the mean is calculated across all parameters that are not
#'   \code{x.var}. Instead a *constant* can be set (for instance \code{0}) to be
#'   applied to the output.
#' @return A [data.frame] with the created partial response.
#' @aliases partial
#' @examples
#' \dontrun{
#'  # Do a partial calculation of a trained model
#'  partial(fit, x.var = "Forest.cover", plot = TRUE)
#' }
#' @keywords partial
#' @export
#' @name partial
methods::setGeneric(
  "partial",
  signature = methods::signature("mod","x.var"),
  function(mod, x.var, constant = NULL, variable_length = 100, values = NULL, newdata = NULL, plot = FALSE, type = "response", ...) standardGeneric("partial"))

#' @name partial
#' @rdname partial
#' @usage
#'   \S4method{partial}{ANY,character,ANY,numeric,ANY,ANY,logical,character}(mod,x.var,constant,variable_length,values,newdata,plot,type,...)
methods::setMethod(
  "partial",
  methods::signature(mod = "ANY", x.var = "character"),
  function(mod, x.var, constant = NULL, variable_length = 100,
           values = NULL, newdata = NULL, plot = FALSE, type = "response",...) {
    assertthat::assert_that(!missing(x.var),msg = 'Specify a variable name in the model!')
    assertthat::assert_that(inherits(mod, "DistributionModel"),
                            is.character(x.var),
                            is.null(constant) || is.numeric(constant),
                            is.numeric(variable_length),
                            is.null(newdata) || is.data.frame(newdata),
                            is.numeric(values) || is.null(values),
                            is.logical(plot)
    )
    # Work around to call partial response directly
    if(inherits(mod,'DistributionModel')){
      partial.DistributionModel(mod, x.var, constant, variable_length, values, newdata, plot, type, ...)
    } else {
      stop('Partial response calculation not supported!')
    }
  }
)

#' @rdname partial
#' @method partial DistributionModel
#' @keywords partial
#' @export
partial.DistributionModel <- function(mod, ...) mod$partial(...)

#' Obtain spatial partial effects of trained model
#'
#' @description Similar as [partial] this function calculates a partial response
#' of a trained model for a given variable. Differently from [partial] in space.
#' However the result is a [`SpatRaster`] showing the spatial magnitude of the
#' partial response.
#' @param mod A [`DistributionModel-class`] object with trained model.
#' @param x.var A [character] indicating the variable for which a partial effect
#'   is to be calculated.
#' @param constant A [numeric] constant to be inserted for all other variables.
#'   Default calculates the [mean] per variable.
#' @param newdata A [`data.frame`] on which to calculate the spartial for. Can be
#' for example created from a raster file (Default: \code{NULL}).
#' @param plot A [logical] indication of whether the result is to be plotted?
#' @param ... Other engine specific parameters.
#' @seealso [partial]
#' @details By default the [mean] is calculated across all parameters that are
#'   not \code{x.var}. Instead a *constant* can be set (for instance \code{0})
#'   to be applied to the output.
#' @returns A [SpatRaster] containing the mapped partial response of the
#'   variable.
#' @aliases spartial
#' @examples
#' \dontrun{
#'  # Create and visualize the spartial effect
#'  spartial(fit, x.var = "Forest.cover", plot = TRUE)
#' }
#' @keywords partial
#' @export
#' @name spartial
methods::setGeneric(
  "spartial",
  signature = methods::signature("mod","x.var"),
  function(mod, x.var, constant = NULL, newdata = NULL, plot = FALSE, ...) standardGeneric("spartial"))

#' @name spartial
#' @rdname spartial
#' @usage
#'   \S4method{spartial}{ANY,character,ANY,ANY,logical}(mod,x.var,constant,newdata,plot,...)
methods::setMethod(
  "spartial",
  methods::signature(mod = "ANY", x.var = "character"),
  function(mod, x.var, constant = NULL, newdata = NULL, plot = FALSE, ...) {
    assertthat::assert_that(!missing(x.var),msg = 'Specify a variable name in the model!')
    assertthat::assert_that(inherits(mod, "DistributionModel"),
                            is.character(x.var),
                            is.null(constant) || is.numeric(constant),
                            is.null(newdata) || is.data.frame(newdata),
                            is.logical(plot)
    )
    # Work around to call partial response directly
    if(inherits(mod,'DistributionModel')){
      spartial.DistributionModel(mod, x.var, constant, newdata, plot, ...)
    } else {
      stop('Spatial partial response calculation not supported!')
    }
  }
)

#' @rdname spartial
#' @method spartial DistributionModel
#' @keywords partial
#' @export
spartial.DistributionModel <- function(mod, ...) mod$spartial(...)

#' Visualize the density of the data over the environmental data
#'
#' @description Based on a fitted model, plot the density of observations over
#' the estimated variable and environmental space. Opposed to the [partial] and
#' [spartial] functions, which are rather low-level interfaces, this function
#' provides more detail in the light of the data. It is also able to contrast
#' different variables against each other and show the used data.
#'
#' @details This functions calculates the observed density of presence and
#' absence points over the whole surface of a specific variable. It can be used
#' to visually inspect the fit of the model to data.
#'
#' @note By default all variables that are not \code{x.var} are hold constant at
#' the mean.
#' @param mod A trained `DistributionModel` object. Requires a fitted model and
#'   inferred prediction.
#' @param x.var A [character] indicating the variable to be investigated. Can be
#'   a [`vector`] of length \code{1} or \code{2}.
#' @param df [`logical`] if plotting data should be returned instead (Default:
#'   \code{FALSE}).
#' @param ... Other engine specific parameters.
#'
#' @seealso [partial]
#' @concept Visual style emulated from ENMTools package.
#' @references
#' * Warren, D.L., Matzke, N.J., Cardillo, M., Baumgartner, J.B., Beaumont, L.J., Turelli, M., Glor, R.E., Huron, N.A., Simões, M., Iglesias, T.L. Piquet, J.C., and Dinnage, R. 2021. ENMTools 1.0: an R package for comparative ecological biogeography. Ecography, 44(4), pp.504-511.
#' @returns A [`ggplot2`] object showing the marginal response in light of the
#'   data.
#' @aliases partial_density
#' @examples
#' \dontrun{
#'  # Do a partial calculation of a trained model
#'  partial_density(fit, x.var = "Forest.cover")
#'  # Or with two variables
#'  partial_density(fit, x.var = c("Forest.cover", "bio01"))
#' }
#' @keywords partial
#' @export
#' @name partial_density
methods::setGeneric(
  "partial_density",
  signature = methods::signature("mod","x.var"),
  function(mod, x.var, df = FALSE,...) standardGeneric("partial_density"))

#' @name partial_density
#' @rdname partial_density
#' @usage \S4method{partial_density}{ANY,character,logical}(mod,x.var,df,...)
methods::setMethod(
  "partial_density",
  methods::signature(mod = "ANY", x.var = "character"),
  function(mod, x.var, df = FALSE,...) {
    assertthat::assert_that(!missing(x.var),msg = 'Specify a variable name in the model!')
    assertthat::assert_that(inherits(mod, "DistributionModel"),
                            all(is.character(x.var)),
                            length(x.var)>=1, length(x.var) <=2,
                            is.logical(df)
    )
    # Check that mod has all the information necessary
    assertthat::assert_that(!is.Waiver(mod$get_data()),
                            !is.Waiver(mod$get_data("fit_best"))
                            )

    # Ensure that x.var variables are among the coefficients
    assertthat::assert_that(all(x.var %in% mod$get_coefficients()[[1]]),
                            msg = "Provided variables not found among the coefficients!")

    # --- #
    model <- mod$model # Get model
    pred <- mod$get_data() # Get prediction
    # Get target variables from the model object
    vars <- model$predictors_object$get_data()[[x.var]]
    assertthat::assert_that(is.Raster(pred), is.Raster(vars),
                            msg = "Either no prediction or predictor object found.")

    # Collect observations from model object
    obs <- collect_occurrencepoints(model = model,
                             include_absences = TRUE,
                             addName = TRUE,
                             tosf = FALSE)

    # Get variable bounds
    vars_lims <- terra::minmax(vars)
    if(any(is.na(vars_lims))){
      # Quick check in case range can not be correctly extracted
      vars <- terra::setMinMax(vars)
      vars_lims <- terra::minmax(vars)
    }

    # --- #
    normalise <- function(z) (z - min(z))/(max(z)-min(z))
    # Extract variables and make partial prediction
    out <- data.frame()
    for(v in x.var){
      # Create a partial prediction with the model
      pp <- partial.DistributionModel(mod, v,
                                      values = seq(vars_lims["min",v], vars_lims["max",v], length = 100),
                                      variable_length = 100, plot = FALSE, type = "response")
      # Normalize and save
      pp[['mean']] <- normalise(pp[['mean']])

      # For each data type
      for(ty in unique(obs$type)){
        sub <- subset(obs, type == ty)
        # First extract for each observation the predicted covariate
        ss <- guess_sf(sub)
        suppressWarnings(
          ss <- ss |> sf::st_set_crs(value = sf::st_crs(model$background))
        )
        if(is.factor(ss$observed)) ss$observed <- as.numeric(as.character(ss$observed))

        # Extract for all presence and absence points
        pres <- get_rastervalue(ss |> dplyr::filter(observed >= 1), vars, ngb_fill = TRUE)[,v]
        pres.dens <- density(pres,
                             from = vars_lims["min",v],
                             to = vars_lims["max",v],
                             n = 100,
                             na.rm = TRUE)$y
        # Save in object
        pp$presence.density <- pres.dens/max(pres.dens)

        # Absence
        if(any(ss$observed==0)){
          abs <- get_rastervalue(ss |> dplyr::filter(observed == 0), vars, ngb_fill = TRUE)[,v]
          abs.dens <- density(abs,
                               from = vars_lims["min",v],
                               to = vars_lims["max",v],
                               n = 100,
                               na.rm = TRUE)$y
          pp$absence.density <- abs.dens/max(abs.dens)
        } else {
          # FIXME: Maybe do a randomized background extraction?
          pp$absence.density <- NA
        }
        if(!utils::hasName(pp, "variable")) pp$variable <- v
        # (Re)set type
        pp$type <- ty
        # Attach
        out <- rbind(out, pp)
      }
    }
    assertthat::assert_that(nrow(out)>0,
                            any(out$presence.density>0),
                            all( unique(out$variable) %in% x.var )
                            )

    # --- #
    # Format the data for the plot into long format
    # Could be done easier with tidyr, but avoiding another dependency here
    plot.df <- data.frame()
    for(v in out$variable){
      o <- data.frame(variable = v,
                      layer = c(out$partial_effect[which(out$variable==v)],
                                out$partial_effect[which(out$variable==v)],
                                out$partial_effect[which(out$variable==v)]),
                      value = c(out$mean[which(out$variable==v)],
                                out$presence.density[which(out$variable==v)],
                                out$absence.density[which(out$variable==v)]),
                      source = c(rep("Suitability", 100),
                                 rep("Presence", 100),
                                 rep("Background/Absence", 100)))
      # Also add the actual variable density from the layers
      o$rug <- terra::spatSample(vars[[v]], nrow(o), na.rm= TRUE)[,1]
      plot.df <- rbind(plot.df, o)
    }

    if(df){ return(plot.df)}

    # Otherwise create the plot
    if(length(x.var) == 1){

      density.plot <- ggplot2::ggplot(data = plot.df,
                                       ggplot2::aes(x = layer, y = value)) +
        ggplot2::theme_bw() +
          # Add the lines
          ggplot2::geom_line(ggplot2::aes(colour = source, linetype = source)) +
          ggplot2::scale_color_manual(values = c("green4", "red", "blue")) +
          ggplot2::scale_linetype_manual(values = c( "dashed", "twodash", "solid")) +
        # Add a rug
        ggplot2::geom_rug(data = dplyr::sample_n(plot.df,100), ggplot2::aes(x = rug),sides = "b",
                          alpha = .7,length = ggplot2::unit(0.05, "npc")) +
        ggplot2::labs(x = x.var, y = "Value") +
        # ggplot2::facet_wrap(~variable,scales = "free_x") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::theme(legend.title = ggplot2::element_blank())

    } else {
      stop("NOT yet implemented! Need to adapt partial calculation to allow multiple values.")
      # Build new data.frame container
      new <- data.frame(
        var1 = seq(vars_lims[1,1],vars_lims[2,1], length.out = 100), # plot.df$rug[plot.df$variable==x.var[1]]
        var2 = seq(vars_lims[1,2],vars_lims[2,2], length.out = 100)
        )

      # Bivariate plot instead
      density.plot <- ggplot2::ggplot(data = new, ggplot2::aes(x = var1, y = var2)) +
        ggplot2::theme_bw() +
        # Add the raster
        ggplot2::geom_raster(ggplot2::aes(fill = pred)) +
        ggplot2::scale_fill_viridis_c(option = "B", guide = ggplot2::guide_colourbar(title = "Suitability")) +
        ggplot2::labs(x = x.var[1], y = x.var[2], title = "Predicted suitability in environment space") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::theme(legend.title = ggplot2::element_blank())
    }
    return(density.plot)
  }
)

