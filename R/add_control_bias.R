#' Add a specified variable which should be controlled for somehow.
#'
#' @description
#' Sampling and other biases are pervasive drivers of the spatial location of
#' biodiversity datasets. While the integration of other, presumably less biased
#' data can be one way of controlling for sampling biases, another way is to control
#' directly for the bias in the model. Currently supported methods are:
#'
#' [*] \code{"partial"} - An approach described by Warton et al. (2013) to control the biases in a model, by
#' including a specified variable ("layer") in the model, but "partialling" it out during the projection phase.
#' Specifically the variable is set to a specified value ("bias_value"), which is by default the minimum value observed
#' across the background.
#' [*] \code{"offset"} - Dummy method that points to the [`add_offset_bias()`] functionality (see note). Makes use of
#' offsets to factor out a specified bias variable.
#'
#' @note
#' **Covariate transformations applied to other predictors need to **
#' Another option to consider biases particular in Poisson-point process models is to remove them
#' through an offset. Functionality to do so is available through the [`add_offset_bias()`] method. Setting the
#' method to \code{"offset"} will automatically point to this option.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param layer A [`sf`] or [`SpatRaster`] object with the range for the target feature. Specify a variable that is not
#' already added to \code{"x"} to avoid issues with duplications.
#' @param method A [`character`] vector describing the method used for bias control. Available
#' options are \code{"partial"} (Default) and \code{"offset"}.
#' @param bias_value A [`numeric`] with a value for \code{"layer"}. Specifying a [`numeric`] value here sets \code{layer}
#' to the target value during projection. By default the value is set to the minimum value found in the layer (Default: \code{NULL}).
#' @param add [`logical`] specifying whether a new offset is to be added. Setting
#' this parameter to \code{FALSE} replaces the current offsets with the new one (Default: \code{TRUE}).
#' @param ... Other parameters or arguments (currently not supported).
#' @references
#' * Warton, D.I., Renner, I.W. and Ramp, D., 2013. Model-based control of observer bias for the analysis of presence-only data in ecology. PloS one, 8(11), p.e79168.
#' * Merow, C., Allen, J.M., Aiello-Lammens, M., Silander, J.A., 2016. Improving niche and range estimates with Maxent and point process models by integrating spatially explicit information. Glob. Ecol. Biogeogr. 25, 1022–1036. https://doi.org/10.1111/geb.12453
#' * Botella, C., Joly, A., Bonnet, P., Munoz, F., & Monestiez, P. (2021). Jointly estimating spatial sampling effort and habitat suitability for multiple species from opportunistic presence‐only data. Methods in Ecology and Evolution, 12(5), 933-945.
#' @returns Adds bias control option to a [`distribution`] object.
#' @keywords bias, offset
#' @examples
#' \dontrun{
#'  x <- distribution(background) |>
#'    add_predictors(covariates) |>
#'    add_control_bias(biasvariable, bias_value = NULL)
#' }
#' @name add_control_bias
NULL

#' @name add_control_bias
#' @rdname add_control_bias
#' @exportMethod add_control_bias
#' @export
methods::setGeneric(
  "add_control_bias",
  signature = methods::signature("x", "layer"),
  function(x, layer, method = "partial", bias_value = NULL, add = TRUE) standardGeneric("add_control_bias"))

#' @name add_control_bias
#' @rdname add_control_bias
#' @usage \S4method{add_control_bias}{BiodiversityDistribution, SpatRaster}(x, layer)
methods::setMethod(
  "add_control_bias",
  methods::signature(x = "BiodiversityDistribution", layer = "SpatRaster"),
  function(x, layer, method = "partial", bias_value = NULL, add = TRUE) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.Raster(layer),
                            is.character(method),
                            is.null(bias_value) || is.numeric(bias_value),
                            is.logical(add)
    )
    # Match method
    method <- match.arg(method, c("partial", "offset"), several.ok = FALSE)

    # Check that background and range align, otherwise raise error
    if(is_comparable_raster(layer, x$background)) {
      warning('Supplied layer does not align with background! Aligning them now...')
      layer <- alignRasters(layer, x$background, method = 'bilinear', func = mean, cl = FALSE)
    }

    # Calculate a default bias value if not already set
    if(is.null(bias_value)) bias_value <- terra::global(layer, stat = "min", na.rm = TRUE)

    # Check for infinite values
    assertthat::assert_that(
      terra::nlyr(layer) == length(bias_value),
      all( is.finite( terra::global(layer, "range", na.rm = TRUE)) ),
      msg = "Infinite values found in the layer (maybe log of 0?)."
    )

    # Now precede depending on method
    if(method == "partial"){
      if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding bias controlled variable...')

      if(!is.Waiver(x$get_biascontrol())){
        if(getOption('ibis.setupmessages')) myLog('[Setup]','yellow','Overwriting existing bias variable...')
      }
      # Add to bias control
      x <- x$set_biascontrol(layer, method, bias_value)

    } else if(method == "offset") {
      x <- x |> add_offset_bias(layer = layer, add = add)
    }
    return(x)
  }
)
