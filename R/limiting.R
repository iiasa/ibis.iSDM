#' @include utils.R
NULL

#' Identify local limiting factor
#'
#' @description Calculates a [`SpatRaster`] of locally limiting factors from a
#' given projected model. To calculate this first the [`spartial`] effect of
#' each individual covariate in the model is calculated.
#'
#' The effect is estimated as that variable most responsible for decreasing
#' suitability at that cell. The decrease in suitability is calculated, for each
#' predictor in turn, relative to thesuitability that would be achieved if that
#' predictor took the value equal to the mean The predictor associated with the
#' largest decrease in suitability is the most limiting factor.
#'
#' @param mod A fitted \code{'DistributionModel'} object from which limited
#' factors are to be identified.
#' @param plot Should the result be plotted? (Default: \code{TRUE}).
#'
#' @return A `terra` object of the most important variable for a given grid cell.
#'
#' @concept Partly inspired by the rmaxent package.
#'
#' @references
#' * Elith, J., Kearney, M. and Phillips, S. (2010), The art of modelling range-shifting species.
#' Methods in Ecology and Evolution, 1: 330-342. doi: 10.1111/j.2041-210X.2010.00036.x
#'
#' @keywords partial
#'
#' @examples
#' \dontrun{
#' o <- limiting(fit)
#' plot(o)
#' }
#'
#' @name limiting
NULL

#' @rdname limiting
#' @export
methods::setGeneric("limiting",
                    signature = methods::signature("mod"),
                    function(mod, plot = TRUE) standardGeneric("limiting"))

#' @rdname limiting
methods::setMethod(
  "limiting",
  methods::signature("ANY"),
  function(mod, plot = TRUE){
    assertthat::assert_that(!missing(mod) || inherits(mod,'DistributionModel'),
                            is.logical(plot))

    # Assert that prediction is there
    assertthat::assert_that(is.Raster(mod$get_data()))

    # Get model object and variables respectively
    model <- mod$model
    pred <- mod$get_data()[[1]]
    # If threshold is set, mask
    if(!is.Waiver(mod$get_thresholdvalue())){
      tr <- mod$get_thresholdvalue()
      pred[pred < tr] <- NA
    }
    vars <- model$predictors_names

    # Output spatRaster
    ras <- terra::rast()
    pb <- progress::progress_bar$new(total = length(vars),
                                     format = "Processing :what")
    for(v in vars){
      pb$tick(tokens = list(what = v))

      # Build spartial model
      o <- try({spartial(mod, x.var = v, plot = FALSE)},silent = TRUE)
      # If error it is assumed that variable is not in model (regularized out)
      if(inherits(o, "try-error")) next()
      o <- o[[1]]# Get only the first layer (mean)
      names(o) <- v
      suppressWarnings( ras <- c(ras,o) )
      rm(o)
    }
    assertthat::assert_that(terra::nlyr(ras)>1,
                            msg = "Only a single or no spartial coefficients could be computed!")

    # Determine the maximum value
    out <- terra::which.max(ras - pred)
    names(out) <- "Limiting.variable"
    out <- terra::mask(out, model$background)

    # Define factor levels
    cls <- data.frame(id = unique(out)[,1])
    cls$cover <- names(ras)[cls$id]
    levels(out) <- cls
    assertthat::assert_that(is.factor(out))

    # Plot the result
    if(plot){
      cols <- ibis_colours$distinct_random[1:nrow(cls)]
      terra::plot(out, col = cols,
                  grid = FALSE,
                  smooth = FALSE,
                  all_levels = FALSE,
                  axes = FALSE)
    }
    return(out)
  }
)
