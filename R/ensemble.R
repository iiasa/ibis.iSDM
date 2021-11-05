#' Function to create an ensemble of multiple fitted models
#'
#' @description This function creates an ensemble of multiple provided distribution models
#' fitted with the [`ibis.iSDM-package`]. Each model has to have a predictions with a mean layer and
#' optional uncertainty in form of the standard deviation or similar.
#' @details Possible options for creating an ensemble includes:
#' 'mean' - Calculates the mean of several predictions
#' 'weighted.mean' - Calculates a weighted mean. Weights have to be supplied separately (e.g. TSS)
#' 'min.sd' - Ensemble created by minimizing the uncertainty among predictions
#' 'threshold.frequency' - Returns an ensemble based on threshold frequency (simple count). Requires thresholds to be computed
#' @param ... Provided [`DistributionModel`] objects
#' @param method Approach on how the ensemble is to be created. See details for options (Default: 'mean')
#' @param weights (Optional) weights provided to the ensemble function if weighted means are to be constructed. (Default: NULL)
#' @param layername A [`character`] of the layername to be taken from each prediction
#' @returns A [`RasterStack`] containing the ensemble of the provided predictions.
#' Containing both the set 'method' and a coefficient of variation across predictions.
#' Critically this should not be interpreted as measure of uncertainty as it cannot
#' capture parameter uncertainty of individual model; Rather it reflects prediction variation.
#' @export

#' @name ensemble
#' @aliases ensemble
#' @keywords train
#' @exportMethod ensemble
#' @export
NULL
methods::setGeneric("ensemble",
                    signature = methods::signature("..."),
                    function(..., method = "mean", weights = NULL, layername = "mean") standardGeneric("ensemble"))

#' @name ensemble
#' @rdname ensemble
#' @usage \S4method{ensemble}{ANY}(...)
methods::setMethod(
  "ensemble",
  methods::signature("ANY"),
  function(..., method = "mean", weights = NULL, layername = "mean"){
    # Collate provided models
    mc <- list(...)
    # Get all those that are DistributionModels
    mods <- mc[ sapply(mc, function(x) inherits(x, "DistributionModel") ) ]

    # Further checks
    assertthat::assert_that(
      length(mods)>=2, # Need at least 2 otherwise this does not make sense
      is.character(method),
      is.character(layername),
      is.null(weights) || is.vector(weights)
    )

    # Check the method
    method <- match.arg(method, c('mean', 'weighted.mean', 'threshold.frequency', 'min.sd'), several.ok = FALSE)

    # Check that layers all have a prediction layer
    assertthat::assert_that(
      all( sapply(mods, function(x) !is.Waiver(x$get_data('prediction')) ) ),
      msg = "All distribution models need a fitted prediction object!"
    )
    # Check that weight lengths is equal to provided distribution objects
    if(!is.null(weights)) assertthat::assert_that(length(weights) == length(mods))
    # If weights vector is numeric, standardize the weights
    if(is.numeric(weights)) weights <- weights / sum(weights)

    # Get mean prediction from all mods
    ras <- raster::stack( sapply(mods, function(x) x$get_data('prediction')[[layername]]))

    # Now create the ensemble depending on the option
    if(method == 'mean'){
      out <- mean( ras, na.rm = TRUE)
      names(out) <- "ensemble"
    } else if(method == 'weighted.mean'){
      out <- weighted.mean( ras, w = weights, na.rm = TRUE)
      names(out) <- "ensemble"
    } else if(method == 'threshold.frequency'){
      # Check that thresholds are available
      assertthat::assert_that(
        all( sapply(mods, function(x) length( grep("threshold", x$show_rasters()) )>0 ) ),
        msg = "This function requires thresholds to be computed!"
      )
      n_tr <- sapply(mods, function(x) grep("threshold", x$show_rasters(),value = TRUE) )
      n_tr <- grep("mean", n_tr, value = TRUE) # Get the means
      ras_tr <- raster::stack()
      for(l in 1:length(n_tr)){ ras_tr <- raster::addLayer(ras_tr, mods[[l]]$get_data(n_tr[l])) }
      # Calculate frequency
      out <- sum(ras_tr, na.rm = TRUE)
      out <- raster::mask(out, ras_tr[[1]])
      names(out) <- "ensemble"
    } else if(method == 'min.sd'){
      # If method 'min.sd' furthermore check that there is a sd object for all of them
      assertthat::assert_that(
        all( sapply(mods, function(x) "sd" %in% names(x$get_data('prediction')) ) ),
        msg = "Method \'min.sd\' needs parametrized uncertainty (sd) for all objects."
      )
      # Also get SD prediction from models
      ras_sd <- raster::stack( sapply(mods, function(x) x$get_data('prediction')[['sd']]))
      # Get the id of the layer where standard deviation is lower
      min_sd <- raster::whiches.min(ras_sd)
      out <- emptyraster(ras)
      for(cl in raster::unique(min_sd)){
        out[min_sd == cl] <- ras[[cl]][min_sd == cl]
      }
      names(out) <- "ensemble"
    }
    # Add a coefficient of variation
    ras_cv <- raster::cv(ras, na.rm = TRUE); names(ras_cv) <- "cv"
    out <- raster::stack(out, ras_cv)

    # Also calculate coefficient of variation across predictions
    assertthat::assert_that(is.Raster(out))
    return(out)
  }
)
