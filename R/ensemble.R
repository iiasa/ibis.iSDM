#' Function to create an ensemble of multiple fitted models
#'
#' @description This function creates an ensemble of multiple provided distribution models
#' fitted with the [`ibis.iSDM-package`]. Each model has to have estimated predictions with a given method and
#' optional uncertainty in form of the standard deviation or similar.
#'
#' Also returns a coefficient of variation (cv) as output of the ensemble, but note
#' this should not be interpreted as measure of model uncertainty as it cannot
#' capture parameter uncertainty of individual model; rather it reflects prediction variation.
#' @details
#' Possible options for creating an ensemble includes:
#' * \code{'mean'} - Calculates the mean of several predictions
#' * \code{'weighted.mean'} - Calculates a weighted mean. Weights have to be supplied separately (e.g. TSS)
#' * \code{'min.sd'} - Ensemble created by minimizing the uncertainty among predictions
#' * \code{'threshold.frequency'} - Returns an ensemble based on threshold frequency (simple count). Requires thresholds to be computed
#' @note
#' Take care not to create an ensemble of models constructed with different link functions, e.g. [logistic] vs [log]
#' @param ... Provided [`DistributionModel`] objects
#' @param method Approach on how the ensemble is to be created. See details for options (Default: \code{'mean'}).
#' @param weights (*Optional*) weights provided to the ensemble function if weighted means are to be constructed (Default: \code{NULL}).
#' @param layer A [`character`] of the layer to be taken from each prediction (Default: \code{'mean'}).
#' @returns A [`RasterStack`] containing the ensemble of the provided predictions specified by \code{method} and a
#' coefficient of variation across all models.

#' @name ensemble
#' @aliases ensemble
#' @keywords train
#' @exportMethod ensemble
#' @export
NULL
methods::setGeneric("ensemble",
                    signature = methods::signature("..."),
                    function(..., method = "mean", weights = NULL, layer = "mean") standardGeneric("ensemble"))

#' @name ensemble
#' @rdname ensemble
#' @usage \S4method{ensemble}{ANY}(...)
methods::setMethod(
  "ensemble",
  methods::signature("ANY"),
  function(..., method = "mean", weights = NULL, layer = "mean"){
    # Collate provided models
    mc <- list(...)
    # Get all those that are DistributionModels
    mods <- mc[ sapply(mc, function(x) inherits(x, "DistributionModel") ) ]
    if(length(mods)==0) {
      # Check whether scenario objects were not provided instead
      mods <- mc[ sapply(mc, function(x) inherits(x, "BiodiversityScenario") ) ]
      assertthat::assert_that(length(mods)>0,msg = "Ensemble only works with DistributionModel or BiodiversityScenario objects!")
    }

    # Further checks
    assertthat::assert_that(
      length(mods)>=2, # Need at least 2 otherwise this does not make sense
      is.character(method),
      is.character(layer),
      is.null(weights) || is.vector(weights)
    )

    # Check the method
    method <- match.arg(method, c('mean', 'weighted.mean', 'threshold.frequency', 'min.sd'), several.ok = FALSE)

    # For Distribution model ensembles
    if( inherits(mods[[1]], "DistributionModel") ){
      # Check that layers all have a prediction layer
      assertthat::assert_that(
        all( sapply(mods, function(x) !is.Waiver(x$get_data('prediction')) ) ),
        msg = "All distribution models need a fitted prediction object!"
      )
      # Check that weight lengths is equal to provided distribution objects
      if(!is.null(weights)) assertthat::assert_that(length(weights) == length(mods))
      # If weights vector is numeric, standardize the weights
      if(is.numeric(weights)) weights <- weights / sum(weights)

      # Get prediction stacks from all mods
      ll_ras <- sapply(mods, function(x) x$get_data('prediction')[[layer]])
      # Ensure that the layers have the same resolution, otherwise align
      if(!compareRaster(ll_ras[[1]], ll_ras[[2]], stopiffalse = FALSE)){
        if(getOption('ibis.setupmessages')) myLog('[Ensemble]','red','Rasters need to be aligned. Check.')
        ll_ras[[2]] <- raster::resample(ll_ras[[2]], ll_ras[[1]])
      }
      # Now ensemble per layer entry
      out <- raster::stack()
      for(lyr in layer){
        ras <- stack(sapply(ll_ras, function(x) x[[lyr]]))

        # Now create the ensemble depending on the option
        if(method == 'mean'){
          new <- mean( ras, na.rm = TRUE)
          names(new) <- paste0("ensemble_", lyr)
        } else if(method == 'weighted.mean'){
          new <- weighted.mean( ras, w = weights, na.rm = TRUE)
          names(new) <- paste0("ensemble_", lyr)
        } else if(method == 'threshold.frequency'){
          # Check that thresholds are available
          assertthat::assert_that(
            all( sapply(mods, function(x) length( grep("threshold", x$show_rasters()) )>0 ) ),
            msg = "This function requires thresholds to be computed!"
          )
          n_tr <- sapply(mods, function(x) grep("threshold", x$show_rasters(),value = TRUE) )
          # Get layer of each threshold if there are multiple
          ras_tr <- raster::stack()
          for(i in 1:length(n_tr)){
            o <- mods[[i]]$get_data(n_tr[i])
            # Grep layer name from the stack
            ras_tr <- raster::addLayer(ras_tr, o[[grep(layer, names(o))]] )
          }
          # Calculate frequency
          new <- sum(ras_tr, na.rm = TRUE)
          new <- raster::mask(new, ras_tr[[1]])
          names(new) <- paste0("ensemble_", lyr)
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
          new <- emptyraster(ras)
          for(cl in raster::unique(min_sd)){
            new[min_sd == cl] <- ras[[cl]][min_sd == cl]
          }
          names(new) <- paste0("ensemble_", lyr)
        }
        # Add a coefficient of variation
        ras_cv <- raster::cv(ras, na.rm = TRUE); names(ras_cv) <- paste0("cv_", lyr)

        # Add all layers to out
        out <- raster::stack(out, new, ras_cv)
      }

      assertthat::assert_that(is.Raster(out))
      return(out)
  } else {
    # Scenario objects
    # Check that layers all have a prediction layer
    assertthat::assert_that(
      all( sapply(mods, function(x) !is.Waiver(x$get_scenarios()) ) ),
      msg = "All distribution models need a fitted scenario object!"
    )
    # Check that weight lengths is equal to provided distribution objects
    if(!is.null(weights)) assertthat::assert_that(length(weights) == length(mods))
    # If weights vector is numeric, standardize the weights
    if(is.numeric(weights)) weights <- weights / sum(weights)

    # Get projected suitability from all mods
    lmat <- stars::st_as_stars(
      sapply(mods, function(x) x$get_scenarios()[layer])
    ) |> as.data.frame()
    # Get dimensions
    lmat_dim <- stars:::st_dimensions(mods[[1]]$get_scenarios())

    # Now create the ensemble depending on the option
    if(method == 'mean'){
      out <- apply(lmat[,4:ncol(lmat)], # On the assumption that col 1-3 are coordinates+time
                   1, function(x) mean(x, na.rm = TRUE))
    } else if(method == 'weighted.mean'){
      out <- apply(lmat[,4:ncol(lmat)], # On the assumption that col 1-3 are coordinates+time
                   1, function(x) weighted.mean(x, w = weights, na.rm = TRUE))
    } else if(method == 'threshold.frequency'){
      stop("This has not been reasonably implemented in this context.")
      # Check that thresholds are available
    } else if(method == 'min.sd'){
      stop("This has not been reasonably implemented in this context.")
    }
    out <- cbind( lmat[,1:3], "ensemble" = out ) |> as.data.frame()

    # Convert to stars
    out <- out |> stars:::st_as_stars.data.frame(out,dims = c(1,2,3),coords = 1:2)
    # Rename dimension names
    out <- out |> stars:::st_set_dimensions(names = c("x", "y", "band"))
    # Also calculate coefficient of variation across predictions
    assertthat::assert_that(inherits(out, "stars"))
    return(out)
    }
  }
)
