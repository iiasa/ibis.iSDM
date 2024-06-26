#' Function to create an ensemble of multiple fitted models
#'
#' @description Ensemble models calculated on multiple models have often been
#' shown to outcompete any single model in comparative assessments (Valavi et
#' al. 2022).
#'
#' This function creates an ensemble of multiple provided distribution models
#' fitted with the [`ibis.iSDM-package`]. Each model has to have estimated
#' predictions with a given method and optional uncertainty in form of the
#' standard deviation or similar. Through the `layer` parameter it can be
#' specified which part of the prediction should be averaged in an ensemble.
#' This can be for instance the *mean* prediction and/or the standard deviation
#' *sd*. See Details below for an overview of the different methods.
#'
#' Also returns a coefficient of variation (cv) as output of the ensemble, but
#' note this should not be interpreted as measure of model uncertainty as it
#' cannot capture parameter uncertainty of individual models; rather it reflects
#' variation among predictions which can be due to many factors including simply
#' differences in model complexity.
#'
#' @param ... Provided [`DistributionModel`] or [`SpatRaster`] objects.
#' @param method Approach on how the ensemble is to be created. See details for
#' available options (Default: \code{'mean'}).
#' @param weights (*Optional*) weights provided to the ensemble function if
#' weighted means are to be constructed (Default: \code{NULL}).
#' @param min.value A optional [`numeric`] stating a minimum value that needs
#' to be surpassed in each layer before calculating and ensemble (Default: \code{NULL}).
#' @param layer A [`character`] of the layer to be taken from each prediction
#' (Default: \code{'mean'}). If set to \code{NULL} ignore any of the layer
#' names in ensembles of `SpatRaster` objects.
#' @param normalize [`logical`] on whether the inputs of the ensemble should be
#' normalized to a scale of 0-1 (Default: \code{FALSE}).
#' @param uncertainty A [`character`] indicating how the uncertainty among
#' models should be calculated. Available options include \code{"none"}, the
#' standard deviation (\code{"sd"}), the average of all PCA axes except the
#' first \code{"pca"}, the coefficient of variation (\code{"cv"}, Default) or
#' the range between the lowest and highest value (\code{"range"}).
#' @param apply_threshold A [`logical`] flag (Default: \code{TRUE}) specifying
#' whether threshold values should also be created via \code{"method"}. Only
#' applies and works for [`DistributionModel`] and thresholds found.
#'
#' @details Possible options for creating an ensemble includes:
#' * \code{'mean'} - Calculates the mean of several predictions.
#' * \code{'median'} - Calculates the median of several predictions.
#' * \code{'max'} - The maximum value across predictions.
#' * \code{'min'} - The minimum value across predictions.
#' * \code{'weighted.mean'} - Calculates a weighted mean. Weights have to be supplied separately (e.g. TSS).
#' * \code{'min.sd'} - Ensemble created by minimizing the uncertainty among predictions.
#' * \code{'threshold.frequency'} - Returns an ensemble based on threshold frequency (simple count). Requires thresholds to be computed.
#' * \code{'pca'} - Calculates a PCA between predictions of each algorithm and then extract the first axis (the one explaining the most variation).
#'
#' In addition to the different ensemble methods, a minimal threshold
#' (\code{min.value}) can be set that needs to be surpassed for averaging. By
#' default this option is not used (Default: \code{NULL}).
#'
#' Note by default only the band in the \code{layer} parameter is composited. If
#' supported by the model other summary statistics from the posterior (e.g.
#' \code{'sd'}) can be specified.
#'
#' @note If a list is supplied, then it is assumed that each entry in the list
#' is a fitted [`DistributionModel`] object. Take care not to create an ensemble
#' of models constructed with different link functions, e.g. logistic vs [log].
#' In this case the \code{"normalize"} parameter has to be set.
#'
#' @returns A [`SpatRaster`] object containing the ensemble of the provided
#'   predictions specified by \code{method} and a coefficient of variation
#'   across all models.
#'
#' @references
#' * Valavi, R., Guillera‐Arroita, G., Lahoz‐Monfort, J. J., & Elith, J. (2022).
#' Predictive performance of presence‐only species distribution models: a benchmark
#' study with reproducible code. Ecological Monographs, 92(1), e01486.
#'
#' @keywords train
#'
#' @examples
#' # Method works for fitted models as well as as rasters
#' r1 <- terra::rast(nrows = 10, ncols = 10, res = 0.05, xmin = -1.5,
#'  xmax = 1.5, ymin = -1.5, ymax = 1.5, vals = rnorm(3600,mean = .5,sd = .1))
#' r2 <- terra::rast(nrows = 10, ncols = 10, res = 0.05, xmin = -1.5,
#'  xmax = 1.5, ymin = -1.5, ymax = 1.5, vals = rnorm(3600,mean = .5,sd = .5))
#' names(r1) <- names(r2) <- "mean"
#'
#' # Assumes previously computed predictions
#' ex <- ensemble(r1, r2, method = "mean")
#'
#' terra::plot(ex)
#'
#' @name ensemble
NULL

#' @rdname ensemble
#' @export
methods::setGeneric("ensemble",
                    signature = methods::signature("..."),
                    function(..., method = "mean", weights = NULL, min.value = NULL, layer = "mean",
                             normalize = FALSE, uncertainty = "cv", apply_threshold = TRUE) standardGeneric("ensemble"))

#' @rdname ensemble
methods::setMethod(
  "ensemble",
  methods::signature("ANY"),
  function(..., method = "mean", weights = NULL, min.value = NULL, layer = "mean",
           normalize = FALSE, uncertainty = "cv", apply_threshold = TRUE){
    if(length(list(...))>1) {
      mc <- list(...)
    } else {
      # Collate provided models
      if(!is.list(...)){
        mc <- list(...)
      } else mc <- c(...)
    }

    # Get all those that are DistributionModels
    mods <- mc[ sapply(mc, function(x) inherits(x, "DistributionModel") ) ]
    if(length(mods)==0) {
      # Check whether scenario objects were not provided instead
      mods1 <- mc[ sapply(mc, function(x) inherits(x, "BiodiversityScenario") ) ]
      mods2 <- mc[ sapply(mc, function(x) is.Raster(x) ) ]
      mods3 <- mc[  sapply(mc, function(x) inherits(x, "stars") ) ]
      assertthat::assert_that(length(mods1)>0 || length(mods2)>0 || length(mods3)>0,
                              msg = "Ensemble only works with DistributionModel or BiodiversityScenario objects! Alternativly supply SpatRaster or stars objects.")
      if(length(mods1)>0) mods <- mods1 else if(length(mods2)>0) mods <- mods2 else mods <- mods3
    }

    # Further checks
    assertthat::assert_that(
      is.character(method),
      is.null(min.value) || is.numeric(min.value),
      is.null(layer) || ( is.character(layer) && length(layer) == 1 ),
      is.null(weights) || is.vector(weights),
      is.logical(normalize),
      is.character(uncertainty),
      is.logical(apply_threshold)
    )

    # Check the method
    method <- match.arg(method, c('mean', 'weighted.mean', 'median', 'max', 'min',
                                  'threshold.frequency', 'min.sd', 'pca'), several.ok = FALSE)
    # Uncertainty calculation
    uncertainty <- match.arg(uncertainty, c('none','sd', 'cv', 'range', 'pca'), several.ok = FALSE)

    # For Distribution model ensembles
    if( all( sapply(mods, function(z) inherits(z, "DistributionModel")) ) ){
      assertthat::assert_that(length(mods)>=2, # Need at least 2 otherwise this does not make sense
                              msg = "No use calculating an ensemble on one object only..."
      )
      # Check that layers all have a prediction layer
      assertthat::assert_that(
        all( sapply(mods, function(x) !is.Waiver(x$get_data('prediction')) ) ),
        msg = "All distribution models need a fitted prediction object!"
      )
      # Check that layer is present in supplied mods
      assertthat::assert_that(
        all( sapply(mods, function(x) layer %in% names(x$get_data('prediction')) ) ),
        msg = paste("Layer", text_red(layer), "not found in supplied objects! Prediction missing?")
      )

      # Get prediction stacks from all mods
      ll_ras <- sapply(mods, function(x) x$get_data('prediction')[[layer]])
      # Ensure that the layers have the same resolution, otherwise align
      if(!terra::compareGeom(ll_ras[[1]], ll_ras[[2]], stopOnError = FALSE)){
        if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Ensemble]','red','Rasters need to be aligned. Check.')
        ll_ras[[2]] <- terra::resample(ll_ras[[2]], ll_ras[[1]], method = "bilinear")
      }

      # Now ensemble per layer entry
      out <- terra::rast()
      for(lyr in layer){
        ras <- terra::rast(sapply(ll_ras, function(x) x[[lyr]]))
        # Apply threshold if set. Set to 0 thus reducing the value of the
        # ensembled layer.
        if(!is.null(min.value)) ras[ras < min.value] <- 0

        # If normalize before running an ensemble if parameter set
        if(normalize) ras <- predictor_transform(ras, option = "norm")

        # Now create the ensemble depending on the option
        if(method == 'mean'){
          new <- terra::mean( ras, na.rm = TRUE)
        } else if(method == 'median'){
          new <- terra::median(ras, na.rm = TRUE)
        } else if(method == 'max'){
          new <- max(ras, na.rm = TRUE)
        } else if(method == 'min'){
          new <- min(ras, na.rm = TRUE)
        } else if(method == 'weighted.mean'){
          new <- terra::weighted.mean( ras, w = weights, na.rm = TRUE)
        } else if(method == 'threshold.frequency'){
          # Check that thresholds are available
          assertthat::assert_that(
            all( sapply(mods, function(x) length( grep("threshold", x$show_rasters()) )>0 ) ),
            msg = "This function requires thresholds to be computed!"
          )
          n_tr <- sapply(mods, function(x) grep("threshold", x$show_rasters(),value = TRUE) )
          # Get layer of each threshold if there are multiple
          ras_tr <- terra::rast()
          for(i in 1:length(n_tr)){
            o <- mods[[i]]$get_data(n_tr[i])
            # Grep layer name from the stack
            suppressWarnings(
              ras_tr <- c(ras_tr, o[[grep(layer, names(o))]] )
            )
          }
          # Calculate frequency
          new <- sum(ras_tr, na.rm = TRUE)
          new <- terra::mask(new, ras_tr[[1]])
        } else if(method == 'min.sd'){
          # If method 'min.sd' furthermore check that there is a sd object for
          # all of them
          assertthat::assert_that(
            all( sapply(mods, function(x) "sd" %in% names(x$get_data('prediction')) ) ),
            msg = "Method \'min.sd\' needs parametrized uncertainty (sd) for all objects."
          )
          # Also get SD prediction from models
          ras_sd <- terra::rast( sapply(mods, function(x) x$get_data('prediction')[['sd']]))
          # Normalize the sds for each
          ras_sd <- predictor_transform(ras_sd, option = "norm")
          # Get the id of the layer where standard deviation is lower
          min_sd <- terra::where.min(ras_sd)
          new <- emptyraster(ras)
          for(cl in terra::unique(min_sd)){
            new[min_sd == cl] <- ras[[cl]][min_sd == cl]
          }
        } else if(method == 'pca'){
          # Calculate a pca on the layers and return the first axes
          new <- predictor_transform(ras, option = "pca", pca.var = 1)[[1]]
        }

        # Rename
        names(new) <- paste0("ensemble_", lyr)
        # Add attributes on the method of ensembling
        attr(new, "method") <- method
        if(uncertainty!='none'){
          if(uncertainty == "pca") {
            # If PCA selected, calculate uncertainty based as average of all PCA
            # axes (except first)
            rasp <- predictor_transform(ras, option = "pca")
            rasp <- subset(rasp, 2:terra::nlyr(rasp))
          } else rasp <- NULL
          # Add uncertainty
          ras_uncertainty <- switch (uncertainty,
                                     "sd" = terra::app(ras, stats::sd, na.rm = TRUE),
                                     "cv" = terra::app(ras, stats::sd, na.rm = TRUE) / terra::mean(ras, na.rm = TRUE),
                                     "range" = max(ras, na.rm = TRUE) - min(ras, na.rm = TRUE),
                                     "pca" = terra::mean(rasp, na.rm = TRUE)
          )
          names(ras_uncertainty) <- paste0(uncertainty, "_", lyr)
          # Add attributes on the method of ensembling
          attr(ras_uncertainty, "method") <- uncertainty

          # Add all layers to out
          suppressWarnings( out <- c(out, new, ras_uncertainty) )
        } else {
          suppressWarnings( out <- c(out, new) )
        }
      }

      # Check for threshold values and collate
      if(apply_threshold){
        ll_val <- sapply(mods, function(x) x$get_thresholdvalue())
        # Incase no thresholds are found, ignore entirely
        if(!all(any(sapply(ll_val, is.Waiver)))){
          # Respecify weights as otherwise call below fails
          if(any(sapply(ll_val, is.Waiver))){
            if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Ensemble]','yellow','Threshold values not found for all objects')
            ll_val <- ll_val[-which(sapply(ll_val, is.Waiver))]
            ll_val <- ll_val |> as.numeric()
          }
          if(is.null(weights)) weights <- rep(1, length(ll_val))

          # Composite threshold
          tr <- dplyr::case_when(
            method == "mean" ~ mean(ll_val, na.rm = TRUE),
            method == "median" ~ median(ll_val, na.rm = TRUE),
            method == "max" ~ max(ll_val, na.rm = TRUE),
            method == "min" ~ min(ll_val, na.rm = TRUE),
            method == "weighted.mean" ~ weighted.mean(ll_val, w = weights, na.rm = TRUE),
            .default = mean(ll_val, na.rm = TRUE)
          )

          # Ensemble the first layer
          out <- c(out,
                   threshold(out[[1]], method = "fixed", value = tr)
          )
        }
      }
      assertthat::assert_that(is.Raster(out))

      return(out)
    } else if(is.Raster(mods[[1]])) {
      # Check that layer is present in supplied mods
      if(terra::nlyr(mods[[1]])>1 && length(mods) == 1){
        if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Ensemble]','red','Single multiband raster found. Ignoring parameter layer and ensemble.')
        layer <- 1
        ras <- mods[[1]]
      } else {
        if(!is.null(layer)){
          assertthat::assert_that(
            all( sapply(mods, function(x) layer %in% names(x) ) ),
            msg = paste("Layer", text_red(layer), "not found in supplied objects!")
          )
        } else { layer <- 1 } # Take the first one
        # Get prediction stacks from all mods
        ll_ras <- sapply(mods, function(x) x[[layer]])
        # Ensure that the layers have the same resolution, otherwise align
        if(!terra::compareGeom(ll_ras[[1]], ll_ras[[2]], stopOnError = FALSE)){
          if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Ensemble]','red','Rasters need to be aligned. Check.')
          ll_ras[[2]] <- terra::resample(ll_ras[[2]], ll_ras[[1]], method = "bilinear")
        }
        ras <- terra::rast(ll_ras)
      }
      # Apply threshold if set. Set to 0 thus reducing the value of the ensembled
      # layer.
      if(!is.null(min.value)) ras[ras < min.value] <- 0

      # If normalize before running an ensemble if parameter set
      if(normalize) ras <- predictor_transform(ras, option = "norm")

      # Check that weight lengths is equal to provided distribution objects
      if(!is.null(weights) && !is.numeric(layer)) assertthat::assert_that(length(weights) == length(mods))
      # If weights vector is numeric, standardize the weights
      if(is.numeric(weights)) {
        if(any(weights < 0)) weights[weights < 0] <- 0 # Assume those contribute anything
        weights <- weights / sum(weights)
      }

      # Now ensemble per layer entry
      out <- terra::rast()
      for(lyr in layer){

        # Now create the ensemble depending on the option
        if(method == 'mean'){
          new <- terra::mean( ras, na.rm = TRUE)
        } else if(method == 'median'){
          new <- terra::median( ras, na.rm = TRUE)
        } else if(method == 'max'){
          new <- max(ras, na.rm = TRUE)
        } else if(method == 'min'){
          new <- min(ras, na.rm = TRUE)
        } else if(method == 'weighted.mean'){
          new <- terra::weighted.mean( ras, w = weights, na.rm = TRUE)
        } else if(method == 'threshold.frequency'){
          # Check that thresholds are available
          stop("This function does not (yet) work with directly provided Raster objects.")

        } else if(method == 'min.sd'){
          # If method 'min.sd' furthermore check that there is a sd object for all
          # of them
          assertthat::assert_that(
            all( sapply(mods, function(x) "sd" %in% names(mods) ) ),
            msg = "Method \'min.sd\' needs parametrized uncertainty (sd) for all objects."
          )
          # Also get SD prediction from models
          ras_sd <- c( sapply(mods, function(x) x[['sd']]))
          # Get the id of the layer where standard deviation is lower
          min_sd <- terra::where.min(ras_sd)
          new <- emptyraster(ras)
          for(cl in terra::unique(min_sd)){
            new[min_sd == cl] <- ras[[cl]][min_sd == cl]
          }
        } else if(method == 'pca'){
          # Calculate a pca on the layers and return the first axes
          new <- predictor_transform(ras, option = "pca", pca.var = 1)[[1]]
        }
        # Rename
        names(new) <- paste0("ensemble_", lyr)
        # Add attributes on the method of ensemble
        attr(new, "method") <- method
        if(uncertainty != "none"){
          if(uncertainty == "pca") {
            stop("Currently, uncertainty = 'pca' is not implemented for SpatRaster input.")
          }
          # Add uncertainty
          ras_uncertainty <- switch (uncertainty,
                                     "sd" = terra::app(ras, fun = "sd", na.rm = TRUE),
                                     "cv" = terra::app(ras, fun = "sd", na.rm = TRUE) / terra::mean(ras, na.rm = TRUE),
                                     "range" = max(ras, na.rm = TRUE) - min(ras, na.rm = TRUE)
          )
          names(ras_uncertainty) <- paste0(uncertainty, "_", lyr)
          # Add attributes on the method of ensembling
          attr(ras_uncertainty, "method") <- uncertainty
          # Add all layers to out
          suppressWarnings( out <- c(out, new, ras_uncertainty) )
        } else {
          suppressWarnings( out <- c(out, new) )
        }
      }

      assertthat::assert_that(is.Raster(out))
      return(out)
    } else {
      # Scenario objects as stars or Scenario objects
      if(all(sapply(mods, function(z) inherits(z, "stars")))){
        # Check that layer is in stars
        if(!assertthat::see_if(all( sapply(mods, function(z) layer %in% names(z)) ))){
          if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Ensemble]','red','Provided layer not in objects. Taking first option!')
          layer <- names(mods[[1]])[1]
        }
        # Format to table
        lmat <- do.call("rbind", mods) |> as.data.frame()
        # Get dimensions
        lmat_dim <- stars::st_dimensions(mods[[1]])

      } else {
        # Check that layers all have a prediction layer
        assertthat::assert_that(
          all( sapply(mods, function(x) !is.Waiver(x$get_data()) ) ),
          msg = "All distribution models need a fitted scenario object!"
        )
        # Check that layer is present in supplied mods
        assertthat::assert_that(
          all( sapply(mods, function(x) layer %in% names(x$get_data()) ) ),
          msg = paste("Layer", text_red(layer), "not found in supplied objects!")
        )
        # Get projected suitability from all mods
        lmat <- stars::st_as_stars(
          sapply(mods, function(x) x$get_data()[layer])
        ) |> as.data.frame()
        # Get dimensions
        lmat_dim <- stars::st_dimensions(mods[[1]]$get_data())
      }

      # Normalize stars files
      if(normalize){
        # Get overall means and max values
        ovmin <- min(lmat[,4:ncol(lmat)],na.rm = TRUE)
        ovmax <- max(lmat[,4:ncol(lmat)],na.rm = TRUE)
        lmat[,4:ncol(lmat)] <- apply(lmat[,4:ncol(lmat)], # On the assumption that col 1-3 are coordinates+time
                                     2, function(x) {
                                       (x - ovmin) / (ovmax - ovmin )
                                     })
      }

      # Now create the ensemble depending on the option
      if(method == 'mean'){
        out <- apply(lmat[,4:ncol(lmat)], # On the assumption that col 1-3 are coordinates+time
                     1, function(x) mean(x, na.rm = TRUE))
      } else if(method == 'median'){
        out <- apply(lmat[,4:ncol(lmat)], # On the assumption that col 1-3 are coordinates+time
                     1, function(x) stats::median(x, na.rm = TRUE))
      } else if(method == 'max'){
        out <- apply(lmat[,4:ncol(lmat)], # On the assumption that col 1-3 are coordinates+time
                     1, function(x) max(x, na.rm = TRUE))
      } else if(method == 'min'){
        out <- apply(lmat[,4:ncol(lmat)], # On the assumption that col 1-3 are coordinates+time
                     1, function(x) min(x, na.rm = TRUE))
      } else if(method == 'weighted.mean'){
        out <- apply(lmat[,4:ncol(lmat)], # On the assumption that col 1-3 are coordinates+time
                     1, function(x) weighted.mean(x, w = weights, na.rm = TRUE))
      } else if(method == 'threshold.frequency'){
        out <- apply(lmat[,4:ncol(lmat)], # On the assumption that col 1-3 are coordinates+time
                     1, function(x) sum(x, na.rm = TRUE) / (ncol(lmat)-3) )
        # Check that thresholds are available
      } else if(method == 'min.sd'){
        stop("This has not been reasonably implemented in this context.")
      } else if(method == 'pca'){
        stop("This has not been reasonably implemented in this context.")
      }
      # Add dimensions to output
      out <- cbind( sf::st_coordinates(mods[[1]]$get_data()[layer]), "ensemble" = out ) |> as.data.frame()

      # Convert to stars
      out <- out |> stars:::st_as_stars.data.frame(dims = c(1,2,3), coords = 1:2)
      # Rename dimension names
      out <- out |> stars::st_set_dimensions(names = c("x", "y", "band"))
      # Rename
      names(out) <- paste0("ensemble_", layer)
      # Add attributes on the method of ensemble
      attr(out, "method") <- method

      # Check for threshold values and collate
      if(apply_threshold){
        ll_val <- sapply(mods, function(x) x$get_thresholdvalue())
        # Incase no thresholds are found, ignore entirely
        if(!all(any(sapply(ll_val, is.Waiver)))){
          # Respecify weights as otherwise call below fails
          if(any(sapply(ll_val, is.Waiver))){
            if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Ensemble]','yellow','Threshold values not found for all objects')
            ll_val <- ll_val[-which(sapply(ll_val, is.Waiver))]
            ll_val <- ll_val |> as.numeric()
          }
          if(is.null(weights)) weights <- rep(1, length(ll_val))

          # Composite threshold
          tr <- dplyr::case_when(
            method == "mean" ~ mean(ll_val, na.rm = TRUE),
            method == "median" ~ median(ll_val, na.rm = TRUE),
            method == "max" ~ max(ll_val, na.rm = TRUE),
            method == "min" ~ min(ll_val, na.rm = TRUE),
            method == "weighted.mean" ~ weighted.mean(ll_val, w = weights, na.rm = TRUE),
            .default = mean(ll_val, na.rm = TRUE)
          )

          # reclassify to binary
          new <- out
          new[new < tr[[1]]] <- 0; new[new >= tr[[1]]] <- 1
          names(new) <- 'ensemble_threshold'
          out <- c(out, new)
        }
      }

      # --- #
      if(uncertainty != 'none'){
        if(uncertainty == "pca") {
          stop("Currently, uncertainty = 'pca' is not implemented for stars input.")
        }
        # Add uncertainty
        out_uncertainty <- switch (uncertainty,
                                   "sd" = apply(lmat[,4:ncol(lmat)], 1, function(x) stats::sd(x, na.rm = TRUE)),
                                   "cv" = apply(lmat[,4:ncol(lmat)], 1, function(x) stats::sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)),
                                   "range" = apply(lmat[,4:ncol(lmat)], 1, function(x) (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
        )
        if(any(is.infinite(out_uncertainty))) out_uncertainty[is.infinite(out_uncertainty)] <- NA
        # Add dimensions to output
        out_uncertainty <- cbind( sf::st_coordinates(mods[[1]]$get_data()[layer]), "ensemble" = out_uncertainty ) |> as.data.frame()

        # Convert to stars
        out_uncertainty <- out_uncertainty |> stars:::st_as_stars.data.frame(dims = c(1,2,3), coords = 1:2)
        # Rename dimension names
        out_uncertainty <- out_uncertainty |> stars::st_set_dimensions(names = c("x", "y", "band"))
        # Rename
        names(out_uncertainty) <- paste0(uncertainty, "_", layer)
        # Add attributes on the method of ensembling
        attr(out_uncertainty, "method") <- uncertainty
        # --- #
        # Combine both ensemble and uncertainty
        ex <- stars:::c.stars(out, out_uncertainty)
        # Correct projection is unset
        if(is.na(sf::st_crs(ex))) ex <- sf::st_set_crs(ex, sf::st_crs(mods[[1]]$get_data()))
      } else {
        # Only the output
        ex <- out
      }
      # Correct projection is unset
      if(is.na(sf::st_crs(ex))) ex <- sf::st_set_crs(ex, sf::st_crs(mods[[1]]$get_data()))
      assertthat::assert_that(inherits(ex, "stars"))
      return(ex)
    }
  }
)

#### Ensemble partial ----

#' Function to create an ensemble of partial effects from multiple models
#'
#' @description Similar to the `ensemble()` function, this function creates an
#' ensemble of partial responses of provided distribution models fitted with
#' the [`ibis.iSDM-package`]. Through the `layer` parameter it can be
#' specified which part of the partial prediction should be averaged in an
#' ensemble (if given). This can be for instance the *mean* prediction and/or
#' the standard deviation *sd*. Ensemble partial is also being called if more
#' than one input [`DistributionModel`] object is provided to `partial`.
#'
#' By default the ensemble of partial responses is created as average across
#' all models with the uncertainty being the standard deviation of responses.
#'
#' @param ... Provided [`DistributionModel`] objects from which partial responses
#' can be called. In the future provided data.frames might be supported as well.
#' @param x.var A [`character`] of the variable from which an ensemble is to be
#' created.
#' @param method Approach on how the ensemble is to be created. See details for
#' options (Default: \code{'mean'}).
#' @param layer A [`character`] of the layer to be taken from each prediction
#' (Default: \code{'mean'}). If set to \code{NULL} ignore any of the layer names
#' in ensembles of `SpatRaster` objects.
#' @param newdata A optional [`data.frame`] or [`SpatRaster`] object supplied to
#' the model (DefaultL \code{NULL}). This object needs to have identical names as the original predictors.
#' @param normalize [`logical`] on whether the inputs of the ensemble should be
#' normalized to a scale of 0-1 (Default: \code{TRUE}).
#'
#' @details Possible options for creating an ensemble includes:
#' * \code{'mean'} - Calculates the mean of several predictions.
#' * \code{'median'} - Calculates the median of several predictions.
#'
#' @note If a list is supplied, then it is assumed that each entry in the list
#' is a fitted [`DistributionModel`] object. Take care not to create an ensemble
#' of models constructed with different link functions, e.g. logistic vs [log].
#' By default the response functions of each model are normalized.
#'
#' @returns A [data.frame] with the combined partial effects of the supplied models.
#'
#' @keywords train partial
#'
#' @examples
#' \dontrun{
#'  # Assumes previously computed models
#'  ex <- ensemble_partial(mod1, mod2, mod3, method = "mean")
#' }
#'
#' @name ensemble_partial
NULL

#' @rdname ensemble_partial
#' @export
methods::setGeneric("ensemble_partial",
                    signature = methods::signature("..."),
                    function(..., x.var, method = "mean", layer = "mean", newdata = NULL, normalize = TRUE) standardGeneric("ensemble_partial"))

#' @rdname ensemble_partial
methods::setMethod(
  "ensemble_partial",
  methods::signature("ANY"),
  function(..., x.var, method = "mean", layer = "mean", newdata = NULL, normalize = TRUE){
    assertthat::assert_that(
      is.character(x.var),
      msg = "Partial ensemble requires explicit specification of the parameter x.var."
    )
    if(length(list(...))>1) {
      mc <- list(...)
    } else {
      # Collate provided models
      if(!is.list(...)){
        mc <- list(...)
      } else mc <- c(...)
    }
    assertthat::assert_that(is.null(newdata) || is.data.frame(newdata),
                            msg = "Provide new data as data.frame.")

    # Get all those that are DistributionModels
    mods <- mc[ sapply(mc, function(x) inherits(x, "DistributionModel") ) ]

    if(length(mods)==1){
      # Only one object provided, just return partial results for it
      obj <- mods[[1]]
      return( obj$partial(x.var = x.var) )
    }

    # Further checks
    assertthat::assert_that(
      is.character(method),
      is.null(layer) || is.character(layer),
      is.logical(normalize)
    )

    # Check the method
    method <- match.arg(method, c('mean', 'median'), several.ok = FALSE)

    if(getOption("ibis.setupmessages", default = TRUE)) myLog("[Inference]","green","Creating a partial ensemble...")

    # Get variable range from the first object assuming they have similar variables
    # FIXME: Ideally make a consensus, otherwise assumes that same predictor been used
    rr <- range(mods[[1]]$model$predictors[,x.var], na.rm = TRUE)
    assertthat::assert_that(length(rr)==2, !anyNA(rr))
    rr <- seq(rr[1], rr[2], length.out = 100)

    # Now for each object get the partial values for the target variable
    out <- data.frame()
    for(obj in mods){
      if(length(grep(x.var, summary(obj)[[1]]))==0){
        message(paste("Variable", text_red(x.var), "not found in",class(obj)[1]," Skipping!"))
        next()
      }
      # Get partial with identical variable length
      if(is.null(newdata)){
        o <- try({partial(mod = obj, x.var = x.var, variable_length = 100, values = rr, plot = FALSE)},silent = TRUE)
      } else {
        o <- try({partial(mod = obj, x.var = x.var, newdata = newdata, plot = FALSE)},silent = TRUE)
      }
      if(inherits(o,"try-error")) next() # Skip, variable likely regularized out.
      # Subset to target variable
      o <- o[, c("partial_effect", layer)]
      # Normalize if set
      if(normalize){
        if(length(unique(o[[layer]]))>1){
          o[[layer]] <- (o[[layer]] - min( o[[layer]])) / (max(o[[layer]] ) - min(o[[layer]] ))
        } else {
          o[[layer]] <- 0 # Assumption being the variable has been regularized out
        }
      }
      o$cid <- 1:nrow(o)
      o$id <- as.character(obj$id)
      out <- rbind(out, o)
    }

    # Catch error in case none of them computed
    if(nrow(out)==0){
      if(getOption("ibis.setupmessages", default = TRUE)) myLog("[Inference]","red","None of the models seemed to contain the variable.")
      stop("No estimates found!")
    }

    # Now composite the ensemble depending on the option
    if(method == 'mean'){
      new <- aggregate(out[,layer], by = list(partial_effect = out$partial_effect),
                       FUN = function(x = out[[layer]]) {
                         return(cbind( mean = mean(x,na.rm = TRUE),
                                       sd = stats::sd(x,na.rm = TRUE)))
                       }) |> as.matrix() |> as.data.frame()
      colnames(new) <- c("partial_effect", "mean", "sd")
    } else if(method == 'median'){
      new <- aggregate(out[,layer], by = list(partial_effect = out$partial_effect),
                       FUN = function(x = out[[layer]]) {
                         return(cbind( median = stats::median(x,na.rm = TRUE),
                                       mad = stats::mad(x,na.rm = TRUE)))
                       }) |> as.matrix() |> as.data.frame()
      colnames(new) <- c("partial_effect", "median", "mad")
    }

    # FIXME:
    # A workaround fix specifically if a BART model is found
    # Reason is that BART partial effect calculation currently uses
    # its own step function, thus altering found relationships
    if(any(is.na(new[,3]))){
      # Rematch based on distance and aggregate
      new[,1] <- sapply(new[,1], function(z) rr[which.min(abs(z - rr))] )
      if(method=="mean"){
        new <- new |> dplyr::group_by(partial_effect) |>
          dplyr::summarise(sd = stats::sd(mean,na.rm=TRUE),
                           mean = mean(mean,na.rm=TRUE)
          ) |>
          dplyr::relocate(mean,.before = sd)
      } else {
        new <- new |> dplyr::group_by(partial_effect) |>
          dplyr::summarise(mad = stats::mad(median,na.rm=TRUE),
                           median = median(median,na.rm=TRUE)
          ) |>
          dplyr::relocate(median,.before = mad)
      }
    }
    return(new)
  }
)

#### Ensemble spartial ----

#' Function to create an ensemble of spartial effects from multiple models
#'
#' @inherit ensemble_partial description
#'
#' @inheritParams ensemble_partial
#' @param min.value A optional [`numeric`] stating a minimum value that needs
#'   to be surpassed in each layer before calculating and ensemble
#'   (Default: \code{NULL}).
#'
#' @inherit ensemble_partial details
#' @inherit ensemble_partial note
#'
#' @returns A [SpatRaster] object with the combined partial effects of the supplied models.
#'
#' @keywords train partial
#'
#' @examples
#' \dontrun{
#'  # Assumes previously computed models
#'  ex <- ensemble_spartial(mod1, mod2, mod3, method = "mean")
#' }
#'
#' @name ensemble_spartial
NULL

#' @rdname ensemble_spartial
#' @export
methods::setGeneric("ensemble_spartial",
                    signature = methods::signature("..."),
                    function(..., x.var, method = "mean", layer = "mean",
                             newdata = NULL, min.value = NULL, normalize = TRUE) standardGeneric("ensemble_spartial"))

#' @rdname ensemble_spartial
methods::setMethod(
  "ensemble_spartial",
  methods::signature("ANY"),
  function(..., x.var, method = "mean", layer = "mean", newdata = NULL,
           min.value = NULL, normalize = TRUE){
    assertthat::assert_that(
      is.character(x.var),
      msg = "Spartial ensemble requires explicit specification of the parameter x.var."
    )
    if(length(list(...))>1) {
      mc <- list(...)
    } else {
      # Collate provided models
      if(!is.list(...)){
        mc <- list(...)
      } else mc <- c(...)
    }
    assertthat::assert_that(is.null(newdata) || is.Raster(newdata),
                            msg = "Provide new data as SpatRaster.")

    # Get all those that are DistributionModels
    mods <- mc[ sapply(mc, function(x) inherits(x, "DistributionModel") ) ]

    # Further checks
    assertthat::assert_that(length(mods)>0,
                            msg = "No DistributionModel found among provided objects.")
    assertthat::assert_that(
      is.character(method),
      is.null(layer) || is.character(layer),
      is.null(min.value) || is.numeric(min.value),
      is.logical(normalize)
    )

    # Check the method
    method <- match.arg(method, c('mean', 'median', 'max', 'min'), several.ok = FALSE)

    if(getOption("ibis.setupmessages", default = TRUE)) myLog("[Inference]","green","Creating a spartial ensemble...")

    # If new data is provided
    if(!is.null(newdata)){
      # Set all variables other than the target variables to their mean
      # First check that variables are aligned to used predictor objects
      assertthat::assert_that(all(x.var %in% names(newdata)),
                              msg = "Variable not found in newdata!")
      template <- terra::rast(mods[[1]]$model$predictors[,c("x", "y")],
                              crs = terra::crs(mods[[1]]$model$background),type = "xyz") |>
        emptyraster()
      newdata <- alignRasters(newdata, template,cl = FALSE)
      assertthat::assert_that(is.Raster(template),
                              terra::compareGeom(newdata,template))

      # Calculate global means
      means <- terra::global(newdata, "mean", na.rm = TRUE)
      assertthat::assert_that(x.var %in% rownames(means))
      # Set all variables except x.var to the means
      nd <- newdata |> terra::as.data.frame(xy = TRUE, na.rm = FALSE)
      for(val in names(nd)){
        if(val %in% c(x.var,"x", "y")) next()
        nd[[val]][!is.na(nd[[val]])] <- means[rownames(means)==val,1]
      }
    }

    # Now for each object get the partial values for the target variable
    out <- list()
    for(obj in mods){
      if(length(grep(x.var, summary(obj)[[1]]))==0){
        message(paste("Variable", text_red(x.var), "not found in",class(obj)[1]," Skipping!"))
        next()
      }
      # Get spartial
      if(is.null(newdata)){
        o <- try({spartial(mod = obj, x.var = x.var, plot = FALSE)},silent = TRUE)
      } else {
        o <- try({obj$project(newdata = nd, layer = layer)},silent = TRUE)
      }
      if(inherits(o, "try-error")){
        if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Inference]','red',paste0('Spartial calculation failed for ',class(obj)[1]))
        next()
      }
      assertthat::assert_that(is.Raster(o))
      if(terra::nlyr(o)>1) o <- o[layer]
      # Rename
      names(o) <- layer
      # Append to output object
      out[[as.character(obj$id)]] <- o
    }
    # Catch error in case none of them computed
    if(length(out)==0){
      if(getOption("ibis.setupmessages", default = TRUE)) myLog("[Inference]","red","None of the models seemed to contain the variable.")
      stop("No estimates found!")
    }

    if(length(out)==1){
      if(getOption("ibis.setupmessages", default = TRUE)) myLog("[Inference]","yellow","Only a single model was estimated. Returning output.")
      new <- out[[1]]
    } else {
      # Now construct an ensemble by calling ensemble directly
      new <- ensemble(out, method = method,
                      normalize = normalize,
                      layer = layer,
                      min.value = min.value,
                      uncertainty = "none")
    }

    assertthat::assert_that(
      is.Raster(new),msg = "Something went wrong with the ensemble calculation!"
    )
    names(new) <- paste0("ensemble_spartial__",x.var)
    return(new)
  }
)
