#' @include utils-spatial.R
NULL

#' Function to create an ensemble of multiple fitted models
#'
#' @description
#' Ensemble models calculated on multiple models have often been shown to
#' outcompete any single model in comparative assessments (Valavi et al. 2022).
#'
#' This function creates an ensemble of multiple provided distribution models
#' fitted with the [`ibis.iSDM-package`]. Each model has to have estimated predictions with a given method and
#' optional uncertainty in form of the standard deviation or similar.
#' Through the `layer` parameter it can be specified which part of the prediction
#' should be averaged in an ensemble. This can be for instance the *mean* prediction and/or
#' the standard deviation *sd*. See Details below for an overview of the different methods.
#'
#' Also returns a coefficient of variation (cv) as output of the ensemble, but note
#' this should not be interpreted as measure of model uncertainty as it cannot
#' capture parameter uncertainty of individual models; rather it reflects variation among predictions which
#' can be due to many factors including simply differences in model complexity.
#'
#' @details
#' Possible options for creating an ensemble includes:
#' * \code{'mean'} - Calculates the mean of several predictions.
#' * \code{'median'} - Calculates the median of several predictions.
#' * \code{'weighted.mean'} - Calculates a weighted mean. Weights have to be supplied separately (e.g. TSS).
#' * \code{'min.sd'} - Ensemble created by minimizing the uncertainty among predictions.
#' * \code{'threshold.frequency'} - Returns an ensemble based on threshold frequency (simple count). Requires thresholds to be computed.
#' * \code{'pca'} - Calculates a PCA between predictions of each algorithm and then extract the first axis (the one explaining the most variation).
#'
#' In addition to the different ensemble methods, a minimal threshold (\code{min.value}) can be set that needs to be surpassed for averaging.
#' By default this option is not used (Default: \code{NULL}).
#'
#' Note by default only the band in the \code{layer} parameter is composited. If supported by the model
#' other summary statistics from the posterior (e.g. \code{'sd'}) can be specified.
#'
#' @note
#' If a list is supplied, then it is assumed that each entry in the list is a fitted [`DistributionModel`] object.
#' Take care not to create an ensemble of models constructed with different link functions, e.g. [logistic] vs [log]. In this case
#' the \code{"normalize"} parameter has to be set.
#' @param ... Provided [`DistributionModel`] objects.
#' @param method Approach on how the ensemble is to be created. See details for available options (Default: \code{'mean'}).
#' @param weights (*Optional*) weights provided to the ensemble function if weighted means are to be constructed (Default: \code{NULL}).
#' @param min.value A [`numeric`] stating a minimum threshold value that needs to be surpassed in each layer (Default: \code{NULL}).
#' @param layer A [`character`] of the layer to be taken from each prediction (Default: \code{'mean'}). If set to \code{NULL}
#' ignore any of the layer names in ensembles of `Raster` objects.
#' @param normalize [`logical`] on whether the inputs of the ensemble should be normalized to a scale of 0-1 (Default: \code{FALSE}).
#' @param uncertainty A [`character`] indicating how the uncertainty among models should be calculated. Available options include
#' \code{"none"}, the standard deviation (\code{"sd"}), the coefficient of variation (\code{"cv"}, Default)
#' or the range between the lowest and highest value (\code{"range"}).
#' @references
#' * Valavi, R., Guillera‐Arroita, G., Lahoz‐Monfort, J. J., & Elith, J. (2022). Predictive performance of presence‐only species distribution models: a benchmark study with reproducible code. Ecological Monographs, 92(1), e01486.
#' @examples
#' \dontrun{
#'  # Assumes previously computed predictions
#'  ex <- ensemble(mod1, mod2, mod3, method = "mean")
#'  names(ex)
#'
#'  # Make a bivariate plot (might require other packages)
#'  bivplot(ex)
#' }
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
                    function(..., method = "mean", weights = NULL, min.value = NULL, layer = "mean",
                             normalize = FALSE, uncertainty = "cv") standardGeneric("ensemble"))

#' @name ensemble
#' @rdname ensemble
#' @usage \S4method{ensemble}{ANY}(...)
methods::setMethod(
  "ensemble",
  methods::signature("ANY"),
  function(..., method = "mean", weights = NULL, min.value = NULL, layer = "mean",
           normalize = FALSE, uncertainty = "cv"){
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
                              msg = "Ensemble only works with DistributionModel or BiodiversityScenario objects! Alternativel supply raster or stars objects.")
      if(length(mods1)>0) mods <- mods1 else if(length(mods2)>0) mods <- mods2 else mods <- mods3
    }

    # Further checks
    assertthat::assert_that(length(mods)>=2, # Need at least 2 otherwise this does not make sense
                            msg = "No use calculating an ensemble on one object only..."
    )
    assertthat::assert_that(
      is.character(method),
      is.null(min.value) || is.numeric(min.value),
      is.null(layer) || is.character(layer),
      is.null(weights) || is.vector(weights),
      is.logical(normalize),
      is.character(uncertainty)
    )

    # Check the method
    method <- match.arg(method, c('mean', 'weighted.mean', 'median', 'threshold.frequency', 'min.sd', 'pca'), several.ok = FALSE)
    # Uncertainty calculation
    uncertainty <- match.arg(uncertainty, c('none','sd', 'cv', 'range'), several.ok = FALSE)

    # Check that weight lengths is equal to provided distribution objects
    if(!is.null(weights)) assertthat::assert_that(length(weights) == length(mods))
    # If weights vector is numeric, standardize the weights
    if(is.numeric(weights)) weights <- weights / sum(weights)

    # For Distribution model ensembles
    if( all( sapply(mods, function(z) inherits(z, "DistributionModel")) ) ){
      # Check that layers all have a prediction layer
      assertthat::assert_that(
        all( sapply(mods, function(x) !is.Waiver(x$get_data('prediction')) ) ),
        msg = "All distribution models need a fitted prediction object!"
      )
      # Check that layer is present in supplied mods
      assertthat::assert_that(
        all( sapply(mods, function(x) layer %in% names(x$get_data('prediction')) ) ),
        msg = paste("Layer", text_red(layer), "not found in supplied objects!")
      )

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
        ras <- raster::stack(sapply(ll_ras, function(x) x[[lyr]]))

        # If normalize before running an ensemble if parameter set
        if(normalize) ras <- predictor_transform(ras, option = "norm")

        # Apply threshold if set. Set to 0 thus reducing the value of the ensembled layer.
        if(!is.null(min.value)) ras[ras < min.value] <- 0

        # Now create the ensemble depending on the option
        if(method == 'mean'){
          new <- mean( ras, na.rm = TRUE)
        } else if(method == 'median'){
          new <- raster::calc(ras, fun = median, na.rm = TRUE)
        } else if(method == 'weighted.mean'){
          new <- weighted.mean( ras, w = weights, na.rm = TRUE)
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
        } else if(method == 'min.sd'){
          # If method 'min.sd' furthermore check that there is a sd object for all of them
          assertthat::assert_that(
            all( sapply(mods, function(x) "sd" %in% names(x$get_data('prediction')) ) ),
            msg = "Method \'min.sd\' needs parametrized uncertainty (sd) for all objects."
          )
          # Also get SD prediction from models
          ras_sd <- raster::stack( sapply(mods, function(x) x$get_data('prediction')[['sd']]))
          # Normalize the sds for each
          ras_sd <- predictor_transform(ras_sd, option = "norm")
          # Get the id of the layer where standard deviation is lower
          min_sd <- raster::whiches.min(ras_sd)
          new <- emptyraster(ras)
          for(cl in raster::unique(min_sd)){
            new[min_sd == cl] <- ras[[cl]][min_sd == cl]
          }
        } else if(method == 'pca'){
          # Calculate a pca on the layers and return the first axes
          new <- predictor_transform(ras, option = "pca",pca.var = 1)[[1]]
        }

        # Rename
        names(new) <- paste0("ensemble_", lyr)
        # Add attributes on the method of ensembling
        attr(new, "method") <- method
        if(uncertainty!='none'){
          # Add uncertainty
          ras_uncertainty <- switch (uncertainty,
                                     "sd" = raster::calc(ras, sd, na.rm = TRUE),
                                     "cv" = raster::cv(ras, na.rm = TRUE),
                                     "range" = max(ras, na.rm = TRUE) - min(ras, na.rm = TRUE)
          )
          names(ras_uncertainty) <- paste0(uncertainty, "_", lyr)
          # Add attributes on the method of ensembling
          attr(ras_uncertainty, "method") <- uncertainty

          # Add all layers to out
          out <- raster::stack(out, new, ras_uncertainty)
        } else {
          out <- raster::stack(out, new)
        }
      }

      assertthat::assert_that(is.Raster(out))

      return(out)
  } else if(is.Raster(mods[[1]])) {
    # Check that layer is present in supplied mods
    if(!is.null(layer)){
      assertthat::assert_that(
        all( sapply(mods, function(x) layer %in% names(x) ) ),
        msg = paste("Layer", text_red(layer), "not found in supplied objects!")
      )
    } else { layer <- 1 } # Take the first one
    # TODO:
    if(length(layer)>1) stop("Not implemented yet")
    # Get prediction stacks from all mods
    ll_ras <- sapply(mods, function(x) x[[layer]])
    # Ensure that the layers have the same resolution, otherwise align
    if(!compareRaster(ll_ras[[1]], ll_ras[[2]], stopiffalse = FALSE)){
      if(getOption('ibis.setupmessages')) myLog('[Ensemble]','red','Rasters need to be aligned. Check.')
      ll_ras[[2]] <- raster::resample(ll_ras[[2]], ll_ras[[1]])
    }
    ras <- raster::stack(ll_ras)
    # If normalize before running an ensemble if parameter set
    if(normalize) ras <- predictor_transform(ras, option = "norm")

    # Apply threshold if set. Set to 0 thus reducing the value of the ensembled layer.
    if(!is.null(min.value)) ras[ras < min.value] <- 0

    # Now ensemble per layer entry
    out <- raster::stack()
    for(lyr in layer){

      # Now create the ensemble depending on the option
      if(method == 'mean'){
        new <- mean( ras, na.rm = TRUE)
      } else if(method == 'median'){
        new <- median( ras, na.rm = TRUE)
      } else if(method == 'weighted.mean'){
        new <- weighted.mean( ras, w = weights, na.rm = TRUE)
      } else if(method == 'threshold.frequency'){
        # Check that thresholds are available
        stop("This function does not (yet) work with directly provided Raster objects.")

      } else if(method == 'min.sd'){
        # If method 'min.sd' furthermore check that there is a sd object for all of them
        assertthat::assert_that(
          all( sapply(mods, function(x) "sd" %in% names(mods) ) ),
          msg = "Method \'min.sd\' needs parametrized uncertainty (sd) for all objects."
        )
        # Also get SD prediction from models
        ras_sd <- raster::stack( sapply(mods, function(x) x[['sd']]))
        # Get the id of the layer where standard deviation is lower
        min_sd <- raster::whiches.min(ras_sd)
        new <- emptyraster(ras)
        for(cl in raster::unique(min_sd)){
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
        # Add uncertainty
        ras_uncertainty <- switch (uncertainty,
                                   "sd" = raster::calc(ras, sd, na.rm = TRUE),
                                   "cv" = raster::cv(ras, na.rm = TRUE),
                                   "range" = max(ras, na.rm = TRUE) - min(ras, na.rm = TRUE)
        )
        names(ras_uncertainty) <- paste0(uncertainty, "_", lyr)
        # Add attributes on the method of ensembling
        attr(ras_uncertainty, "method") <- uncertainty
        # Add all layers to out
        out <- raster::stack(out, new, ras_uncertainty)
      } else {
        out <- raster::addLayer(out, new)
      }
    }

    assertthat::assert_that(is.Raster(out))
    return(out)
  } else {
    # Scenario objects as stars or Scenario objects
    if(all(sapply(mods, function(z) inherits(z, "stars")))){
      # Check that layer is in stars
      if(!assertthat::see_if(all( sapply(mods, function(z) layer %in% names(z)) ))){
        if(getOption('ibis.setupmessages')) myLog('[Ensemble]','red','Provided layer not in objects. Taking first option!')
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
    if(normalize){
      lmat[,4:ncol(lmat)] <- apply(lmat[,4:ncol(lmat)], # On the assumption that col 1-3 are coordinates+time
                                   2, function(x) {
                                     (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE) )
                                   })
    }

    # Now create the ensemble depending on the option
    if(method == 'mean'){
      out <- apply(lmat[,4:ncol(lmat)], # On the assumption that col 1-3 are coordinates+time
                   1, function(x) mean(x, na.rm = TRUE))
    } else if(method == 'median'){
      out <- apply(lmat[,4:ncol(lmat)], # On the assumption that col 1-3 are coordinates+time
                   1, function(x) stats::median(x, na.rm = TRUE))
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

    # --- #
    if(uncertainty != 'none'){
      # Add uncertainty
      out_uncertainty <- switch (uncertainty,
                                 "sd" = apply(lmat[,4:ncol(lmat)], 1, function(x) sd(x, na.rm = TRUE)),
                                 "cv" = apply(lmat[,4:ncol(lmat)], 1, function(x) raster::cv(x, na.rm = TRUE)),
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
      if(is.na(sf::st_crs(ex))) ex <- st_set_crs(ex, st_crs(mods[[1]]$get_data()))
    } else {
      # Only the output
      ex <- out
    }
    # Correct projection is unset
    if(is.na(sf::st_crs(ex))) ex <- st_set_crs(ex, st_crs(mods[[1]]$get_data()))
    assertthat::assert_that(inherits(ex, "stars"))
    return(ex)
    }
  }
)

#' Function to create an ensemble of partial effects from multiple models
#'
#' @description Similar to the `ensemble()` function, this function creates an ensemble of
#' partial responses of provided distribution models fitted with the [`ibis.iSDM-package`].
#' Through the `layer` parameter it can be specified which part of the partial prediction
#' should be averaged in an ensemble (if given). This can be for instance the *mean* prediction and/or
#' the standard deviation *sd*. Ensemble partial is also being called if more than one input
#' [`DistributionModel`] object is provided to `partial`.
#'
#' By default the ensemble of partial responses is created as average across all models with the
#' uncertainty being the standard deviation of responses.
#'
#' @details
#' Possible options for creating an ensemble includes:
#' * \code{'mean'} - Calculates the mean of several predictions.
#' * \code{'median'} - Calculates the median of several predictions.
#'
#' @note
#' If a list is supplied, then it is assumed that each entry in the list is a fitted [`DistributionModel`] object.
#' Take care not to create an ensemble of models constructed with different link functions, e.g. [logistic] vs [log].
#' By default the response functions of each model are normalized.
#' @param ... Provided [`DistributionModel`] objects from which partial responses can be called. In the future provided data.frames might be supported as well.
#' @param x.var A [`character`] of the variable from which an ensemble is to be created.
#' @param method Approach on how the ensemble is to be created. See details for options (Default: \code{'mean'}).
#' @param layer A [`character`] of the layer to be taken from each prediction (Default: \code{'mean'}). If set to \code{NULL}
#' ignore any of the layer names in ensembles of `Raster` objects.
#' @param normalize [`logical`] on whether the inputs of the ensemble should be normalized to a scale of 0-1 (Default: \code{TRUE}).
#' @returns A [`RasterStack`] containing the ensemble of the provided predictions specified by \code{method} and a
#' coefficient of variation across all models.

#' @name ensemble_partial
#' @aliases ensemble_partial
#' @keywords train
#' @exportMethod ensemble_partial
#' @export
NULL
methods::setGeneric("ensemble_partial",
                    signature = methods::signature("..."),
                    function(..., x.var, method = "mean", layer = "mean", normalize = TRUE) standardGeneric("ensemble_partial"))

#' @name ensemble_partial
#' @rdname ensemble_partial
#' @usage \S4method{ensemble_partial}{ANY}(...)
methods::setMethod(
  "ensemble_partial",
  methods::signature("ANY"),
  function(..., x.var, method = "mean", layer = "mean", normalize = TRUE){
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

    if(getOption("ibis.setupmessages")) myLog("[Inference]","green","Creating a partial ensemble...")

    # Get variable range from the first object
    # FIXME: Ideally make a consensus, otherwise assumes that same predictor been used
    rr <- range(mods[[1]]$model$predictors[,x.var], na.rm = TRUE)
    assertthat::assert_that(length(rr)==2, !anyNA(rr))
    rr <- seq(rr[1], rr[2], length.out = 100)

    # Now for each object get the partial values for the target variable
    out <- data.frame()
    for(obj in mods){
      if(length(grep(x.var, summary(obj)[[1]]))==0){
        message(paste("Layer", text_red(layer), "not found in model. Skipping!"))
        next()
      }
      # Get partial with identical variable length
      o <- partial(mod = obj, x.var = x.var, variable_length = 100, values = rr, plot = FALSE)
      assertthat::assert_that(all( o$partial_effect == rr ))
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

    # Now composite the ensemble depending on the option
    if(method == 'mean'){
      new <- aggregate(out[,layer], by = list(partial_effect = out$partial_effect),
                                  FUN = function(x = out[[layer]]) {
                                    return(cbind( mean = mean(x), sd = sd(x)))
                                    }) |> as.matrix() |> as.data.frame()
      colnames(new) <- c("partial_effect", "mean", "sd")
    } else if(method == 'median'){
      new <- aggregate(out[,layer], by = list(partial_effect = out$partial_effect),
                       FUN = function(x = out[[layer]]) {
                         return(cbind( median = stats::median(x), mad = mad(x)))
                       }) |> as.matrix() |> as.data.frame()
      colnames(new) <- c("partial_effect", "median", "mad")
    }
    return(new)
  }
)
