#' @include utils.R utils-spatial.R
NULL

#' Spatial adjustment of environmental predictors and raster stacks
#'
#' @description
#' This function allows the transformation of provided environmental predictors (in [`Raster`] format).
#' A common use case is for instance the standardization (or scaling) of all predictors prior to model fitting.
#' This function works both with [`Raster`] as well as with [`stars`] objects.
#' @details
#' Available options are:
#' * \code{'none'} The original layer(s) are returned.
#' * \code{'scale'} This run the [`scale()`] function with default settings (1 Standard deviation) across all predictors.
#' A sensible default to for most model fitting.
#' * \code{'norm'} This normalizes all predictors to a range from \code{0-1}.
#' * \code{'windsor'} This applies a 'windsorization' to an existing raster layer by setting the lowest, respectively
#' largest values to the value at a certain percentage level (e.g. 95%). Those can be set via the parameter \code{"windsor_props"}.
#' * \code{'windsor_thresh'} Same as option 'windsor', however in this case values are clamped to a thresholds
#' rather than certain percentages calculated on the data.
#' * \code{'percentile'} This converts and bins all values into percentiles, e.g. the top 10% or lowest 10% of values and so on.
#' * \code{'pca'} This option runs a principal component decomposition of all predictors (via [`prcomp()`]).
#' It returns new predictors resembling all components in order of the most important ones. Can be useful to
#' reduce collinearity, however note that this changes all predictor names to 'PCX', where X is the number of the component.
#' The parameter \code{'pca.var'} can be modified to specify the minimum variance to be covered by the axes.
#' * \code{'revjack'} Removes outliers from the supplied stack via a reverse jackknife procedure.
#' Identified outliers are by default set to \code{NA}.
#'
#' @param env A [`Raster`] object.
#' @param option A [`vector`] stating whether predictors should be preprocessed in any way (Options: \code{'none'},
#' \code{'scale'}, \code{'norm'}, \code{'windsor'}, \code{'windsor_thresh'}, \code{'percentile'} \code{'pca'}, \code{'revjack'}). See Details.
#' @param windsor_props A [`numeric`] vector specifying the proportions to be clipped for windsorization (Default: \code{c(.05,.95)}).
#' @param pca.var A [`numeric`] value between \code{>0} and \code{1} stating the minimum amount of variance to be covered (Default: \code{0.8}).
#' @param method As \code{'option'} for more intuitive method setting. Can be left empty (in this case option has to be set).
#' @param ... other options (Non specified).
#' @returns Returns a adjusted [`Raster`] object of identical resolution.
#' @seealso predictor_derivate
#' @examples
#' \dontrun{
#' # Where x is a rasterstack
#' new_x <- predictor_transform(x, option = 'scale')
#' }
#' @keywords utils
#' @export
predictor_transform <- function(env, option, windsor_props = c(.05,.95), pca.var = 0.8, method = NULL, ...){
  assertthat::assert_that(
    inherits(env,'Raster') || inherits(env, 'stars'),
    # Support multiple options
    is.numeric(windsor_props) & length(windsor_props)==2,
    is.numeric(pca.var)
  )
  # Convenience function
  if(missing(option)){
    assertthat::assert_that(!is.null(method))
    option <- method
  }
  assertthat::assert_that(
    is.character(option),
    base::length(option) == 1
  )
  # Match option
  option <- match.arg(option, c('none','pca', 'scale', 'norm','windsor', 'windsor_thresh', 'revjack', 'percentile'), several.ok = FALSE)

  # Nothing to be done
  if(option == 'none') return(env)

  # If stars see if we can convert it to a stack
  if(inherits(env, 'stars')){
    lyrs <- names(env) # Names of predictors
    times <- stars::st_get_dimension_values(env, which = 3) # Assume this being the time attribute
    dims <- stars::st_dimensions(env)
    # Convert to list
    env_list <- list()
    for(name in lyrs) env_list[[name]] <- methods::as(env[name], 'Raster')
  } else {
    # Get times in case a stack is supplied (this can get lost depending on transformation)
    times <- raster::getZ(env)
  }

  # Normalization
  if(option == 'norm'){
    if(is.Raster(env)){
      out <- (env - raster::cellStats(env, stat="min")) /
        (raster::cellStats(env, stat="max") -
           raster::cellStats(env, stat="min"))
    } else {
      out <- lapply(env_list, function(x) {
        (x - raster::cellStats(x, stat="min")) /
          (raster::cellStats(x, stat="max") -
             raster::cellStats(x, stat="min"))
      })
    }
  }
  # Scaling
  if(option == 'scale'){
    if(is.Raster(env)){
      out <- raster::scale(env, center = TRUE, scale = TRUE)
    } else {
      out <- lapply(env_list, function(x) raster::scale(x, center = TRUE, scale = TRUE))
    }
  }

  # Percentile cutting
  if(option == 'percentile'){
    if(is.Raster(env)){
      perc <- raster::quantile(env, seq(0,1, length.out = 11))
      perc <- unique(perc)
      out <- raster::cut(env, perc)
    } else {
      out <- lapply(env_list, function(x) {
        perc <- raster::quantile(x, seq(0,1, length.out = 11))
        perc <- unique(perc)
        raster::cut(x, perc)
      })
    }
  }

  # Windsorization
  if(option == 'windsor'){
    win <- function(x, windsor_props){
      xq <- stats::quantile(x = x[], probs = windsor_props, na.rm = TRUE)
      min.value <- xq[1]
      max.value <- xq[2]
      if(is.vector(env)) out <- units::drop_units(env) else out <- env
      out[out < min.value] <- min.value
      out[out > max.value] <- max.value
      out
    }
    if(is.Raster(env)){
      out <- win(env, windsor_props )
    } else {
      out <- lapply(env_list, function(x) win(x, windsor_props))
    }
  } else if(option == 'windsor_thresh'){
    win <- function(x, windsor_thresh){
      if(is.vector(env)) out <- units::drop_units(env) else out <- env
      out[out < windsor_thresh[1]] <- windsor_thresh[1]
      out[out > windsor_thresh[2]] <- windsor_thresh[2]
      out
    }
    if(is.Raster(env)){
      out <- win(env, windsor_props )
    } else {
      out <- lapply(env_list, function(x) win(x, windsor_props))
    }
  }

  # Reverse jackknife removal of outliers
  if(option == 'revjack'){
    rj <- function(x){
      o <- emptyraster(x)
      o[] <- rm_outlier_revjack(x[], procedure = "missing")
      return(o)
    }
    if(is.Raster(env)){
      out <- raster::stack()
      for(n in 1:nlayers(env)){
        out <- raster::addLayer(out, rj(env[[n]]) )
      }
    } else {
      out <- lapply(env_list, function(x) rj(x))
    }
  }

  # Principle component separation of variables
  # Inspiration taken from RSToolbox package
  if(option == 'pca'){
    if(is.Raster(env)){
      assertthat::assert_that(raster::nlayers(env)>=2,msg = 'Need at least two predictors to calculate PCA.')

      # FIXME: Allow a reduction to few components than nr of layers?
      nComp <- nlayers(env)
      # Construct mask of all cells
      envMask <- !sum(raster::calc(env, is.na))
      assertthat::assert_that(cellStats(envMask, sum)>0,msg = 'A predictor is either NA only or no valid values across all layers')
      env <- raster::mask(env, envMask, maskvalue = 0)

      # Sample covariance from stack and fit PCA
      covMat <- raster::layerStats(env, stat = "cov", na.rm = TRUE)
      pca <- stats::princomp(covmat = covMat[[1]], cor = FALSE)
      # Add means and grid cells
      pca$center <- covMat$mean
      pca$n.obs <- raster::ncell(env)

      # Check how many components are requested:
      if(pca.var<1){
        sums <- stats::loadings( summary(pca) )[]
        props <- cumsum(colSums(sums^2) / nrow(sums)) # Cumulative explained variance
        nComp <- length( which(props <= pca.var) )
      }
      # Predict principle components
      out <- raster::predict(env, pca,na.rm = TRUE, index = 1:nComp)
      names(out) <- paste0("PC", 1:nComp)

      return(out)
    } else {
      # TODO:
      stop("Principal component transformation for stars objects is not yet implemented. Pre-process externally!")
    }
  }

  # If stars convert back to stars object
  if(inherits(env, 'stars')){
    # Convert list back to stars
    out <- do.call(
      stars:::c.stars,
      lapply(out, function(x) stars::st_as_stars(x))
    )
    # Reset names of attributes
    names(out) <- lyrs
    # FIXME: Hacky solution, but breaks other scenarios otherwise
    out2 <- try({stars::st_set_dimensions(out, which = 3, values = times, names = "time")},silent = TRUE)
    if(inherits(out2, "try-error")){
      # This happens when a stars provided layer has only a single time band
      out <- stars::st_redimension(out, new_dims = dims) # use the previously saved dimensions
    } else { out <- out2; out2}
  } else {
    # Final security checks
    assertthat::assert_that(
      raster::nlayers(env) == raster::nlayers(out),
      is_comparable_raster(out, env)
    )
    # Reset times
    if(!is.null(times)) out <- raster::setZ(out, times)

    return(out)
  }
}

#' Create spatial derivative of raster stacks
#'
#' @description
#' This function creates derivatives of existing covariates and returns them in Raster format.
#' Derivative variables can in the machine learning literature commonly be understood as one aspect of feature
#' engineering. They can be particularly powerful in introducing non-linearities in otherwise linear models,
#' for example is often done in the popular Maxent framework.
#' @details
#' Available options are:
#' * \code{'none'} - The original layer(s) are returned.
#' * \code{'quadratic'} - A quadratic transformation (\eqn{x^{2}}) is created of the provided layers.
#' * \code{'hinge'} - Creates hinge transformation of covariates, which set all values lower than a set threshold to \code{0}
#' and all others to a range of \eqn{[0,1]}. The number of thresholds and thus new derivates is specified
#' via the parameter \code{'nknots'} (Default: \code{4}).
#' * \code{'interaction'} - Creates interactions between variables. Target variables have to be specified via \code{"int_variables"}.
#' * \code{'thresh'} - A threshold transformation of covariates, which sets all values lower than a set threshold ot
#' \code{0} and those larger to \code{1}.
#' The number of thresholds and thus new derivates is specified via the parameter \code{'nknots'} (Default: \code{4}).
#' * \code{'bin'} - Creates a factor representation of a covariates by cutting the range of covariates by their percentiles.
#' The number of percentile cuts and thus new derivates is specified via the parameter \code{'nknots'} (Default: \code{4}).
#' @param env A [`Raster`] object.
#' @param option A [`vector`] stating whether predictors should be preprocessed in any way
#' (Options: \code{'none'}, \code{'quadratic'}, \code{'hinge'}, \code{'thresh'}, \code{'bin'}).
#' @param nknots The number of knots to be used for the transformation (Default: \code{4}).
#' @param deriv A [`vector`] with [`characters`] of specific derivates to create (Default: \code{NULL}).
#' @param int_variables A [`vector`] with length greater or equal than \code{2} specifying the covariates  (Default: \code{NULL}).
#' @param method As \code{'option'} for more intuitive method setting. Can be left empty (in this case option has to be set).
#' @param ... other options (Non specified).
#' @return Returns the derived adjusted [`Raster`] objects of identical resolution.
#' @seealso predictor_derivate
#' @examples
#' \dontrun{
#' # Create a hinge transformation of one or multiple RasterLayers.
#' predictor_derivate(covs, option = "hinge", knots = 4)
#' }
#' @keywords utils
#' @export
predictor_derivate <- function(env, option, nknots = 4, deriv = NULL, int_variables = NULL, method = NULL, ...){
  assertthat::assert_that(
    inherits(env,'Raster') || inherits(env, "stars"),
    !missing(env),
    is.numeric(nknots) && nknots > 1,
    is.null(deriv) || is.character(deriv),
    is.null(int_variables) || is.vector(int_variables)
  )
  # Convenience function
  if(missing(option)){
    assertthat::assert_that(!is.null(method))
    option <- method
  }
  assertthat::assert_that(
    is.character(option),
    base::length(option) == 1
  )
  # Match argument.
  option <- match.arg(option, c('none','quadratic', 'hinge', 'thresh', 'bin', 'interaction'), several.ok = FALSE)

  # None, return as is
  if(option == 'none') return(env)

  # If stars see if we can convert it to a stack
  if(inherits(env, 'stars')){
    assertthat::assert_that(!is.null(deriv),msg = "Derivate names could not be found!")
    # Decompose derivate variable names if set
    deriv <- grep(paste0(option, "__"), deriv, value = TRUE)
    if(length(deriv)==0){
      if(getOption('ibis.setupmessages')) myLog('[Setup]','red','Predictors with derivates not found!')
      return(NULL)
    }
    cutoffs <- do.call(rbind,strsplit(deriv, "__")) |> as.data.frame()
    cutoffs$deriv <- deriv

    lyrs <- names(env) # Names of predictors
    times <- stars::st_get_dimension_values(env, which = 3) # Time attribute
    # Create a list to house the results
    env_list <- list()
    for(name in cutoffs$deriv){
      env_list[[name]] <- methods::as(env[cutoffs[which(cutoffs$deriv==name),2]], 'Raster') # Specify original raster
    }
    assertthat::assert_that(length(env_list) > 0)
  } else {cutoffs <- NULL}

  # Simple quadratic transformation
  if(option == 'quadratic'){
    if(is.Raster(env)){
      if(raster::nlayers(env)==1){
        new_env <- env^2
      } else {
        new_env <- raster::calc(env, function(x) I(x^2))
      }
      names(new_env) <- paste0('quad__', names(env))
    } else {
      # Stars processing
      new_env <- lapply(env_list, function(x) {
        raster::calc(x, function(z) I(z^2))
      })
    }
  }

  # Hinge transformation
  # From`maxnet` package
  if(option == 'hinge'){
    if(is.Raster(env)){
      # Build new stacks
      new_env <- raster::stack()
      for(val in names(env)){
        o <- makeHinge(env[[val]], n = val, nknots = nknots, cutoffs = cutoffs)
        if(is.null(o)) next()
        new_env <- raster::addLayer(new_env,
                                    fill_rasters(o, emptyraster(env) )
        )
        rm(o)
      }
    } else {
      # Stars object
      for(val in names(env_list)){
        # Format cutoffs
        cu <- cutoffs[which(cutoffs$deriv == val), 3]
        cu <- strsplit(cu, "_") |> unlist()
        # Remove any leading points
        if(any(substr(cu,1, 1)==".")){
          cu[which(substr(cu,1, 1)==".")] <- gsub("^.","",cu[which(substr(cu,1, 1)==".")])
        }
        cu <- as.numeric(cu)
        assertthat::assert_that(!anyNA(cu), is.numeric(cu))
        for(k in 1:nlayers(env_list[[val]])){
          o <- emptyraster(env_list[[val]][[k]])
          o[] <- hingeval(env_list[[val]][[k]][], cu[1], cu[2])
          env_list[[val]][[k]] <- o
          rm(o)
        }
      }
      invisible(gc())
    }
  }

  # For thresholds
  # Take functionality in maxnet package
  if(option == 'thresh'){
    if(is.Raster(env)){
      new_env <- raster::stack()
      for(val in names(env)){
        o <- makeThresh(env[[val]],n = val,nknots = nknots, cutoffs = cutoffs)
        if(is.null(o)) next()
        new_env <- raster::addLayer(new_env,
                                    fill_rasters(o, emptyraster(env))
        )
        rm(o)
      }
    } else {
      # For stats layers
      for(val in names(env_list)){
        # Format cutoffs
        cu <- cutoffs[which(cutoffs$deriv == val), 3]
        cu <- strsplit(cu, "_") |> unlist()
        # Remove any leading points
        if(any(substr(cu,1, 1)==".")){
          cu[which(substr(cu,1, 1)==".")] <- gsub("^.","",cu[which(substr(cu,1, 1)==".")])
        }
        cu <- as.numeric(cu)
        assertthat::assert_that(!anyNA(cu), is.numeric(cu))
        for(k in 1:nlayers(env_list[[val]])){
          o <- emptyraster(env_list[[val]][[k]])
          o[] <- thresholdval(env_list[[val]][[k]][], cu)
          env_list[[val]][[k]] <- o
          rm(o)
        }
      }
      invisible(gc())
    }
  }

  # For binning, calculate cuts of thresholds
  if(option == 'bin'){
    if(is.Raster(env)){
      new_env <- raster::stack()
      for(val in names(env)){
        o <- makeBin(env[[val]],n = val,nknots = nknots, cutoffs = cutoffs)
        if(is.null(o)) next()
        new_env <- raster::addLayer(new_env, o)
        rm(o)
      }
    } else {
      # For stats layers
      for(val in names(env_list)){
        # Format cutoffs
        cu <- cutoffs[which(cutoffs$deriv == val), 3]
        cu <- strsplit(cu, "_") |> unlist()
        # Remove any leading points
        if(any(substr(cu,1, 1)==".")){
          cu[which(substr(cu,1, 1)==".")] <- gsub("^.","",cu[which(substr(cu,1, 1)==".")])
        }
        cu <- as.numeric(cu)
        assertthat::assert_that(!anyNA(cu), is.numeric(cu))
        for(k in 1:nlayers(env_list[[val]])){
          o <- emptyraster(env_list[[val]][[k]])
          o <- raster::cut(env_list[[val]][[k]], cu)
          o[is.na(o)] <- 0
          o <- raster::mask(o, env_list[[val]][[k]] )
          env_list[[val]][[k]] <- o
          rm(o)
        }
      }
      invisible(gc())
    }
  }

  # Create interaction variables
  if(option == 'interaction'){
    # Check whether interaction is provided or an attribute
    if(is.null(int_variables)){
      int_variables <- attr(env, "int_variables")
    }
    assertthat::assert_that(is.vector(int_variables))

    if(is.Raster(env)){
      # Make unique combinations
      ind <- utils::combn(int_variables, 2)

      # Now for each combination build new variable
      new_env <- raster::stack()

      for(i in 1:ncol(ind)){
        # Multiply first with second entry
        o <- env[[ind[1,1]]] * env[[ind[2,1]]]
        names(o) <- paste0('inter__', names(env)[ind[1,1]],".",names(env)[ind[2,1]])
        new_env <- raster::addLayer(new_env,o)
        rm(o)
      }
    } else {
      # Stars processing
      stop("Not yet implemented!")
    }
  }

  # If stars convert back to stars object
  if(inherits(env, 'stars')){
    # Add the original layers back
    for(name in names(env)){
      env_list[[name]] <- methods::as(env[name], 'Raster') # Specify original raster
    }

    # Convert list back to stars
    new_env <- do.call(
      stars:::c.stars,
      lapply(env_list, function(x) stars::st_as_stars(x))
    )
    # Reset names of attributes
    names(new_env) <- c( cutoffs$deriv, names(env))
    new_env <- stars::st_set_dimensions(new_env, which = 3, values = times, names = "time")
  }
  return(new_env)
}

#' Homogenize NA values across a set of predictors.
#'
#' @description This method allows the homogenization of missing data across a set of environmental predictors.
#' It is by default called when predictors are added to [´BiodiversityDistribution´] object. Only grid cells with NAs that contain
#' values at some raster layers are homogenized.
#' Additional parameters allow instead of homogenization to fill the missing data with neighbouring values
#' @param env A [`Raster`] object with the predictors
#' @param fill A [`logical`] value indicating whether missing data are to be filled (Default: FALSE).
#' @param fill_method A [`character`] of the method for filling gaps to be used (Default: 'ngb')
#' @param return_na_cells A [`logical`] value of whether the ids of grid cells with NA values is to be returned instead (Default: FALSE)
#' @returns A [`Raster`] object with the same number of layers as the input.
#' @keywords utils
#' @export
predictor_homogenize_na <- function(env, fill = FALSE, fill_method = 'ngb', return_na_cells = FALSE){
  assertthat::assert_that(
    is.Raster(env) || inherits(env, 'stars'),
    is.logical(fill),
    is.character(fill_method), fill_method %in% c('ngb'),
    is.logical(return_na_cells)
  )
  # Workflow for raster layers
  if(is.Raster(env)){
    nl <- raster::nlayers(env)
    # If the number of layers is 1, no need for homogenization
    if(nl > 1){
      # Calculate number of NA grid cells per stack
      mask_na <- sum( is.na(env) )
      # Remove grid cells that are equal to the number of layers (all values NA)
      none_area <- mask_na == nl
      none_area[none_area == 0 ] <- NA
      mask_na <- raster::mask(mask_na,
                              mask = none_area,inverse = TRUE)

      # Should any fill be conducted?
      if(fill){
        stop('Not yet implemented!')
      } else {
        # Otherwise just homogenize NA values across predictors
        if(cellStats(mask_na,'max')>0){
          mask_all <- mask_na == 0; mask_all[mask_all == 0] <- NA
          env <- raster::mask(env, mask = mask_all)
        }
      }
      # Should NA coordinates of cells where 1 or more predictor is NA be returned?
      # FIXME: One could directly return a data.frame with the predictor names to allow easier lookup.
      if(return_na_cells){
        vals <- which((mask_na>0)[])
        env <- list(cells_na = vals, env = env)
      }
      rm(mask_na, none_area) # Cleanup
    }
  } else if(inherits(env, 'stars')){
    stop('Not implemented yet.')
  }
  # Security checks
  assertthat::assert_that(
    is.Raster(env) || is.list(env) || inherits(env, 'stars')
  )
  # Return the result
  return(env)
}
