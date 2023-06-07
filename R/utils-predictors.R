#' @include utils.R utils-spatial.R
NULL

#' Spatial adjustment of environmental predictors and raster stacks
#'
#' @description
#' This function allows the transformation of provided environmental predictors (in [`SpatRaster`] format).
#' A common use case is for instance the standardization (or scaling) of all predictors prior to model fitting.
#' This function works both with [`SpatRaster`] as well as with [`stars`] objects.
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
#' @param env A [`SpatRaster`] object.
#' @param option A [`vector`] stating whether predictors should be preprocessed in any way (Options: \code{'none'},
#' \code{'scale'}, \code{'norm'}, \code{'windsor'}, \code{'windsor_thresh'}, \code{'percentile'} \code{'pca'}, \code{'revjack'}). See Details.
#' @param windsor_props A [`numeric`] vector specifying the proportions to be clipped for windsorization (Default: \code{c(.05,.95)}).
#' @param pca.var A [`numeric`] value between \code{>0} and \code{1} stating the minimum amount of variance to be covered (Default: \code{0.8}).
#' @param method As \code{'option'} for more intuitive method setting. Can be left empty (in this case option has to be set).
#' @param ... other options (Non specified).
#' @returns Returns a adjusted [`SpatRaster`] object of identical resolution.
#' @seealso predictor_derivate
#' @examples
#' \dontrun{
#' # Where x is a SpatRaster
#' new_x <- predictor_transform(x, option = 'scale')
#' }
#' @keywords utils
#' @export
predictor_transform <- function(env, option, windsor_props = c(.05,.95), pca.var = 0.8, method = NULL, ...){
  assertthat::assert_that(
    is.Raster(env) || inherits(env, 'stars'),
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
    # Drop units as it causes bugs with terra (5/5/2023)
    env <- units::drop_units(env)
    # Convert to list
    env_list <- list()
    for(name in lyrs) env_list[[name]] <- Reduce(c, stars_to_raster( env[name] ))

    # Make quick checks
    assertthat::assert_that( all( diff(sapply(env_list, terra::nlyr))==0 ) )
  } else {
    # Get times in case a stack is supplied (this can get lost depending on transformation)
    times <- terra::time(env)
  }

  # Normalization
  if(option == 'norm'){
    if(is.Raster(env)){
      nx <- terra::minmax(env)
      out <- (env - nx[1,]) / (nx[2,] - nx[1,])
    } else {
      out <- lapply(env_list, function(x) {
        nx <- terra::minmax(x)
        (x - nx[1,]) / (nx[2,] - nx[1,])
      })
    }
  }
  # Scaling
  if(option == 'scale'){
    if(is.Raster(env)){
      out <- terra::scale(env, center = TRUE, scale = TRUE)
    } else {
      out <- lapply(env_list, function(x) terra::scale(x, center = TRUE, scale = TRUE))
    }
  }

  # Percentile cutting
  if(option == 'percentile'){
    if(is.Raster(env)){
      perc <- terra::global(env, fun = function(z) terra::quantile(z, probs = seq(0,1, length.out = 11), na.rm = TRUE))
      perc <- unique(perc)
      out <- terra::classify(env, t(perc))
    } else {
      out <- lapply(env_list, function(x) {
        perc <- terra::global(x, fun = function(z) terra::quantile(z, probs = seq(0,1, length.out = 11), na.rm = TRUE))
        perc <- unique(perc)
        # For terra need to loop here as classify does not support multiple columns
        o <- terra::rast()
        for(i in 1:nrow(perc)) o <- suppressWarnings( c(o, terra::classify(x[[i]], rcl = t(perc)[,i]) ))
        return(o)
      })
      assertthat::assert_that( all( sapply(out, function(z) all(is.factor(z))) ))
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
      out <- terra::rast()
      for(n in 1:terra::nlyr(env)){
        suppressWarnings( out <- c(out, rj(env[[n]]) ) )
      }
    } else {
      out <- lapply(env_list, function(x) rj(x))
    }
  }

  # Principle component separation of variables
  # Inspiration taken from RSToolbox package
  if(option == 'pca'){
    if(is.Raster(env)){
      assertthat::assert_that(terra::nlyr(env)>=2,msg = 'Need at least two predictors to calculate PCA.')

      # Check that there are no duplicates in the layer names, if so append numbers to them
      if(anyDuplicated(names(env))>0) names(env) <- make.unique(names(env))

      # FIXME: Allow a reduction to few components than nr of layers?
      nComp <- terra::nlyr(env)
      # Construct mask of all cells
      envMask <- !sum(terra::app(env, is.na))
      assertthat::assert_that(terra::global(envMask, "sum")>0,msg = 'A predictor is either NA only or no valid values across all layers')
      env <- terra::mask(env, envMask, maskvalues = 0)

      # Sample covariance from stack and fit PCA
      covMat <- terra::layerCor(env, fun = "cov", na.rm = TRUE)
      pca <- stats::princomp(covmat = covMat[[1]], cor = FALSE)
      # Add means and grid cells
      pca$center <- covMat$mean
      pca$n.obs <- terra::ncell(env)

      # Check how many components are requested:
      if(pca.var<1){
        sums <- stats::loadings( summary(pca) )[]
        props <- cumsum(colSums(sums^2) / nrow(sums)) # Cumulative explained variance
        nComp <- length( which(props <= pca.var) )
      }
      # Predict principle components
      out <- terra::predict(env, pca,na.rm = TRUE, index = 1:nComp)
      names(out) <- paste0("PC", 1:nComp)

      return(out)
    } else {
      # TODO:
      stop("Principal component transformation for stars objects is not yet implemented. Pre-process externally!")
    }
  }

  # If stars convert back to stars object
  if(inherits(env, 'stars')){
    # Set crs to all
    if(!is.na( sf::st_crs(env) )){
      for(i in 1:length(out)) terra::set.crs(out[[i]], crs(sf::st_crs(env)$wkt))
    }
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
      terra::nlyr(env) == terra::nlyr(out),
      is_comparable_raster(out, env)
    )
    # Reset times
    if(!all(is.na(terra::time(env))) ){
      if(!is.null(times)) terra::time(out) <- times
    }

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
#' @param env A [`SpatRaster`] object.
#' @param option A [`vector`] stating whether predictors should be preprocessed in any way
#' (Options: \code{'none'}, \code{'quadratic'}, \code{'hinge'}, \code{'thresh'}, \code{'bin'}).
#' @param nknots The number of knots to be used for the transformation (Default: \code{4}).
#' @param deriv A [`vector`] with [`characters`] of specific derivates to create (Default: \code{NULL}).
#' @param int_variables A [`vector`] with length greater or equal than \code{2} specifying the covariates  (Default: \code{NULL}).
#' @param method As \code{'option'} for more intuitive method setting. Can be left empty (in this case option has to be set).
#' @param ... other options (Non specified).
#' @return Returns the derived adjusted [`SpatRaster`] objects of identical resolution.
#' @seealso predictor_derivate
#' @examples
#' \dontrun{
#' # Create a hinge transformation of one or multiple SpatRaster.
#' predictor_derivate(covs, option = "hinge", knots = 4)
#' }
#' @keywords utils
#' @export
predictor_derivate <- function(env, option, nknots = 4, deriv = NULL, int_variables = NULL, method = NULL, ...){
  assertthat::assert_that(
    is.Raster(env) || inherits(env, "stars"),
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
      env_list[[name]] <- terra::rast(env[cutoffs[which(cutoffs$deriv==name),2]]) # Specify original raster
    }
    assertthat::assert_that(length(env_list) > 0)
  } else {cutoffs <- NULL}

  # Simple quadratic transformation
  if(option == 'quadratic'){
    if(is.Raster(env)){
      if(terra::nlyr(env)==1){
        new_env <- env^2
      } else {
        new_env <- terra::app(env, function(x) I(x^2))
      }
      names(new_env) <- paste0('quad__', names(env))
    } else {
      # Stars processing
      new_env <- lapply(env_list, function(x) {
        terra::app(x, function(z) I(z^2))
      })
    }
  }

  # Hinge transformation
  # From`maxnet` package
  if(option == 'hinge'){
    if(is.Raster(env)){
      # Build new stacks
      new_env <- terra::rast()
      for(val in names(env)){
        o <- makeHinge(env[[val]], n = val, nknots = nknots, cutoffs = cutoffs)
        if(is.null(o)) next()
        new_env <- c(new_env, fill_rasters(o, emptyraster(env) ) )
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
        for(k in 1:terra::nlyr(env_list[[val]])){
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
      new_env <- terra::rast()
      for(val in names(env)){
        o <- makeThresh(env[[val]],n = val,nknots = nknots, cutoffs = cutoffs)
        if(is.null(o)) next()
        new_env <- c(new_env, fill_rasters(o, emptyraster(env)) )
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
        for(k in 1:terra::nlyr(env_list[[val]])){
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
      new_env <- terra::rast()
      for(val in names(env)){
        o <- makeBin(env[[val]], n = val, nknots = nknots, cutoffs = cutoffs)
        if(is.null(o)) next()
        new_env <- c(new_env, o)
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
        for(k in 1:terra::nlyr(env_list[[val]])){
          o <- emptyraster(env_list[[val]][[k]])
          o <- terra::classify(env_list[[val]][[k]], cu)
          o[is.na(o)] <- 0
          o <- terra::mask(o, env_list[[val]][[k]] )
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
      new_env <- terra::rast()

      for(i in 1:ncol(ind)){
        # Multiply first with second entry
        o <- env[[ind[1,i]]] * env[[ind[2,i]]]
        names(o) <- paste0('inter__', ind[1, i],".", ind[2, i])
        new_env <- c(new_env, o)
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
      env_list[[name]] <- terra::rast(env[name]) # Specify original raster
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
#' @param env A [`SpatRaster`] object with the predictors.
#' @param fill A [`logical`] value indicating whether missing data are to be filled (Default: \code{FALSE}).
#' @param fill_method A [`character`] of the method for filling gaps to be used (Default: \code{'ngb'}).
#' @param return_na_cells A [`logical`] value of whether the ids of grid cells with NA values is to be returned instead (Default: \code{FALSE}).
#' @returns A [`SpatRaster`] object with the same number of layers as the input.
#' @examples
#' \dontrun{
#'  # Harmonize predictors
#'  env <- predictor_homogenize_na(env)
#' }
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
    nl <- terra::nlyr(env)

    if(anyDuplicated(names(env))>0) names(env) <- make.unique(names(env))

    # If the number of layers is 1, no need for homogenization
    if(nl > 1){
      # Calculate number of NA grid cells per stack
      mask_na <- sum( is.na(env) )
      # Remove grid cells that are equal to the number of layers (all values NA)
      none_area <- mask_na == nl
      none_area[none_area == 0 ] <- NA
      mask_na <- terra::mask(mask_na, mask = none_area, inverse = TRUE)

      # Should any fill be conducted?
      if(fill){
        stop('Not yet implemented!')
      } else {
        # Otherwise just homogenize NA values across predictors
        if(terra::global(mask_na,'max',na.rm = TRUE)[,1] > 0){
          mask_all <- mask_na == 0; mask_all[mask_all == 0] <- NA
          env <- terra::mask(env, mask = mask_all)
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

#' Hinge transformation of a given predictor
#'
#' @description
#' This function transforms a provided predictor variable with a hinge transformation,
#' e.g. a new range of values where any values lower than a certain knot are set to \code{0},
#' while the remainder is left at the original values.
#' @param v A [`SpatRaster`] object.
#' @param n A [`character`] describing the name of the variable. Used as basis for new names.
#' @param nknots The number of knots to be used for the transformation (Default: \code{4}).
#' @param cutoffs A [`numeric`] vector of optionally used cutoffs to be used instead (Default: \code{NULL}).
#' @keywords utils, internal
#' @concept Concept taken from the [maxnet] package.
#' @returns A hinge transformed [`data.frame`].
#' @noRd
makeHinge <- function(v, n, nknots = 4, cutoffs = NULL){
  assertthat::assert_that(is.Raster(v),
                          is.character(n),
                          is.numeric(nknots),
                          is.numeric(cutoffs) || is.null(cutoffs))
  # Get stats
  v.min <- terra::global(v, "min", na.rm = TRUE)[,1]
  v.max <- terra::global(v, "max", na.rm = TRUE)[,1]
  assertthat::assert_that(is.numeric(v.min), is.numeric(v.max))

  if(is.null(cutoffs)){
    k <- seq(v.min, v.max, length = nknots)
  } else {
    k <- cutoffs
  }
  if(length(k)<=1) return(NULL)

  # Hinge up to max
  lh <- outer(v[] |> as.vector(), utils::head(k, -1), function(w, h) hingeval(w,h, v.max))
  # Hinge starting from min
  rh <- outer(v[] |> as.vector(), k[-1], function(w, h) hingeval(w, v.min, h))
  colnames(lh) <- paste0("hinge__",n,'__', round( utils::head(k, -1), 2),'_', round(v.max, 2))
  colnames(rh) <- paste0("hinge__",n,'__', round( v.min, 2),'_', round(k[-1], 2))
  o <- as.data.frame(
    cbind(lh, rh)
  )
  # Kick out first (min) and last (max) col as those are perfectly correlated
  o <- o[,-c(1,ncol(o))]
  attr(o, "deriv.hinge") <- k
  return(o)
}

#' Threshold transformation of a given predictor
#'
#' @description
#' This function transforms a provided predictor variable with a threshold transformation,
#' e.g. a new range of values where any values lower than a certain knot are set to \code{0},
#' while the remainder is set to \code{1}.
#' @param v A [`Raster`] object.
#' @param n A [`character`] describing the name of the variable. Used as basis for new names.
#' @param nknots The number of knots to be used for the transformation (Default: \code{4}).
#' @param cutoffs A [`numeric`] vector of optionally used cutoffs to be used instead (Default: \code{NULL}).
#' @keywords utils, internal
#' @concept Concept taken from the [maxnet] package.
#' @returns A threshold transformed [`data.frame`].
#' @noRd
makeThresh <- function(v, n, nknots = 4, cutoffs = NULL){
  assertthat::assert_that(is.Raster(v),
                          is.character(n),
                          is.numeric(nknots),
                          is.numeric(cutoffs) || is.null(cutoffs))
  if(is.null(cutoffs)){
    # Get min max
    v.min <- terra::global(v, "min", na.rm= TRUE)[,1]
    v.max <- terra::global(v, "max", na.rm= TRUE)[,1]
    k <- seq(v.min, v.max, length = nknots + 2)[2:nknots + 1]
  } else {
    k <- cutoffs
  }
  if(length(k)<=1) return(NULL)
  f <- outer(v[] |> as.vector(), k, function(w, t) ifelse(w >= t, 1, 0))
  colnames(f) <- paste0("thresh__", n, "__",  round(k, 2))
  f <- as.data.frame(f)
  attr(f, "deriv.thresh") <- k
  return(f)
}

#' Binned transformation of a given predictor
#'
#' @description
#' This function takes predictor values and 'bins' them into categories based on a
#' percentile split.
#' @param v A [`SpatRaster`] object.
#' @param n A [`character`] describing the name of the variable. Used as basis for new names.
#' @param nknots The number of knots to be used for the transformation (Default: \code{4}).
#' @param cutoffs A [`numeric`] vector of optionally used cutoffs to be used instead (Default: \code{NULL}).
#' @keywords utils, internal
#' @returns A binned transformed [`data.frame`] with columns representing each bin.
#' @noRd
makeBin <- function(v, n, nknots, cutoffs = NULL){
  assertthat::assert_that(is.Raster(v),
                          is.character(n),
                          is.numeric(nknots),
                          is.numeric(cutoffs) || is.null(cutoffs))
  if(is.null(cutoffs)){
    # Calculate cuts
    cu <- terra::quantile(v[], probs = seq(0, 1, by = 1/nknots), na.rm = TRUE )
    assertthat::assert_that(all(is.numeric(cu)), length(cu)>1)
  } else { cu <- cutoffs }

  if(anyDuplicated(cu)){
    # If duplicated quantiles (e.g. 0, 0, 0.2..), sample from a larger number
    cu <- terra::quantile(v[], probs = seq(0, 1, by = 1/(nknots*2)) )
    cu <- cu[-which(duplicated(cu))] # Remove duplicated cuts
    if(length(cu)<=2) return( NULL )
    if(length(cu) > nknots){
      cu <- cu[(length(cu)-(nknots)):length(cu)]
    }
  }
  # Make cuts and explode
  out <- explode_factorized_raster(
      terra::classify(v, cu)
  )

  # Format threshold names
  cu.brk <- as.character(cut(cu[-1], cu))
  cu.brk <- gsub(",","_",cu.brk)
  cu.brk <- gsub("\\(|\\]", "", cu.brk)
  # names(out) <- paste0("bin__",n, "__", gsub(x = names(cu)[-1], pattern = "\\D", replacement = ""),"__", cu.brk )
  names(out) <- paste0("bin__",n, "__", cu.brk )
  for(i in 1:terra::nlyr(out)){
    attr(out[[i]], "deriv.bin") <- cu[i:(i+1)]
  }
  return(out)
}

#### Filter predictor functions ----

#' Filter a set of correlated predictors to fewer ones
#'
#' @description
#' This function helps to remove highly correlated variables from a set of predictors. It supports multiple options
#' some of which require both environmental predictors and observations, others only predictors.
#'
#' Some of the options require different packages to be pre-installed, such as [ranger] or [Boruta].
#'
#' @details
#' Available options are:
#'
#' * \code{"none"} No prior variable removal is performed (Default).
#' * \code{"pearson"}, \code{"spearman"} or \code{"kendall"} Makes use of pairwise comparisons to identify and
#' remove highly collinear predictors (Pearson's \code{r >= 0.7}).
#' * \code{"abess"} A-priori adaptive best subset selection of covariates via the [abess] package (see References). Note that this
#' effectively fits a separate generalized linear model to reduce the number of covariates.
#' * \code{"boruta"} Uses the [Boruta] package to identify non-informative features.
#'
#' @note
#' Using this function on predictors effectively means that a separate model is fitted on the data
#' with all the assumptions that come with in (e.g. linearity, appropriateness of response, normality, etc).
#'
#' @param env A [`data.frame`] or [`matrix`] with extracted environmental covariates for a given species.
#' @param obs A [`vector`] with observational records to use for determining variable importance. Can be \code{NULL}.
#' @param keep A [`vector`] with variables to keep regardless. These are usually variables for which prior
#' information is known.
#' @param method Which method to use for constructing the correlation matrix (Options: \code{'pearson'} (Default),
#'  \code{'spearman'}| \code{'kendal'}), \code{"abess"}, or \code{"boruta"}.
#' @param ... Other options for a specific method
#'
#' @keywords utils
#' @return A [`character`] [`vector`] of variable names to be excluded.
#' If the function fails due to some reason return \code{NULL}.
#' @examples
#' \dontrun{
#'  # Remove highly correlated predictors
#'  env <- predictor_filter( env, option = "pearson")
#' }
#' @export
predictor_filter <- function( env, keep = NULL, method = "pearson", ...){
  assertthat::assert_that(
    is.data.frame(env) || is.matrix(env),
    ncol(env) >2,
    is.null(keep) || is.vector(keep),
    is.character(method)
  )
  # Match the predictor names
  method <- match.arg(method,
                      c("none", "pearson", "spearman", "kendall", "abess", "boruta"),
                      several.ok = FALSE)

  # Now apply the filter depending on the option
  if(method == "none"){
    co <- NULL
  } else if(method %in% c("pearson", "spearman", "kendall")){
    # Simply collinearity check based on colinear predictors
    co <- predictors_filter_collinearity(
      env, keep = keep, method = method, ...
    )
  } else if(method == "abess"){
    if(getOption('ibis.setupmessages')) myLog('[Estimation]','yellow','Applying abess method to reduce predictors...')
    co <- predictors_filter_abess(
      env = env, keep = keep, method = method, ...
    )
  } else if(method == "boruta"){
    check_package("Boruta")
    co <- predictors_filter_boruta(
      env = env, keep = keep, method = method, ...
    )

  } else {
    stop("Method not yet implemented!")
  }

  # Security checks and return
  assertthat::assert_that(is.null(co) || is.character(co))
  return(co)
}

#' Identify collinear predictors
#'
#' @inheritParams predictor_filter
#' @param cutoff A [`numeric`] variable specifying the maximal correlation cutoff.
#' @concept Code inspired from the [`caret`] package
#' @keywords utils, internal
#' @returns [`vector`] of variable names to exclude
predictors_filter_collinearity <- function( env, keep = NULL, cutoff = getOption('ibis.corPred'), method = 'pearson', ...){
  # Security checks
  assertthat::assert_that(is.data.frame(env),
                          is.character(method),
                          is.numeric(cutoff),
                          is.null(keep) || is.vector(keep)
  )
  keep <- keep[keep %in% names(env)] # Remove those not in the data.frame. For instance if a spatial effect is selected
  if(!is.null(keep) || length(keep) == 0) x <- env |> dplyr::select(-keep) else x <- env

  # Removing non-numeric columns
  non.numeric.columns <- colnames(x)[!sapply(x, is.numeric)]
  x <- x[, !(colnames(x) %in% non.numeric.columns)]

  # Get all variables that are singular or unique in value
  singular_var <- which(round( apply(x, 2, var),4) == 0)
  if(length(singular_var)>0) x <- x[,-singular_var]

  # Calculate correlation matrix
  cm <- stats::cor(x, method = method)

  # Copied from the \code{caret} package to avoid further dependencies
  if (any(!stats::complete.cases(cm))) stop("The correlation matrix has some missing values.")
  averageCorr <- colMeans(abs(cm))
  averageCorr <- as.numeric(as.factor(averageCorr))
  cm[lower.tri(cm, diag = TRUE)] <- NA

  # Determine combinations over cutoff
  combsAboveCutoff <- which(abs(cm) > cutoff)
  colsToCheck <- ceiling(combsAboveCutoff/nrow(cm))
  rowsToCheck <- combsAboveCutoff%%nrow(cm)

  # Exclude columns with variables over average correlation
  colsToDiscard <- averageCorr[colsToCheck] > averageCorr[rowsToCheck]
  rowsToDiscard <- !colsToDiscard

  # Get columns to discard
  deletecol <- c(colsToCheck[colsToDiscard], rowsToCheck[rowsToDiscard])
  deletecol <- unique(deletecol)

  # Which variables to discard
  o <- names(env)[deletecol]
  if(length(singular_var)>0) o <- unique( c(o,  names(singular_var) ) )
  o
}

#' Apply the adaptive best subset selection framework on a set of predictors
#'
#' @description
#' This is a wrapper function to fit the adaptive subset selection procedure outlined
#' in Zhu et al. (2021) and Zhu et al. (2020).
#' @inheritParams predictor_filter
#'
#' @param family A [`character`] indicating the family the observational data originates from.
#' @param tune.type [`character`] indicating the type used for subset evaluation.
#' Options are \code{c("gic", "ebic", "bic", "aic", "cv")} as listed in [abess].
#' @param lambda A [`numeric`] single lambda value for regularized best subset selection (Default: \code{0}).
#' @param weight Observation weights. When weight = \code{NULL}, we set weight = \code{1} for each observation as default.
#' @references
#' * abess: A Fast Best Subset Selection Library in Python and R. Jin Zhu, Liyuan Hu, Junhao Huang, Kangkang Jiang, Yanhang Zhang, Shiyun Lin, Junxian Zhu, Xueqin Wang (2021). arXiv preprint arXiv:2110.09697.
#' * A polynomial algorithm for best-subset selection problem. Junxian Zhu, Canhong Wen, Jin Zhu, Heping Zhang, Xueqin Wang. Proceedings of the National Academy of Sciences Dec 2020, 117 (52) 33117-33123; doi: 10.1073/pnas.2014241117
#' @keywords utils, internal
#' @returns A [`vector`] of variable names to exclude
predictors_filter_abess <- function( env, observed, method, family, tune.type = "cv", lambda = 0,
                                       weight = NULL, keep = NULL, ...){
  # Security checks
  assertthat::assert_that(is.data.frame(env) || is.matrix(env),
                          is.vector(observed),
                          is.numeric(lambda),
                          is.character(tune.type),
                          is.null(keep) || is.vector(keep),
                          is.null(weight) || is.vector(weight)
  )
  assertthat::assert_that(
    length(observed) == nrow(env), msg = "Number of observation unequal to number of covariate rows."
  )
  # Match family and type
  family <- match.arg(family, c("gaussian", "binomial", "poisson", "cox", "mgaussian", "multinomial",
                                "gamma"), several.ok = FALSE)
  tune.type <- match.arg(tune.type, c("gic", "ebic", "bic", "aic", "cv"), several.ok = FALSE)

  # Check that abess package is available
  check_package("abess")
  if(!isNamespaceLoaded("abess")) { attachNamespace("abess");requireNamespace('abess') }

  # Build model
  abess_fit <- abess::abess(x = env,
                            y = observed,
                            family = family,
                            tune.type = tune.type,
                            weight = weight,
                            lambda = lambda,
                            always.include = keep,
                            nfolds = 100, # Increase from default 5
                            num.threads = 0
  )

  if(anyNA(stats::coef(abess_fit)[,1]) ) {
    # Refit with minimum support size
    abess_fit <- abess::abess(x = env,
                              y = observed,
                              family = family,
                              lambda = lambda,
                              tune.type = tune.type,
                              weight = weight,
                              always.include = keep,
                              nfolds = 100, # Increase from default 5
                              # Minimum support site of 10% of number of covariates
                              support.size = ceiling(ncol(env) * 0.1),
                              num.threads = 0
    )

  }
  # Get best vars
  co <- stats::coef(abess_fit, support.size = abess_fit[["best.size"]])
  co <- names( which(co[,1] != 0))
  co <- co[grep("Intercept", co, ignore.case = TRUE, invert = TRUE)]
  # Make some checks on the list of reduced variables
  if(length(co) <= 2) {
    warning("Abess was likely to rigours. Likely to low signal-to-noise ratio.")
    return(NULL)
  } else {
    co
  }
}

#' All relevant feature selection using Boruta
#'
#' @description
#' This function uses the [Boruta] package to identify predictor variables with little information content. It iteratively
#' compares importances of attributes with importances of shadow attributes, created by shuffling original ones.
#' Attributes that have significantly worst importance than shadow ones are being consecutively dropped.
#'
#' @note
#' This package depends on the [ranger] package to iteratively fit randomForest models.
#'
#' @inheritParams predictor_filter
#' @param iter [`numeric`] on the number of maximal runs (Default: \code{100}). Increase if too many tentative left.
#' @param verbose [`logical`] whether to be chatty.
#' @references
#' * Miron B. Kursa, Witold R. Rudnicki (2010). Feature Selection with the Boruta Package. Journal of Statistical Software, 36(11), 1-13. URL https://doi.org/10.18637/jss.v036.i11.
#' @keywords utils, internal
#' @returns A [`vector`] of variable names to exclude.
predictors_filter_boruta <- function( env, observed, method, keep = NULL,
                                      iter = 100, verbose = getOption('ibis.setupmessages'), ...){
  # Security checks
  assertthat::assert_that(is.data.frame(env) || is.matrix(env),
                          is.null(observed) || is.vector(observed),
                          is.null(keep) || is.vector(keep),
                          is.numeric(iter), iter>10,
                          is.logical(verbose)
  )
  check_package("Boruta")

  # Get all variable names to test
  vars <- names(env)

  # Remove kept variables
  if(!is.null(keep)){
    keep <- keep[keep %in% vars] # Remove those not in the data.frame. For instance if a spatial effect is selected
    if(!is.null(keep) || length(keep) == 0) {
      env <- env |> dplyr::select(-keep)
      vars <- names(env)
    }
  }

  # Check for other common variables unlikely to be important
  if("Intercept" %in% vars) vars <- vars[-which(vars == "Intercept")]
  if("intercept" %in% vars) vars <- vars[-which(vars == "intercept")]
  if("ID" %in% vars) vars <- vars[-which(vars == "ID")]

  # Apply boruta
  bo_test <- Boruta::Boruta(env, y = observed,
                            maxRuns = iter,
                            # Verbosity
                            doTrace = ifelse(verbose, 1, 0))

  # Get from the bo_test object all variables that are clearly rejected
  res <- bo_test$finalDecision
  co <- names(res)[which(res == "Rejected")]
  if(length(co)==0) co <- NULL
  return(co)
}
