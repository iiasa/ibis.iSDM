#' Built formula for BART model
#'
#' @description
#' This function built a formula for a `engine_bart()` model.
#' @param obj A [`list()`] object containing the prepared model data for a given biodiversity dataset.
#' @author Martin Jung
#' @note Function is not meant to be run outside the train() call.
#' @keywords internal
#' @noRd
built_formula_bart <- function(obj){
  assertthat::assert_that(
    is.list(obj),
    length(obj) > 0,
    assertthat::has_name(obj, "observations"),
    assertthat::has_name(obj, "equation"),
    assertthat::has_name(obj, "predictors_names"),
    msg = "Error in model object. This function is not meant to be called outside ouf train()."
  )

  # Default equation found
  if(is.Waiver(obj$equation) || obj$equation == '<Default>'){
    # Construct formula with all variables
    form <- paste( 'observed' ,
                   ifelse(obj$family=='poisson', '/w', ''), '~ ',
                   paste(obj$predictors_names, collapse = " + "))
    # Convert to formula
    form <- to_formula(form)
  } else {
    if(getOption('ibis.setupmessages')) myLog('[Estimation]','yellow','Use custom model equation')
    form <- to_formula(obj$equation)
    # If response is missing, add manually
    if(attr(stats::terms(form), "response")==0){
      if(obj$type == "poipo"){
        form <- stats::update.formula(form, "observed/w ~ .")
      } else {
        form <- stats::update.formula(form, "observed ~ .")
      }
    }
    assertthat::assert_that(
      is.formula(form),
      attr(stats::terms(form), "response")==1, # Has Response
      all( all.vars(form) %in% c('observed','w', model[['predictors_names']]) )
    )
  }
  return(form)
}

#' Variable importance for dbarts models
#'
#' @description
#' Variable importance measured in the proportion of total branches used for a given variable.
#' Explicitly dropped variables are included as \code{0}.
#' @param model A fitted [dbarts] model.
#' @concept Taken from the \pkg{embarcadero} package.
#' @return A [`data.frame`] with the variable importance information.
#' @keywords utils, internal
#' @noRd
varimp.bart <- function(model){
  assertthat::assert_that(class(model) == 'bart',
                          ("fit" %in% names(model)),
                          msg = 'Model not correctly specified or keeptrees set to FALSE.' )

  # Get basenames and summarize variable counts in trees
  basenames <- unlist(attr(model$fit$data@x, "drop"))
  names <- names(which(basenames == FALSE))
  varimp <- colMeans(model$varcount/rowSums(model$varcount))
  fitobj <- model$fit

  var.df <- data.frame(names, varimp)

  missing <- attr(fitobj$data@x, "term.labels")[!(attr(fitobj$data@x,
                                                       "term.labels") %in% names(unlist(attr(fitobj$data@x,
                                                                                               "drop"))))]
  # If any variables were removed, still add them to the data.frame
  if (length(missing) > 0) {
    missing.df <- data.frame(names = missing, varimp = 0)
    var.df <- rbind(var.df, missing.df)
  }
  var.df <- var.df[order(var.df$varimp,decreasing = TRUE),]
  return(var.df)
}

#' Partial effects for bart models adapted from embarcadero package
#'
#' @param model A fitted [dbarts::bart] model.
#' @param envs A [`SpatRaster`] stack of predictors used in the model.
#' @param x.vars The predictor variables to be mapped (Default: \code{"All"}).
#' @param equal Whether equal spacing on x breaks or quantiles is applied (Default: \code{FALSE}).
#' @param smooth Smoothing factor for the x breaks (works like partials). (Default: \code{1}).
#' @param transform Backtransform using pnorm or not. Set to \code{FALSE} if response was not binomial.
#' @param values Either a [`numeric`] vector of supplied value ranges or \code{NULL} (Default).
#' @param plot Whether a model should be created (Default: \code{TRUE}).
#' @concept Function taken and adapted from the [embarcadero] package.
#' @references
#' * Carlson, CJ. embarcadero: Species distribution modelling with Bayesian additive regression trees in r. Methods Ecol Evol. 2020; 11: 850– 858. https://doi.org/10.1111/2041-210X.13389
#' @examples
#' \dontrun{
#' mod <- distribution(background) |>
#'           add_biodiversity(speciesdata) |>
#'           add_predictors(covariates) |>
#'           engine_bart()
#'
#' bart_partial_effect(mod, "bio01mean.tif")
#'
#' }
#' @return A [`SpatRaster`] layer containing the partial effect
#' @keywords utils
#' @noRd
bart_partial_effect <- function (model, x.vars = NULL, equal = FALSE,
                                 smooth = 1, transform = TRUE, values = NULL, plot = TRUE) {

  assertthat::assert_that(
    inherits(model,'bart'),
    is.null(x.vars) || is.character(x.vars),
    is.logical(transform),
    is.logical(equal),
    is.numeric(smooth),
    is.null(values) || is.numeric(values),
    is.logical(plot)
  )

  # Get Fit object
  if (inherits(model,"bart")) {
    fitobj <- model$fit
  }
  # If no x.vars are specified, use all
  if(!is.null(values)){
    raw <- list()
    raw[[x.vars]] <- values
    raw <- raw |> as.data.frame()
  } else {
    if (is.null(x.vars)) raw <- fitobj$data@x else raw <- fitobj$data@x[, x.vars]
  }

  # Define binning in equal area width or not
  if(equal) {
    if(!is.null(x.vars) && length(x.vars) == 1) {
      minmax <- data.frame(mins = min(raw), maxs = max(raw))
    }
    else {
      minmax <- data.frame(mins = apply(raw, 2, min), maxs = apply(raw,
                                                                   2, max))
    }
    lev <- lapply(c(1:nrow(minmax)), function(i) {
      seq(minmax$mins[i], minmax$maxs[i], (minmax$maxs[i] -
                                             minmax$mins[i])/(10 * smooth))
    })
    for (i in 1:length(lev)) {
      if (length(lev) == 1) {
        if (length(unique(raw)) == 2) {
          lev[[i]] <- unique(raw)
        }
      }
      else {
        if (length(unique(raw[, i])) == 2) {
          lev[[i]] <- unique(raw[, i])
        }
      }
    }
    pd <- dbarts::pdbart(model, xind = x.vars, levs = lev, keepevery = 10, pl = FALSE)

  } else {
    levq = c(0.05, seq(0.1, 0.9, 0.1/smooth),
             0.95)
    pd <- dbarts::pdbart(model, xind = x.vars, levquants = levq,
                         keepevery = 10,
                         pl = FALSE)
  }

  out <- data.frame()
  # Summarize the posterior
  for(i in 1:length(pd$fd)){
    df <- data.frame(partial_effect = pd$levs[[i]])
    # Summarize quantiles and sd from posterior
    ms <- as.data.frame(
      cbind( apply(pd$fd[[i]], 2, function(x) mean(x, na.rm = TRUE)),
             matrixStats::colSds(pd$fd[[i]]),
             matrixStats::colQuantiles(pd$fd[[i]], probs = c(.05,.5,.95)),
             apply(pd$fd[[i]], 2, mode)
      )
    )
    names(ms) <- c("mean","sd", "q05", "q50", "q95", "mode")
    if(transform) ms[,c("mean","q05","q50","q95","mode")] <- apply(ms[,c("mean","q05","q50","q95","mode")], 2, stats::pnorm)
    ms$cv <- ms$sd / ms$mean
    ms$variable <- pd$xlbs[[i]]
    df <- cbind(df, ms)
    out <- rbind(out, df);rm(df)
  }

  # Create plot if specified
  if(plot){
    # Construct overall plot. Might fail for multiple variables
    g <- ggplot2::ggplot(out, ggplot2::aes(x = partial_effect, y = q50)) +
      ggplot2::labs(title = "", y = expression(y^q50), x = unique(out$variable)) +
      ggplot2::theme_light(base_size = 20) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = q05, ymax = q95),
                           fill = "deepskyblue1", alpha = 0.3) +
      ggplot2::geom_line(size = 1.25) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     axis.title.y = ggplot2::element_text(vjust = 1.7))
    # If multiple variables, add facets
    if(length(unique(out$variable))>1){
      g <- g + ggplot2::facet_wrap(~variable)
    }
    g
  }
  # Return the partial results
  return(out)
}

#' Spatial partial effects for bart models adapted from embarcadero package
#'
#' @param model A fitted [dbarts::bart] model.
#' @param envs A [`SpatRaster`] stack of predictors used in the model.
#' @param x.vars The predictor variables to be mapped (Default: All).
#' @param equal Whether equal spacing on x breaks or quantiles is applied (Default: \code{FALSE}).
#' @param smooth Smoothing factor for the x breaks (works like partials). (Default: \code{1}).
#' @param transform Backtransform using pnorm or not. Set to FALSE if response was not Binomial.
#' @concept Taken and adapted from embarcadero package.
#' @references
#' * Carlson, CJ. embarcadero: Species distribution modelling with Bayesian additive regression trees in r. Methods Ecol Evol. 2020; 11: 850– 858. https://doi.org/10.1111/2041-210X.13389
#' @return A [`SpatRaster`] layer containing the partial effect.
#' @keywords utils
#' @noRd
bart_partial_space <- function(model, envs, x.vars = NULL, equal = FALSE, smooth = 1, transform = TRUE){
  # Input checks
  assertthat::assert_that(
    inherits(model,'bart'),
    is.Raster(envs),
    is.null(x.vars) || is.character(x.vars),
    is.logical(equal), is.numeric(smooth),
    is.logical(transform)
  )
  # No x.vars chosen, take all variables
  if (is.null(x.vars)) raw <- model$fit$data@x else raw <- model$fit$data@x[, x.vars]

  if (equal == TRUE) {
    if (!is.null(x.vars) && length(x.vars) == 1) {
      minmax <- data.frame(mins = min(raw), maxs = max(raw))
    }
    else {
      minmax <- data.frame(mins = apply(raw, 2, min), maxs = apply(raw, 2, max))
    }
    lev <- lapply(c(1:nrow(minmax)), function(i) {
      seq(minmax$mins[i], minmax$maxs[i], (minmax$maxs[i] -
                                             minmax$mins[i])/(10 * smooth))
    })
    for (i in 1:length(lev)) {
      if (length(lev) == 1) {
        if (length(unique(raw)) == 2) {
          lev[[i]] <- unique(raw)
        }
      }
      else {
        if (length(unique(raw[, i])) == 2) {
          lev[[i]] <- unique(raw[, i])
        }
      }
    }
  } else {
    lev = c(0.05, seq(0.1, 0.9, 0.1/smooth), 0.95)
  }
  # Use barts to get partial effects
  pd <- dbarts::pdbart(model, xind = x.vars, levquants = lev, pl = FALSE)
  # Loop through
  for(i in 1:length(pd$fd)) {
    # Get first rasterlayer class
    envi <- envs[[pd$xlbs[[i]]]]
    if (length(unique(pd$fit$data@x[, pd$xlbs[[i]]])) ==2) {
      print(paste("WARNING: ", " is a binary variable; the plot will look bad/be uninformative",
                  sep = pd$xlbs[[i]]))

      dfbin <- data.frame(pd$fd[[i]])
      colnames(dfbin) <- c(0, 1)
      dfbin <- reshape2::melt(dfbin)
      if (transform == TRUE) {
        dfbin$value <- stats::pnorm(dfbin$value)
      }
      # FIXME: To replace with base::aggregate to get rid of dplyr dependency
      dfbin <- dfbin |> dplyr::group_by(variable) |> dplyr::summarize(value = stats::median(value)) |>
        data.frame()
      colnames(dfbin) <- c("is", "becomes")
      dfbin$is <- as.numeric(as.character(dfbin$is))
      if (is.Raster(envs) && (terra::nlyr(envs)>1) ) {
        lyrtmp <- envs[[pd$xlbs[[i]]]]
        lyrtr <- terra::reclassify(lyrtmp, as.matrix(dfbin))
      } else if (inherits(envs, "list")) {
        lyrtr <- lapply(envs, function(x) {
        lyrtmp <- x[[pd$xlbs[[i]]]]
          return(terra::reclassify(lyrtmp, as.matrix(dfbin)))
        })
      }
      if (exists("pdstack")) {
        pdstack <- c(pdstack, lyrtr)
      }
      else {
        pdstack <- c(lyrtr)
      }
    } else {
      # Nothing binary, calculate median
      q50 <- stats::pnorm(apply(pd$fd[[i]], 2, median))
      if (transform == TRUE) { q50 <- stats::pnorm(q50) }
      df <- data.frame(x = pd$levs[[i]], med = q50)
      nmax <- length(df$x)
      xmeds <- (df$x[2:nmax] - df$x[1:(nmax - 1)])/2 + df$x[1:(nmax - 1)]

      if(is.Raster(envs) && terra::nlyr(envs)>1) {
        lyrtmp <- envs[[pd$xlbs[[i]]]]
        xmat <- data.frame(from = c(min( terra::global(lyrtmp, "min", na.rm = TRUE)[,1], min(df$x)), xmeds),
                           to = c(xmeds, max( terra::global(lyrtmp, "max", na.rm = TRUE)[,1], max(df$x))), becomes = df$med)
        lyrtr <- terra::reclassify(lyrtmp, xmat, include.lowest = TRUE)
      } else if (inherits(x = envs, what = "list")) {
        lyrtr <- lapply(envs, function(x) {
          lyrtmp <- x[[pd$xlbs[[i]]]]
          xmat <- data.frame(from = c(min(terra::global(lyrtmp, "min", na.rm = TRUE)[,1], min(df$x)), xmeds),
                             to = c(xmeds, max(terra::global(lyrtmp, "max", na.rm = TRUE)[,1], max(df$x))), becomes = df$med)
          return(terra::reclassify(lyrtmp, xmat, include.lowest = TRUE))
        })
      }
      # Check if stack exists, otherwise create
        if (exists("pdstack")) {
          pdstack <- c(pdstack, lyrtr)
        } else {
          pdstack <- c(lyrtr)
        }
    }
  }
  # Return the output
  if (exists("pdstack")) {
    return(pdstack)
  }
}
