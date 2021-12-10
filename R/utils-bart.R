#' Variable importance for dbarts models
#'
#' Variable importance measured in the proportion of total branches used for a given variable. Explicitly dropped variables are included as 0
#' @param model A fitted [dbarts] model
#' @concept Taken from embarcadero package
#' @return A [`data.frame`] with the variable importance information
#' @keywords utils, internals
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
#' @param model A fitted [dbarts::bart] model
#' @param envs A [`raster`] stack of predictors used in the model
#' @param x.vars The predictor variables to be mapped (Default: All)
#' @param equal Whether equal spacing on x breaks or quantiles is applied (Default: FALSE)
#' @param smooth Smoothing factor for the x breaks (works like partials). (Default: 1)
#' @param trace Add traces to plot (Default: TRUE)
#' @param ci Should credible intervals be shown (Default: TRUE)
#' @param ciwidth Width of credible intervals (numeric)
#' @param transform Backtransform using pnorm or not. Set to FALSE if response was not binomial
#' @concept Taken and adapted from [embarcadero] package
#' @references  Carlson, CJ. embarcadero: Species distribution modelling with Bayesian additive regression trees in r. Methods Ecol Evol. 2020; 11: 850– 858. https://doi.org/10.1111/2041-210X.13389
#' @return A [`Raster`] layer containing the partial effect
#' @keywords utils
#' @noRd
bart_partial_effect <- function (model, x.vars = NULL, equal = TRUE, smooth = 1, ci = TRUE,
                     ciwidth = 0.95, trace = TRUE, transform = TRUE) {

  assertthat::assert_that(
    inherits(model,'bart'),
    is.null(x.vars) || is.character(x.vars),
    # is.logical(equal), is.numeric(smooth),
    is.logical(ci) && is.numeric(ciwidth),
    is.logical(transform)
  )

  # Get Fit object
  if (class(model) == "bart") {
    fitobj <- model$fit
  }
  # If no x.vars are specified, use all
  if (is.null(x.vars)) raw <- fitobj$data@x else raw <- fitobj$data@x[, x.vars]

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
    pd <- dbarts::pdbart(model, xind = x.vars, levs = lev, pl = FALSE)

  } else {
    levq = c(0.5 - ciwidth/2, seq(0.1, 0.9, 0.1/smooth),
             0.5 + ciwidth/2)
    pd <- pdbart(model, xind = x.vars, levquants = levq,
                 pl = FALSE)
  }

  plots <- list()
  for (i in 1:length(pd$fd)) {
    if (length(unique(pd$fit$data@x[, pd$xlbs[[i]]])) == 2) {
      dfbin <- data.frame(pd$fd[[i]])
      colnames(dfbin) <- c(0, 1)
      dfbin <- reshape2::melt(dfbin)
      if (transform == TRUE) dfbin$value <- pnorm(dfbin$value)
      if (ci == FALSE) {
        g <- ggplot(dfbin, aes(x = variable, y = value)) +
          geom_boxplot() + labs(title = pd$xlbs[[i]],
                                y = "Response", x = "") + theme_light(base_size = 20) +
          theme(plot.title = element_text(hjust = 0.5),
                axis.title.y = element_text(vjust = 1.7))
      }
      else {
        g <- ggplot(dfbin, aes(x = variable, y = value)) +
          geom_boxplot(fill = "deepskyblue1") +
          labs(title = pd$xlbs[[i]], y = "Response",
               x = "") + theme_light(base_size = 20) +
          theme(plot.title = element_text(hjust = 0.5),
                axis.title.y = element_text(vjust = 1.7))
      }
      if (panels == FALSE) {
        g <- g + theme(plot.margin = unit(c(0.5, 0.5,
                                            0.5, 0.5), "cm"))
      }
      else {
        g <- g + theme(plot.margin = unit(c(0.15, 0.15,
                                            0.15, 0.15), "cm"))
      }
      plots[[i]] <- g
    } else {
      q50 <- apply(pd$fd[[i]], 2, median)
      if (transform == TRUE) q50 <- pnorm(q50)
      df <- data.frame(x = pd$levs[[i]], med = q50)
      if(ci) {
        # Lower bound
        q05 <- apply(pd$fd[[i]], 2, quantile, probs = 0.5 -
                       ciwidth/2)
        if (transform == TRUE) q05 <- pnorm(q05)
        # Upper bound
        q95 <- apply(pd$fd[[i]], 2, quantile, probs = 0.5 +
                       ciwidth/2)
        if(transform) q95 <- pnorm(q95)

        df$q05 <- q05;df$q95 <- q95

      }
      if(trace) {
        f <- data.frame(t(pd$fd[[i]]))
        df <- cbind(df, f)
      }
      g <- ggplot(df, aes(x = x, y = med)) + labs(title = pd$xlbs[[i]],
                                                  y = "Response", x = "") + theme_light(base_size = 20) +
        theme(plot.title = element_text(hjust = 0.5),
              axis.title.y = element_text(vjust = 1.7))
      if(ci) {
        alpha2 <- 0.05
        k <- 4
      } else {
        alpha2 <- 0.025 * (fitobj$control@n.trees/200)
        k <- 2
      }
      if(trace) {
        if(transform) {
          for (j in 1:nrow(pd$fd[[i]])) {
            g <- g + geom_line(aes_string(y = pnorm(df[,
                                                       j + k])), alpha = alpha2)
          }
        } else {
          for (j in 1:nrow(pd$fd[[i]])) {
            g <- g + geom_line(aes_string(y = df[, j +
                                                   k]), alpha = alpha2)
          }
        }
      }
      # Add credible interval
      if(ci) {
        g <- g + geom_ribbon(aes(ymin = q05, ymax = q95),
                             fill = "deepskyblue1", alpha = 0.3)
      }
      g <- g + geom_line(size = 1.25)
      if(panels == FALSE) {
        g <- g + theme(plot.margin = unit(c(0.5, 0.5,
                                            0.5, 0.5), "cm"))
      }
      else {
        g <- g + theme(plot.margin = unit(c(0.15, 0.15,
                                            0.15, 0.15), "cm"))
      }
      plots[[i]] <- g
    }
  }
  if (panels == TRUE) {
    return(gridExtra::grid.arrange(plotlist = plots))
  } else {
    return(plots)
  }
}

#' Spatial partial effects for bart models adapted from embarcadero package
#'
#' @param model A fitted [dbarts::bart] model
#' @param envs A [`raster`] stack of predictors used in the mode]l
#' @param x.vars The predictor variables to be mapped (Default: All)
#' @param equal Whether equal spacing on x breaks or quantiles is applied (Default: FALSE)
#' @param smooth Smoothing factor for the x breaks (works like partials). (Default: 1)
#' @param transform Backtransform using pnorm or not. Set to FALSE if response was not binomial
#' @concept Taken and adapted from embarcadero package]
#' @references  Carlson, CJ. embarcadero: Species distribution modelling with Bayesian additive regression trees in r. Methods Ecol Evol. 2020; 11: 850– 858. https://doi.org/10.1111/2041-210X.13389
#' @return A [`raster`] layer containing the partial effect
#' @keywords utils
#' @noRd
bart_partial_space <- function(model, envs, x.vars = NULL, equal = FALSE, smooth = 1, transform = TRUE){
  # Input checks
  assertthat::assert_that(
    inherits(model,'bart'),
    inherits(envs, 'Raster'),
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
        dfbin$value <- pnorm(dfbin$value)
      }
      # FIXME: To replace with base::aggregate to get rid of dplyr dependency
      dfbin <- dfbin %>% group_by(variable) %>% summarize(value = median(value)) %>%
        data.frame()
      colnames(dfbin) <- c("is", "becomes")
      dfbin$is <- as.numeric(as.character(dfbin$is))
      if (class(envs) %in% c("RasterStack", "RasterBrick")) {
        lyrtmp <- envs[[pd$xlbs[[i]]]]
        lyrtr <- raster::reclassify(lyrtmp, as.matrix(dfbin))
      } else if (class(envs) == "list") {
        lyrtr <- lapply(envs, function(x) {
        lyrtmp <- x[[pd$xlbs[[i]]]]
          return(raster::reclassify(lyrtmp, as.matrix(dfbin)))
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
      q50 <- pnorm(apply(pd$fd[[i]], 2, median))
      if (transform == TRUE) { q50 <- pnorm(q50) }
      df <- data.frame(x = pd$levs[[i]], med = q50)
      nmax <- length(df$x)
      xmeds <- (df$x[2:nmax] - df$x[1:(nmax - 1)])/2 + df$x[1:(nmax - 1)]

      if (class(envs) %in% c("RasterStack", "RasterBrick")) {
        lyrtmp <- envs[[pd$xlbs[[i]]]]
        xmat <- data.frame(from = c(min(cellStats(lyrtmp,
                                                  min), min(df$x)), xmeds), to = c(xmeds, max(cellStats(lyrtmp,
                                                                                                        max), max(df$x))), becomes = df$med)
        lyrtr <- raster::reclassify(lyrtmp, xmat, include.lowest = TRUE)
      } else if (class(envs) == "list") {
        lyrtr <- lapply(envs, function(x) {
          lyrtmp <- x[[pd$xlbs[[i]]]]
          xmat <- data.frame(from = c(min(cellStats(lyrtmp,
                                                    min), min(df$x)), xmeds), to = c(xmeds, max(cellStats(lyrtmp,
                                                                                                          max), max(df$x))), becomes = df$med)
          return(raster::reclassify(lyrtmp, xmat, include.lowest = TRUE))
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
    return(stack(pdstack))
  }
}
