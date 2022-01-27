#' Validation of distribution object
#'
#' @description This function conducts a comprehensive model evaluation based on
#' either on the fitted point data or any supplied independent.
#' **Currently only supporting point datasets. For validation of integrated models more work is needed.**
#' @param mod A fitted [`BiodiversityDistribution`] object with set predictors.
#' @param method Should the validation be conducted on continuous metrics or thresholded? See Details.
#' @param layer In case multiple layers exist, which one to use? (Default: \code{'mean'}).
#' @param point A [`sf`] object with type `POINT` or `MULTIPOINT`.
#' @param point_column A [`character`] vector with the name of the column containing the independent observations.
#' (Default: \code{'observed'}).
#' @param ... Other parameters that are passed on. Currently unused.
#' @returns Return a tidy [`tibble`] with validation results.
#' @details The \code{validate} function does not work for all datasets equally.
#' @note If you use the Boyce Index, cite the original Hirzel et al. (2006) paper.
#'
#' @references
#' * Liu, C., White, M., Newell, G., 2013. Selecting thresholds for the prediction of species occurrence with presence-only data. J. Biogeogr. 40, 778–789. https://doi.org/10.1111/jbi.12058
#' * Hirzel, A. H., Le Lay, G., Helfer, V., Randin, C., & Guisan, A. (2006). Evaluating the ability of habitat suitability models to predict species presences. Ecological modelling, 199(2), 142-152.
#' @examples
#' \dontrun{
#'  # Assuming that mod is a distribution object and has a thresholded layer
#'  mod <- threshold(mod, method = "TSS")
#'  validate(mod, method = "discrete")
#'  }
#' @name validate
#' @aliases validate
#' @keywords train
#' @exportMethod validate
#' @export
NULL
methods::setGeneric("validate",
                    signature = methods::signature("mod"),
                    function(mod, method = 'continuous', layer = "mean",
                             point = NULL, point_column = 'observed', ...) standardGeneric("validate"))

#' @name validate
#' @rdname validate
#' @usage \S4method{validate}{ANY, character, sf, character, character}(mod, method, point, layer, point_column)
methods::setMethod(
  "validate",
  methods::signature(mod = "ANY"),
  function(mod, method = 'continuous', layer = "mean",
           point = NULL, point_column = 'observed', ...){
    assertthat::assert_that(
      inherits(mod, "DistributionModel"),
      "prediction" %in% mod$show_rasters(),
      inherits(point, 'sf') || is.null(point),
      is.null(point_column) || is.character(point_column),
      is.character(layer),
      is.character(method)
    )
    # Check that independent data is provided and if so that the used column is there
    if(!is.null(point)){
      assertthat::assert_that(is.character(point_column),
                              hasName(point, point_column),
                              anyNA(point[[point_column]])==FALSE
                              )
    }
    # Match method to be sure
    method <- match.arg(method, c('continuous', 'discrete'), several.ok = FALSE)

    # Get prediction and threshold if available
    prediction <- mod$get_data('prediction')[[layer]]
    if( any(grep('threshold', mod$show_rasters())) ){
      tr_lyr <- grep('threshold', mod$show_rasters(),value = TRUE)
      if(length(tr_lyr)>1) warning("There appear to be multiple thresholds. Using the first one.")
      threshold <- mod$get_data(tr_lyr[1])
      # Get mean layer if there are multiple
      if( grep(layer, names(threshold),value = TRUE ) != "") threshold <- threshold[[grep(layer, names(threshold),value = TRUE )]]
    } else { threshold <- NULL }

    # Check that threshold and method match
    if(is.null(threshold) && method == 'discrete'){
      if(getOption('ibis.setupmessages')) myLog('[Validation]','red','No threshold data found. Switching to continuous validation metrics.')
      method <- 'continuous'
    }

    # Get/check point data
    if(!is.null(point)){
      assertthat::assert_that(
        unique(sf::st_geometry_type(point)) %in% c('POINT', 'MULTIPOINT'),
        # Check that the point data has presence-absence information
        hasName(point, point_column),
        length( unique( point[[point_column]] ) ) >1
      )
      # If sf is different, reproject to prediction
      if(sf::st_crs(point)!= sf::st_crs(prediction)){
        point <- sf::st_transform(point, crs = sf::st_crs(prediction) )
      }
      if(!hasName(point, "name")) point$name <- "Validation data" # Assign a name for validation. Assuming only one dataset is present
    } else {
      # TODO: Think about how to do validation with non-point data
      # Get all point datasets and combine them
      point <- do.call(sf:::rbind.sf,
                        lapply(mod$model$biodiversity, function(y){
                         o <-  guess_sf(y$observations)
                         o$name <- y$name; o$type <- y$type
                         subset(o, select = c(point_column, "name", "type", "geometry"))
                        } )
                        ) %>% tibble::remove_rownames()
      # Add ID
    }
    assertthat::assert_that(nrow(point)>0,
                            hasName(point, point_column))
    # --- #
    # Do the extraction
    df <- as.data.frame(point)
    df$pred <- raster::extract(prediction, point)
    if(!is.null(threshold)) df$pred_tr <- raster::extract(threshold, point)
    # Remove any sfc column if present
    if(!is.null(attr(df, "sf_column"))) df[[attr(df, "sf_column")]] <- NULL
    # Remove any NAs
    df <- subset(df, complete.cases(df))
    if(nrow(df) < 2) stop("Validation was not possible owing to missing data.")
    # --- #
    # Messenger
    if(getOption('ibis.setupmessages')) myLog('[Validation]','green','Calculating validation statistics')

    # Output container
    results <- data.frame()

    for(dataset in unique(df$name)){
      # Subset to name
      df2 <- subset.data.frame(df, name == dataset)

      if(method == 'continuous'){
        # continuous evaluation
        assertthat::assert_that(hasName(df2, 'pred'),
                                hasName(df2, point_column)
        )
        #### Calculating Boyce index as in Hirzel et al. 2006
        # fit: A vector or Raster-Layer containing the predicted suitability values
        # obs: A vector containing the predicted suitability values or xy-coordinates (if fit is a Raster-Layer) of the validation points (presence records)
        # nclass : number of classes or vector with classes threshold. If nclass=0, Boyce index is calculated with a moving window (see next parameters)
        # windows.w : width of the moving window (by default 1/10 of the suitability range)
        # res : resolution of the moving window (by default 101 focals)
        # PEplot : if True, plot the predicted to expected ratio along the suitability class
        ecospat.boyce <-
          function(fit,
                   obs,
                   nclass = 0,
                   window.w = "default",
                   res = 100,
                   PEplot = TRUE){
            boycei <- function(interval, obs, fit) {
              fit.bin <- fit
              obs.bin <- obs
              fit.bin[fit[] >= interval[1] & fit[] <= interval[2]] <- "i"
              fit.bin[fit.bin != "i"] <- 0
              obs.bin[obs[] >= interval[1] & obs[] <= interval[2]] <- "i"
              obs.bin[obs.bin != "i"] <- 0
              pi <- length(which(obs.bin == "i")) / length(obs)
              ei <- length(which(fit.bin == "i")) / length(fit.bin)
              fi <- pi / ei
              return(fi)
            }

            if (window.w == "default") {
              window.w <- (max(fit, na.rm = TRUE) - min(fit, na.rm = TRUE)) / 10
            }

            interval <- c(min(fit, na.rm = TRUE), max(fit, na.rm = TRUE))
            mini <- interval[1]
            maxi <- interval[2]

            if (nclass == 0) {
              vec.mov <-
                seq(
                  from = mini,
                  to = maxi - window.w,
                  by = (maxi - mini - window.w) / res
                )

              vec.mov[res + 1] <-
                vec.mov[res + 1] + 1  #Trick to avoid error with closed interval in R

              interval <- cbind(vec.mov, vec.mov + window.w)
            } else if (length(nclass) > 1) {
              vec.mov <- c(mini, nclass)
              interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
            } else if (nclass > 0 & length(nclass) < 2) {
              vec.mov <- seq(from = mini,
                             to = maxi,
                             by = (maxi - mini) / nclass)
            }

            f <- apply(interval, 1, boycei, obs, fit)
            to.keep <- which(f != "NaN")  # index to keep no NaN data
            f <- f[to.keep]

            if (length(f) < 2) {
              b <- NA  #at least two points are necessary to draw a correlation
            } else {
              r <-
                c(1:length(f))[f != c(f[-1], FALSE)]  #index to remove successive duplicates
              b <-
                stats::cor(f[r], vec.mov[to.keep][r], method = "spearman")  # calculation of the spearman correlation (i.e. Boyce index) after removing successive duplicated values
            }

            HS <- apply(interval, 1, sum) / 2  # mean habitat suitability in the moving window
            HS[length(HS)] <- HS[length(HS)] - 1  #Correction of the 'trick' to deal with closed interval
            HS <- HS[to.keep] # exlude the NaN

            if (PEplot == TRUE) {
              plot(
                HS,
                f,
                xlab = "Habitat suitability",
                ylab = "Predicted/Expected ratio",
                col = "grey",
                cex = 0.75
              )
              points(HS[r], f[r], pch = 19, cex = 0.75)

            }

            results <- list(F.ratio = f,
                            Spearman.cor = round(b, 3),
                            HS = HS)
            return(results)
          }

        # Function for Root-mean square error
        RMSE <- function(pred, obs, na.rm = TRUE) {
          sqrt(mean((pred - obs)^2, na.rm = na.rm))
        }
        # Mean absolute error
        MAE <- function(pred, obs, na.rm = TRUE) {
          mean(abs(pred - obs), na.rm = na.rm)
        }
        # Function for log loss/cross-entropy loss.
        Poisson_LogLoss <- function(y_pred, y_true) {
          eps <- 1e-15
          y_pred <- pmax(y_pred, eps)
          Poisson_LogLoss <- mean(log(gamma(y_true + 1)) + y_pred - log(y_pred) * y_true)
          return(Poisson_LogLoss)
        }
        # Normalized Gini Coefficient
        NormalizedGini <- function(y_pred, y_true) {
          SumGini <- function(y_pred, y_true) {
            y_true_sort <- y_true[order(y_pred, decreasing = TRUE)]
            y_random <- 1:length(y_pred) / length(y_pred)
            y_Lorentz <- cumsum(y_true_sort) / sum(y_true_sort)
            SumGini <- sum(y_Lorentz - y_random)
            return(SumGini)
          }
          NormalizedGini <- SumGini(y_pred, y_true) / SumGini(y_true, y_true)
          return(NormalizedGini)
        }
        # Create output container
        out <- data.frame(
          modelid = as.character(mod$id),
          name = dataset,
          method = method,
          metric = c('n','rmse', 'mae',
                     'logloss','normgini',
                     'cont.boyce'),
          value = NA
        )
        # - #
        out$value[out$metric=='n'] <- nrow(df2) # Number of records
        out$value[out$metric=='rmse'] <- RMSE(pred = df2$pred, obs = df2[[point_column]]) # RMSE
        out$value[out$metric=='mae'] <- MAE(pred = df2$pred, obs = df2[[point_column]]) # Mean absolute error
        out$value[out$metric=='logloss'] <- Poisson_LogLoss(y_pred = df2$pred, y_true = df2[[point_column]])
        out$value[out$metric=='normgini'] <- NormalizedGini(y_pred = df2$pred, y_true = df2[[point_column]])
        # Boyce index. Wrap in try since is known to crash
        try({
          boi <- ecospat.boyce(obs = df2[[point_column]], fit = df2$pred, PEplot = FALSE)
        },silent = TRUE)
        if(exists('boi')) out$value[out$metric=='cont.boyce'] <- boi$Spearman.cor

        # Append to results
        results <- rbind(results, out);rm(df2, out)
      } else {
        # discrete evaluation
        assertthat::assert_that(hasName(df2, 'pred_tr'),
                                length(unique(df2[[point_column]])) > 1,
                                msg = "It appears as either the observed data or the threshold does not allow discrete validation.")
        # For discrete functions to work correctly, ensure that all values are 0/1
        df2[[point_column]] <- ifelse(df2[[point_column]] > 0, 1, 0 )
        # Build the confusion matrix
        ta  <-  sum((df2["pred_tr"] == 0) & (df2[point_column] == 0))
        fp  <-  sum((df2["pred_tr"] == 1) & (df2[point_column] == 0))
        fa  <-  sum((df2["pred_tr"] == 0) & (df2[point_column] == 1))
        tp  <-  sum((df2["pred_tr"] == 1) & (df2[point_column] == 1))

        # Output data.frame
        out <- data.frame(
          modelid =  as.character(mod$id),
          name = dataset,
          method = method,
          metric = c('n','auc','overall.accuracy', 'true.presence.ratio',
                     'precision','sensitivity', 'specificity',
                     'tss', 'f1', 'logloss',
                     'expected.accuracy', 'kappa'),
          value = NA
        )

        # Accuracy indices
        out$value[out$metric=='n'] <- N <- ta + fp + fa + tp # Total number of records
        out$value[out$metric=='overall.accuracy'] <- OA <- (tp + ta) / N # Overall accuracy
        out$value[out$metric=='true.presence.ratio'] <- FOM <- tp / (tp + fp + fa) # True presence classifications
        out$value[out$metric=='precision'] <- precision <- tp / (tp + fp) # Precision
        out$value[out$metric=='sensitivity'] <- Sensitivity <- tp / (tp + fa) # Sensitivity
        out$value[out$metric=='specificity'] <- Specificity <- ta / (ta + fp) # Specificity
        out$value[out$metric=='tss'] <- TSS <- Sensitivity + Specificity - 1 # True Skill statistic
        out$value[out$metric=='f1'] <- 2 * (precision * Sensitivity) / (precision + Sensitivity) # F1 score
        Prob_1and1 <- ((tp + fp) / N) * ((tp + fa) / N) # Probability presence
        Prob_0and0 <- ((ta + fa) / N) * ((ta + fp) / N) # Probability absence
        out$value[out$metric=='expected.accuracy'] <- Expected_accuracy <- Prob_1and1 + Prob_0and0 # Expected accuracy
        out$value[out$metric=='kappa'] <- (OA - Expected_accuracy) / (1 - Expected_accuracy)

        if("modEvA" %in% installed.packages()[,1]){
          # Calculate AUC
          out$value[out$metric=='auc'] <- modEvA::AUC(obs = df2[[point_column]], pred = df2[['pred_tr']], simplif = TRUE, plot = FALSE)
        }

        # Evaluate Log loss / Cross-Entropy Loss for a predicted probability measure
        # FIXME: Hacky. This likely won't work with specific formulations
        if(mod$model$biodiversity[[which(sapply(mod$model$biodiversity, function(x) x$name) == dataset)]]$family == 'binomial'){
          LogLoss <- function(y_pred, y_true) {
            eps <- 1e-15
            y_pred <- pmax(pmin(y_pred, 1 - eps), eps)
            LogLoss <- -mean(y_true * log(y_pred) + (1 - y_true) * log(1 - y_pred))
            return(LogLoss)
          }
          out$value[out$metric=='logloss'] <- LogLoss(y_pred = df2$pred_tr, y_true = df2[[point_column]])
        }

        results <- rbind(results, out); rm(df2, out)
      } # End of discrete clause

    }
    # Return result
    if(exists("results")) return(results) else warning('Something went wrong during the validation...')
  }
)
