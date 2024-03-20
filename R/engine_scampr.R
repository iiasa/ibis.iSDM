#' @include class-engine.R class-distributionmodel.R
NULL

#' Engine for process models using scampr
#'
#' @description
#' Similar to others, this engine enables the fitting and prediction of
#' log-Gaussian Cox process (LGCP) and Inhomogeneous Poisson process (IPP) processes.
#' It uses the \code{scampr} package, which uses maximum likelihood estimation
#' fitted via \code{TMB} (Template Model Builder).
#'
#' It also support the addition of spatial latent effects which can be added via
#' Gaussian fields and approximated by 'FRK' (Fixed Rank Kriging) and are
#' integrated out using either variational or Laplace approximation.
#'
#' The main use case for this engine is as an alternative to [engine_inlabru()] and
#' [engine_inla()] for fitting iSDMs, e.g. those combining both presence-only
#' and presence-absence point occurrence data.
#'
#' @details
#' This engine may only be used to predict for one or two datasets at most. It
#' supports only presence-only PPMs and presence/absence Binary GLMs, or 'IDM'
#' (for an integrated data model).
#'
#' @note
#' * The package can currently be installed from github directly only \code{"ElliotDovers/scampr"}
#' * Presence-absence models in SCAMPR currently only support cloglog link functions!
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param type The mode used for creating (posterior or prior) predictions. Either stting
#' \code{"link"} or \code{"response"} (Default: \code{"response"}).
#' @param dens A [`character`] on how predictions are made, either from the \code{"posterior"} (Default)
#' or \code{"prior"}.
#' @param maxit A [`numeric`] on the number of iterations for the optimizer (Default: \code{500}).
#' @param ... Other parameters passed on.
#'
#' @returns An [Engine].
#'
#' @references
#' * Dovers, E., Popovic, G. C., & Warton, D. I. (2024). A fast method for fitting integrated species distribution models. Methods in Ecology and Evolution, 15(1), 191-203.
#' * Dovers, E., Stoklosa, D., and Warton D. I. (2024). Fitting log-Gaussian Cox processes using generalized additive model software. The American Statistician, 1-17.
#' @family engine
#'
#' @examples
#' \dontrun{
#' # Load background
#' background <- terra::rast(system.file('extdata/europegrid_50km.tif',
#' package='ibis.iSDM',mustWork = TRUE))
#'
#' # Add GLM as an engine
#' x <- distribution(background) |> engine_scampr()
#' }
#' @name engine_scampr
NULL

#' @rdname engine_scampr
#' @export
engine_scampr <- function(x,
                       type = "response",
                       dens = "posterior",
                       maxit = 500,
                       ...) {

  # Check whether package is available
  check_package('scampr')
  if(!("scampr" %in% loadedNamespaces()) || ('scampr' %notin% utils::sessionInfo()$otherPkgs) ) {
    try({requireNamespace('scampr');attachNamespace("scampr")},silent = TRUE)
  }

  # assert that arguments are valid
  assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                          inherits(x$background,'sf'),
                          is.character(type),
                          is.character(dens),
                          is.numeric(maxit), maxit >=1
  )
  type <- match.arg(type, choices = c("predictor","link", "response"),several.ok = FALSE)
  if(type=="predictor") type <- "link" # Convenience conversion

  # Expectation form for prediction
  dens <- match.arg(dens, choices = c("posterior", "prior"), several.ok = FALSE)

  # Create a background raster
  if(is.Waiver(x$predictors)){
    # Create from background
    template <- terra::rast(
      ext = terra::ext(x$background),
      crs = terra::crs(x$background),
      res = c(diff( (sf::st_bbox(x$background)[c(1,3)]) ) / 100, # Simplified assumption for resolution
              diff( (sf::st_bbox(x$background)[c(1,3)]) ) / 100
      )
    )
  } else {
    # If predictor existing, use them
    template <- emptyraster(x$predictors$get_data() )
  }

  # Burn in the background
  template <- terra::rasterize(x$background, template, field = 0)

  # Set up the parameter list
  params <- list(
    type = type,
    dens = dens,
    maxit = maxit,
    ...
  )

  # Print a message in case there is already an engine object
  if(!is.Waiver(x$engine)) myLog('[Setup]','yellow','Replacing currently selected engine.')

  # Define new engine object of class
  eg <- Engine

  # Dummy function for spatial latent effects
  eg$set("public", "calc_latent_spatial", function(type = NULL, priors = NULL){
    new_waiver()
  }, overwrite = TRUE)

  # Dummy function for getting the equation of latent effects
  eg$set("public", "get_equation_latent_spatial", function(method){
    new_waiver()
  }, overwrite = TRUE)

  # Function to respecify the control parameters
  eg$set("public", "set_control", function(params){
    assertthat::assert_that(is.list(params))
    # Overwrite existing
    self$data$params <- params
    invisible()
  },overwrite = TRUE)

  #### Setup function ----
  eg$set("public", "setup", function(model, settings = NULL, ...){
    # Simple security checks
    assertthat::assert_that(
      assertthat::has_name(model, 'background'),
      assertthat::has_name(model, 'biodiversity'),
      inherits(settings,'Settings') || is.null(settings),
      nrow(model$predictors) == terra::ncell(self$get_data('template')),
      !is.Waiver(self$get_data("params")),
      length(model$biodiversity)>0 && length(model$biodiversity) <=2
    )
    # Messenger
    if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','green','Engine setup.')

    # Get parameters
    params <- self$data$params
    settings$set('type', params$type)

    for(j in 1:length(model$biodiversity)){

      # Distribution specific procedure
      fam <- model$biodiversity[[j]]$family
      form <- model$biodiversity[[j]]$equation

      if(any(model$biodiversity[[j]]$expect!=1)){
        if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','yellow','Custom weights for dataset are ignored in this engine!')
      }

      # If a poisson family is used get pseudo-abs
      if(fam == "poisson"){
        # Get background layer
        bg <- self$get_data("template")
        assertthat::assert_that(!is.na( terra::global(bg, "min", na.rm = TRUE)[,1]))

        # Add pseudo-absence points
        presabs <- add_pseudoabsence(df = model$biodiversity[[j]]$observations,
                                     field_occurrence = 'observed',
                                     template = bg,
                                     settings = model$biodiversity[[j]]$pseudoabsence_settings)
        if(inherits(presabs, 'sf')) presabs <- presabs |> sf::st_drop_geometry()

        # Sample environmental points for absence only points
        abs <- subset(presabs, observed == 0)
        # Re-extract environmental information for absence points
        envs <- get_rastervalue(coords = abs[,c('x','y')],
                                env = model$predictors_object$get_data(df = FALSE),
                                rm.na = FALSE)
        if(assertthat::has_name(model$biodiversity[[j]]$predictors, "Intercept")){ envs$Intercept <- 1}

        # Format out
        df <- rbind(model$biodiversity[[j]]$predictors[,c('x','y','Intercept', model$biodiversity[[j]]$predictors_names)],
                    envs[,c('x','y','Intercept', model$biodiversity[[j]]$predictors_names)] )
        any_missing <- which(apply(df, 1, function(x) any(is.na(x))))
        if(length(any_missing)>0) {
          presabs <- presabs[-any_missing,] # This works as they are in the same order
          model$biodiversity[[j]]$expect <- model$biodiversity[[j]]$expect[-any_missing]
        }
        df <- subset(df, stats::complete.cases(df))
        assertthat::assert_that(nrow(presabs) == nrow(df))

        # Check that expect matches
        if(length(model$biodiversity[[j]]$expect)!=nrow(df)){
          # Fill the absences with 1 as multiplier. This works since absences follow the presences
          model$biodiversity[[j]]$expect <- c( model$biodiversity[[j]]$expect,
                                               rep(1, nrow(presabs)-length(model$biodiversity[[j]]$expect) ))
        }

        # Assign quadrature weights name and coord name check
        w <- ppm_weights(df = presabs,
                         pa = presabs$observed,
                         bg = bg,
                         use_area = TRUE,
                         weight = 1e-6, # Arbitrary small weight
                         type = "DWPR" # Weights for down-weighted Poisson regression
        )
        presabs$scampr.quad.size <- w
        assertthat::assert_that(utils::hasName(presabs,"x"),
                                utils::hasName(presabs,"y"),
                                length(unique(presabs$scampr.quad.size))==2)

        # Overwrite observation data
        model$biodiversity[[j]]$observations <- presabs

        # Preprocessing security checks
        assertthat::assert_that( all( model$biodiversity[[j]]$observations[['observed']] >= 0 ),
                                 any(!is.na(presabs[['observed']])),
                                 nrow(df) == nrow(model$biodiversity[[j]]$observations)
        )

        # Add offset if existent
        if(!is.Waiver(model$offset)){
          ofs <- get_rastervalue(coords = df[,c('x','y')],
                                 env = model$offset_object,
                                 rm.na = FALSE)
          # Rename to spatial offset
          names(ofs)[which(names(ofs)==names(model$offset_object))] <- "spatial_offset"
          model$biodiversity[[j]]$offset <- ofs
        }

        model$biodiversity[[j]]$predictors <- df

      } else if(fam == "binomial"){
        # Check that observations are all <=1
        model$biodiversity[[j]]$observations[['observed']] <- ifelse(model$biodiversity[[j]]$observations[['observed']]>=1,1,0)
        assertthat::assert_that(utils::hasName(model$biodiversity[[j]]$observations,"x"),
                                utils::hasName(model$biodiversity[[j]]$observations,"y")
        )
      }

    }

    # Instead of invisible return the model object
    return( model )
  },overwrite = TRUE)

  # Training function
  eg$set("public", "train", function(model, settings, ...){
    assertthat::assert_that(
      inherits(settings,'Settings'),
      is.list(model),length(model)>1,
      # Check that model id and setting id are identical
      settings$modelid == model$id
    )

    # Messenger
    if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','green',paste0('Starting fitting'))

    # Verbosity
    verbose <- settings$get("verbose")

    # Set prediction type also for later
    settings$set('type', self$get_data("params")$type)

    # seed
    seed <- settings$get("seed")
    if(is.Waiver(seed)) { settings$set('seed', getOption("ibis.seed", default = 1000)) }

    # Get output raster
    prediction <- self$get_data('template')

    # Get parameters control
    params <- self$get_data('params')

    # All other needed data for model fitting
    if(length(model$biodiversity)>1){
      model.type <- "IDM"
    } else {
      model.type <- switch(model$biodiversity[[1]]$family,
        "poisson" = "PO",
        "binomial" = "PA"
      )
      if(getOption('ibis.setupmessages', default = TRUE) && model.type=="PA") myLog('[Estimation]','yellow',paste0('Currently only cloglog link functions are supported for poipa!'))
    }

    # Get formula and combine data dependent on type
    if(model.type == "PO"){
      df <- cbind(model$biodiversity[[1]]$predictors,
                  data.frame(observed = model$biodiversity[[1]]$observations[,'observed', drop = TRUE],
                             scampr.quad.size = model$biodiversity[[1]]$observations[,'scampr.quad.size', drop = TRUE])
      )
      df <- subset(df, select = c(model$biodiversity[[1]]$predictors_names, "observed", "scampr.quad.size","x","y"))
      # Get overall equation
      equation <- model$biodiversity[[1]]$equation
    } else if(model.type == "PA") {
      df <- cbind(model$biodiversity[[1]]$predictors,
                  data.frame(observed = model$biodiversity[[1]]$observations[,'observed', drop = TRUE])
      )
      df <- subset(df, select = c(model$biodiversity[[1]]$predictors_names, "observed", "x","y"))
      # Get overall equation
      equation <- model$biodiversity[[1]]$equation
    } else {
      # Integrated model
      if(model$biodiversity[[1]]$type=="poipo"){
        df.po <- cbind(model$biodiversity[[1]]$predictors,
                       data.frame(observed = model$biodiversity[[1]]$observations[,'observed', drop = TRUE],
                                  scampr.quad.size = model$biodiversity[[1]]$observations[,'scampr.quad.size', drop = TRUE])
        )
        df.pa <- cbind(model$biodiversity[[2]]$predictors,
                       data.frame(observed = model$biodiversity[[2]]$observations[,'observed', drop = TRUE])
        )
      } else {
        df.pa <- cbind(model$biodiversity[[1]]$predictors,
                       data.frame(observed = model$biodiversity[[1]]$observations[,'observed', drop = TRUE],
                                  scampr.quad.size = model$biodiversity[[1]]$observations[,'scampr.quad.size', drop = TRUE])
        )
        df.po <- cbind(model$biodiversity[[2]]$predictors,
                       data.frame(observed = model$biodiversity[[2]]$observations[,'observed', drop = TRUE])
        )
      }

      # Compare terms and raise warning otherwise
      te1 <- attr(stats::terms.formula(model$biodiversity[[1]]$equation), "term.labels")
      te2 <- attr(stats::terms.formula(model$biodiversity[[2]]$equation), "term.labels")
      if(!all(te1 %in% te2) || !(length(te1)==length(te2))){
        if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','red',paste0('Dataset-specific formulas are not yet supported. Combining both objects!'))
        equation <- combine_formulas(model$biodiversity[[1]]$equation,
                                     model$biodiversity[[2]]$equation)
      } else {
        equation <- model$biodiversity[[1]]$equation
      }

      assertthat::assert_that(
        is.formula(equation),
        nrow(df.po)>0,
        utils::hasName(df.po, "x"),utils::hasName(df.po, "y"),
        nrow(df.pa)>0,
        utils::hasName(df.pa, "x"),utils::hasName(df.pa, "y")
      )
    }

    # Check for spatial effects
    if(!is.Waiver(model$latent)){
      # Simply use a flag and rely on FRK::auto_basis() function
      # In the future we could allow more flexibility here if needed
      if(model$latent == 'poly') latent <- TRUE else latent <- FALSE
      latent <- TRUE
    } else {
      latent <- FALSE
    }

    # Get full prediction container
    full <- model$predictors
    assertthat::assert_that(nrow(full)>1) # Security check?

    # Get offset and add it to exposure
    if(!is.Waiver(model$offset)){
      # Add offset to full prediction and load vector
      # Get whichever is poipo
      ind <- which(sapply(model$biodiversity, function(z) z$type)=="poipo")
      ofs <- model$biodiversity[[ind]]$offset[, 'spatial_offset']
      ofs_pred <- model$offset[,'spatial_offset']
      df$spatial_offset <- ofs
      full$spatial_offset <- ofs_pred
    } else { ofs <- NULL; ofs_pred <- NULL }

    # Clamp?
    if( settings$get("clamp") ) full <- clamp_predictions(model, full)

    # -- #
    # Expand predictors if non-linear is specified in settings
    if(settings$get('only_linear') == FALSE){
      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','yellow',
                                                                'Non-linearity to scrampr is best introduced by adding derivates. Ignored!')
    }

    # --- #
    if(getOption("ibis.runparallel")){
      if(!foreach::getDoParRegistered()) ibis_future(cores = getOption("ibis.nthread"),
                                                     strategy = getOption("ibis.futurestrategy"))
    }
    # Depending if hyper-parameter checks should be used, specify this separately
    if( (settings$get('optim_hyperparam')) ){
      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','yellow',
                                                                'No hyperparameter optimization for scampr implemented!')
    } else {
      # Prepare and fit model depending on type
      # For LGCP and binomial
      if(model.type %in% c("PO", "PA")){
        suppressWarnings(
          fit_scampr <- try({
            scampr::scampr(
              formula = equation,
              data = df,
              quad.weights.name = "scampr.quad.size",
              model.type = model.type,
              coord.names = c("x", "y"),
              latent.po.biasing = FALSE,
              include.sre = latent,
              se = TRUE,
              maxit = params$maxit # Number iterations
            )
          },silent = FALSE)
        )
      } else {
        # Integrated model
        suppressWarnings(
          fit_scampr <- try({
            scampr::scampr(
              formula = equation,
              bias.formula = ~ 1, # FIXME: This does not work at the moment.
              pa.data = df.pa,
              data = df.po,
              model.type = model.type,
              latent.po.biasing = FALSE,
              coord.names = c("x", "y"),
              quad.weights.name = "scampr.quad.size",
              include.sre = latent,
              se = TRUE,
              maxit = params$maxit # Number iterations
            )
          },silent = FALSE)
        )
      }
    }

    if(inherits(fit_scampr, "try-error")) stop("Model failed to converge with provided input data!")

    # --- #
    # Predict spatially
    if(!settings$get('inference_only')){
      # Messenger
      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','green','Starting prediction...')
      # Set target variables to bias_value for prediction if specified
      if(!is.Waiver(settings$get('bias_variable'))){
        for(i in 1:length(settings$get('bias_variable'))){
          if(settings$get('bias_variable')[i] %notin% names(full)){
            if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','red','Did not find bias variable in prediction object!')
            next()
          }
          full[[settings$get('bias_variable')[i]]] <- settings$get('bias_value')[i]
        }
      }
      # Make a subset of non-na values
      full$rowid <- 1:nrow(full)
      full_sub <- subset(full, stats::complete.cases(full))

      # Attempt prediction
      out <- try({
        predict(object = fit_scampr,
                newdata = full_sub,
                type = params$type,
                dens = params$dens,
                include.bias.accounting = FALSE)
      },silent = TRUE)

      # FIXME NOTE:
      # Hard transform owing to model unable to do other link functions
      if(model.type=="PA" && params$type == "response") out <- ilink(out, "cloglog")

      if(!inherits(out,"try-error")){
        # Fill output with summaries of the posterior
        prediction[full_sub$rowid] <- out
        names(prediction) <- "mean"
        prediction <- terra::mask(prediction, self$get_data("template"))
      } else {
        stop("SCAMPR prediction failed!")
      }
      try({rm(out, full, full_sub)},silent = TRUE)
    } else {
      # No prediction done
      prediction <- NULL
    }
    # Compute end of computation time
    settings$set('end.time', Sys.time())
    # Also append other parameters to settings
    for(entry in names(params)) settings$set(entry, params[[entry]])
    # Also save model type
    settings$set('model.type', model.type)

    # Definition of SCAMPR Model object ----
    obj <- DistributionModel # Make a copy to set new functions

    # Partial effect functions
    # Partial effects
    obj$set("public", "partial", function(x.var = NULL, constant = NULL, variable_length = 100,
                                          values = NULL, newdata = NULL, plot = FALSE, type = NULL, ...){
      assertthat::assert_that(is.character(x.var) || is.null(x.var),
                              is.null(constant) || is.numeric(constant),
                              is.null(type) || is.character(type),
                              is.null(newdata) || is.data.frame(newdata),
                              is.numeric(variable_length))

      # Settings
      settings <- self$settings
      # Set type
      if(is.null(type)) type <- self$settings$get("type")
      type <- match.arg(type, c("link", "response"), several.ok = FALSE)
      settings$set("type", type)

      mod <- self$get_data('fit_best')
      model <- self$model
      df <- model$biodiversity[[1]]$predictors
      df <- df |> dplyr::select("x","y", dplyr::any_of(names(df)))

      # Match x.var to argument
      if(is.null(x.var)){
        x.var <- colnames(df)
      } else {
        x.var <- match.arg(x.var, colnames(df), several.ok = TRUE)
      }

      # Take newdata or generate partial dummy
      if(is.null(newdata)){

        rr <- sapply(df[, names(df) %in% model$predictors_types$predictors[model$predictors_types$type=="numeric"]],
                     function(x) range(x, na.rm = TRUE)) |> as.data.frame()

        assertthat::assert_that(nrow(rr)>1, ncol(rr)>=1)

        df_partial <- list()
        if(!is.null(values)){ assertthat::assert_that(length(values) >= 1) }
        # Add all others as constant
        if(is.null(constant)){
          for(n in names(rr)) df_partial[[n]] <- rep( mean(df[[n]], na.rm = TRUE), variable_length )
        } else {
          for(n in names(rr)) df_partial[[n]] <- rep( constant, variable_length )
        }
        # Convert list to data.frame (same class as if newdata is provided)
        df_partial <- do.call(cbind, df_partial) |> as.data.frame()
      } else {
        df_partial <- newdata |> dplyr::select(dplyr::any_of(names(df)))
      }
      # Check if x and y present as this is required for scampr
      if(!utils::hasName(df_partial, "x") || !utils::hasName(df_partial, "y")){
        df_partial$x <- seq(min(df$x),max(df$x), length.out = nrow(df_partial))
        df_partial$y <- seq(min(df$y),max(df$y), length.out = nrow(df_partial))
      }

      # create list to store results
      o <- vector(mode = "list", length = length(x.var))
      names(o) <- x.var

      # loop through x.var
      for(v in x.var) {

        df2 <- df_partial

        if(!is.null(values)){
          df2[, v] <- values
        } else {
          df2[, v] <- seq(rr[1, v], rr[2, v], length.out = variable_length)
        }

        # Correct and reset any factor levels if set
        if(any(model$predictors_types$type=="factor")){
          lvl <- levels(model$predictors[[model$predictors_types$predictors[model$predictors_types$type=="factor"]]])
          df2[model$predictors_types$predictors[model$predictors_types$type=="factor"]] <- factor(lvl[1], levels = lvl)
        }

        pred_scampr <- predict(
          mod,
          newdata = df2,
          type = settings$get('type'),
          # use.formula = "presence-only",
          dens = 'prior', # More broader predictions
          include.bias.accounting = FALSE
        )            # Also attach the partial variable

        # Add variable name for consistency
        pred_part <- data.frame("variable" = v, "mean" = pred_scampr,
                                "partial_effect" = df2[[v]])

        o[[v]] <- pred_part
        rm(pred_part)
      }
      # Combine all
      o <- do.call(what = rbind, args = c(o, make.row.names = FALSE))

      if(plot){
        # Make a plot
        g <- ggplot2::ggplot(data = o, ggplot2::aes(x = partial_effect)) +
          ggplot2::theme_classic() +
          ggplot2::geom_line(ggplot2::aes(y = mean)) +
          ggplot2::facet_wrap(. ~ variable, scales = "free") +
          ggplot2::labs(x = "Variable", y = "Partial effect")
        print(g)
      }
      # Return the data
      return(o)
    },overwrite = TRUE)

    # Spatial partial dependence plot
    obj$set("public", "spartial", function(x.var, constant = NULL, newdata = NULL,
                                           plot = TRUE, type = NULL){
      assertthat::assert_that(is.character(x.var) || is.null(x.var),
                              "model" %in% names(self),
                              is.null(constant) || is.numeric(constant),
                              is.logical(plot),
                              is.character(type) || is.null(type)
      )

      # Settings
      settings <- self$settings
      # Set type
      if(is.null(type)) type <- self$settings$get("type")
      type <- match.arg(type, c("link", "response"), several.ok = FALSE)
      settings$set("type", type)

      mod <- self$get_data('fit_best')
      model <- self$model

      # Check if newdata is defined, if yes use that one instead
      if(!is.null(newdata)){
        df <- newdata
        assertthat::assert_that(nrow(df) == nrow(model$biodiversity[[1]]$predictors))
      } else {
        # df <- model$biodiversity[[1]]$predictors
        df <- model$predictors
      }
      df <- df |> dplyr::select("x","y", dplyr::any_of(names(df)))

      # Match x.var to argument
      if(is.null(x.var)){
        x.var <- colnames(df)
      } else {
        x.var <- match.arg(x.var, names(df), several.ok = FALSE)
      }

      # Make spatial container for prediction
      prediction <- model_to_background(model)

      # Add all others as constant
      if(is.null(constant)){
        for(n in names(df)) if(n != x.var) df[[n]] <- suppressWarnings( mean(model$predictors[[n]], na.rm = TRUE) )
      } else {
        for(n in names(df)) if(n != x.var) df[[n]] <- constant
      }

      if(any(model$predictors_types$type=="factor")){
        lvl <- levels(model$predictors[[model$predictors_types$predictors[model$predictors_types$type=="factor"]]])
        df[[model$predictors_types$predictors[model$predictors_types$type=="factor"]]] <-
          factor(lvl[1], levels = lvl)
        # FIXME: Assigning the first level (usually reference) for now. But ideally find a way to skip factors from partial predictions
      }

      # Make a subset of non-na values
      df$rowid <- 1:nrow(df)
      df_sub <- subset(df, stats::complete.cases(df))

      # Attempt prediction
      out <- try({
        predict(object = mod,
                newdata = df_sub,
                type = settings$get('type'),
                dens = settings$get('dens'),
                include.bias.accounting = FALSE)
      },silent = TRUE)

      # FIXME NOTE:
      # Hard transform owing to model unable to do other link functions
      if(settings$get('model.type') =="PA" && settings$get('type') == "response") out <- ilink(out, "cloglog")

      if(!inherits(out,"try-error")){
        # Fill output with summaries of the posterior
        prediction[df_sub$rowid] <- out
        names(prediction) <- "mean"
      } else {
        stop("Spartial prediction of scampr failed...")
      }

      # Do plot and return result
      if(plot){
        terra::plot(prediction, col = ibis_colours$ohsu_palette,
                    main = paste0("Spartial effect of ", x.var, collapse = ","))
      }
      return(prediction)
    },overwrite = TRUE)

    # Convergence check
    obj$set("public", "has_converged", function(){
      obj <- self$get_data("fit_best")
      return( obj$convergence )
    }, overwrite = TRUE)

    # Residual function
    obj$set("public", "get_residuals", function(type = NULL){
      # Get best object
      obj <- self$get_data("fit_best")
      if(is.Waiver(obj)) return(obj)
      settings <- self$settings
      if(is.null(type)){
        type <- settings$get('type')
        if(type == "response"){
          type <- "raw"
        } else {
          type <- "inverse"
        }
      }
      type <- match.arg(type, c("raw", "inverse", "pearson"),several.ok = FALSE)
      # Calculate residuals
      rd <- stats::residuals(obj, type = type)
      return(rd)
    }, overwrite = TRUE)

    # Get coefficients
    obj$set("public", "get_coefficients", function(){
      # Returns a vector of the coefficients with direction/importance
      obj <- self$get_data("fit_best")
      model <- self$model
      suppressWarnings(
        cofs <- data.frame("Feature" = names(obj$coefficients), "Beta" = stats::coef(obj)) |>
          dplyr::bind_cols(stats::confint(obj))
      )
      names(cofs) <- make.names(names(cofs))
      # Subset only to variables in predictor frame
      cofs <- subset(cofs, Feature %in% c('(Intercept)', model$predictors_names))
      return(cofs)
    }, overwrite = TRUE)

    # Function for plotting spartials
    obj$set("public", "plot_spatial", function(){
      # Get model for object
      model <- self$model
      if(!is.Waiver(model$latent)){
        stop("Plotting of spatial field not yet implemented!")
      }
    }, overwrite = TRUE)

    # Engine-specific projection function
    obj$set("public", "project", function(newdata, type = NULL, layer = "mean"){
      assertthat::assert_that("model" %in% names(self),
                              nrow(newdata) > 0,
                              all( c("x", "y") %in% names(newdata) ),
                              is.character(type) || is.null(type)
      )

      # Settings
      settings <- self$settings
      # Set type
      if(is.null(type)) type <- self$settings$get("type")
      type <- match.arg(type, c("link", "response"), several.ok = FALSE)
      settings$set("type", type)

      mod <- self$get_data('fit_best')
      model <- self$model

      # Clamp?
      if( settings$get("clamp") ) newdata <- clamp_predictions(model, newdata)

      # Set target variables to bias_value for prediction if specified
      if(!is.Waiver(settings$get('bias_variable'))){
        for(i in 1:length(settings$get('bias_variable'))){
          if(settings$get('bias_variable')[i] %notin% names(newdata)){
            if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','red','Did not find bias variable in prediction object!')
            next()
          }
          newdata[[settings$get('bias_variable')[i]]] <- settings$get('bias_value')[i]
        }
      }

      # Prediction container
      prediction <- model_to_background(model)

      df <- newdata
      df$rowid <- 1:nrow(df)
      if(!is.Waiver(model$offset)) ofs <- model$offset else ofs <- NULL
      assertthat::assert_that(nrow(df)>0)
      # Make a subset of non-na values
      df_sub <- subset(df, stats::complete.cases(df))

      # Attempt prediction
      out <- try({
        predict(object = mod,
                newdata = df_sub,
                type = settings$get('type'),
                dens = settings$get('dens'),
                include.bias.accounting = FALSE)
      },silent = TRUE)

      # FIXME NOTE:
      # Hard transform owing to model unable to do other link functions
      if(settings$get('model.type') =="PA" && settings$get('type') == "response") out <- ilink(out, "cloglog")

      if(!inherits(out,"try-error")){
        # Fill output with summaries of the posterior
        prediction[df_sub$rowid] <- out
        names(prediction) <- layer
      } else {
        stop("Projection of scampr failed...")
      }

      return(prediction)
    }, overwrite = TRUE)

    # Now generate output
    out <- obj$new(name = "SCAMPR-Model")
    out$id <- model$id
    out$model <- model
    out$settings <- settings
    out$fits <- list(
      "fit_best" = fit_scampr,
      "fit_best_equation" = equation,
      "prediction" = prediction
    )

    return(out)
  },overwrite = TRUE)

  # Define engine object and save
  eg <- eg$new(engine = "SCAMPR-Engine", name = "<SCAMPR>")
  eg$data <- list(
    'template' = template,
    'params' = params
  )

  # Set engine in distribution object
  y <- x$clone(deep = TRUE)
  return( y$set_engine(eg) )
} # End of function
