#' Built formula for GDB model
#'
#' @description
#' This function built a formula for a `engine_gdb()` model.
#' @param model A [`list()`] object containing the prepared model data for a given biodiversity dataset.
#' @param x A [`BiodiversityDistribution`] object.
#' @param id The id for the species formula.
#' @param settings A [`Settings`] object.
#' @author Martin Jung
#' @note Function is not meant to be run outside the train() call.
#' @keywords internal
#' @noRd
built_formula_gdb <- function(model, id, x, settings){
  assertthat::assert_that(
    is.list(model),
    length(model) > 0,
    assertthat::has_name(model, "predictors"),
    inherits(x, "BiodiversityDistribution"),
    inherits(settings, 'Settings'),
    is.character(id) || is.Id(id),
    msg = "Error in model object. This function is not meant to be called outside ouf train()."
  )
  # Get object for id
  obj <- model$biodiversity[[id]]
  # Extract basic stats from the model object
  types <- as.character( sapply( model$biodiversity, function(x) x$type ) )
  fams <- as.character( sapply( model$biodiversity, function(z) z$family ) )
  bionames = sapply(model$biodiversity, function(x) x$name)
  ids <- names(model$biodiversity)
  priors <- model$priors

  # Default equation found
  if(obj$equation == '<Default>' || is.Waiver(obj$equation)){
    # Construct formula with all variables
    form <- "observed ~ "

    # Use only variables that have sufficient covariate range for training
    # Finally check that a minimum of unique numbers are present in the predictor range and if not, remove them
    covariates <- rm_insufficient_covs(model = obj, tr = 4)
    if(is.null(covariates)) stop("Not enough variance in training data to fit a SDM!")

    if(!is.Waiver(priors)){
      # Loop through all provided GDB priors
      supplied_priors <- as.vector(priors$varnames())
      for(v in supplied_priors){
        if(v %notin% covariates) next() # In case the variable has been removed
        # First add linear effects
        form <- paste(form, paste0('bmono(', v,
                                   ', constraint = \'', priors$get(v) ,'\'',
                                   ')', collapse = ' + ' ), ' + ' )
      }
      # Add linear and smooth effects for all missing ones
      miss <- covariates[covariates %notin% supplied_priors]
      if(length(miss)>0){
        # Add linear predictors
        form <- paste(form, paste0('bols(',miss,')', collapse = ' + '))
        if(!settings$get('only_linear')){
          # And smooth effects for all numeric data
          miss <- miss[ miss %in% obj$predictors_types$predictors[which(obj$predictors_types$type=="numeric")] ]
          form <- paste(form, ' + ', paste0('bbs(', miss,', knots = 4)',
                                            collapse = ' + '
          ))
        }
      }
    } else {
      # Add linear predictors
      form <- paste(form, paste0('bols(',covariates,')',collapse = ' + '))
      if(settings$get('only_linear') == FALSE){
        # And smooth effects
        form <- paste(form, ' + ', paste0('bbs(',
                                          covariates[which(covariates %in% obj$predictors_types$predictors[obj$predictors_types$type=="numeric"] )],', knots = 4)',
                                          collapse = ' + '
        ))
      }
      # Add also random effect if there are any factors? THIS currently crashes when there are too few factors
      # if(any(model$predictors_types$type=="factor")){
      #   form <- paste(form, ' + ' ,paste0('brandom(',
      #                              model$biodiversity[[id]]$predictors_types$predictors[which(model$biodiversity[[id]]$predictors_types$type == 'factor')],
      #                              ')',collapse = " + "))
      # }
    }
    # Convert to formula
    form <- to_formula(form)
    # Add offset if specified
    if(!is.Waiver(model$offset) ){ form <- stats::update.formula(form, paste0('~ . + offset(spatial_offset)') ) }
    if( length( grep('Spatial',x$get_latent() ) ) > 0 ){
      # Update with spatial term
      form <- stats::update.formula(form, paste0(" ~ . + ",
                                          x$engine$get_equation_latent_spatial())
      )
    }
  } else{
    if(getOption('ibis.setupmessages')) myLog('[Estimation]','yellow','Use custom model equation')
    form <- to_formula(obj$equation)
    # If response is missing, add manually
    if(attr(stats::terms(form), "response")==0){
      form <- stats::update.formula(form, "observed ~ .")
    }
    # Check that bols/bbs are in terms and if not add them for each variable
    if( settings$get("only_linear") ){
      if( length(grep("bols",attr(stats::terms(form2), "term.labels") ))==0 ){
        # Assume that each variable as none, so add
        form <- stats::as.formula(paste0("observed ~", paste0("bols(", obj$predictors_names, ")",collapse = " + ")))
      }
    } else {
      if( length(grep("bbs",attr(stats::terms(form2), "term.labels") ))==0 ){
        # Assume that each variable as none, so add
        form <- stats::as.formula(paste0("observed ~", paste0("bbs(", obj$predictors_names, ")",collapse = " + ")))
      }
    }

    assertthat::assert_that(
      is.formula(form),
      attr(stats::terms(form), "response")==1, # Has Response
      all( all.vars(form) %in% c('observed', obj[['predictors_names']]) )
    )
  }

  return(form)
}

#' Use a fitted model for creating a new class prediction in raster form
#'
#' @param fit A fitted [`mboost`] model with [`binomial`] distribution.
#' @param nd A new data.frame with all predictiors used in fit.
#' @param template A [`SpatRaster`] object that can be used as spatial template.
#' @returns A [`RasterLayer`] containing a presence-absence prediction.
#' @keywords utils
#' @noRd
predict_gdbclass <- function(fit, nd, template){
  assertthat::assert_that(
    inherits(fit, 'mboost'),
    !is.null( grep('Binomial', fit$family@name,ignore.case = TRUE) ),
    is.data.frame(nd),
    is.Raster(template)
  )
  # Redo a template to be sure
  template <- emptyraster(template)

  # Remove missing data in newdata data.frame
  nd$cellid <- rownames(nd)
  nd <- subset(nd, stats::complete.cases(nd))

  suppressWarnings(
    pred_gdb <- mboost::predict.mboost(object = fit, newdata = nd,
                                       type = 'class', aggregate = 'sum')
  )
  # Fill output
  prediction <- emptyraster(template)
  prediction[as.numeric(nd$cellid)] <- pred_gdb
  prediction[prediction < max(as.numeric(pred_gdb),na.rm = T)] <- 0; prediction[prediction >0] <- 1
  names(prediction) <- 'presence'
  return(prediction)
}

#' Sanitize mboost summary output
#'
#' @description
#' Extracts the coefficients and selection frequencies from a [`mboost`] model.
#' @param obj A fitted [`mboost`] object.
#' @noRd
#' @keywords internal
clean_mboost_summary <- function(obj){
  assertthat::assert_that(
    inherits(obj, "mboost")
  )

  # Get Variable importance
  vi <- mboost::varimp( obj )
  vi <- sort( vi[which(vi>0)],decreasing = TRUE )

  # Get coefficients
  co <- mboost::extract(obj, "coefficient")
  co <- co[names(vi)]
  assertthat::assert_that(all(names(vi) == names(co)))
  # Calculate coefficient. If smooth, calculate mean of the differential between knots
  co <- sapply(co, function(x) {
    if(length(x)>2){
      mean(diff(x))
    } else x[2]
  })

  # Now split up the names and types
  types <- sapply(names(vi), function(z) strsplit(z, "\\(")[[1]][1]) |> as.vector()
  vars <- sapply(names(vi), function(z) strsplit(z, "\\(")[[1]][2]) |> as.vector()
  vars <- gsub(vars, pattern="\\)",replacement="")

  # Construct tibble
  out <- tibble::tibble(
    variable = vars,
    type = types,
    varimp = vi,
    beta = co
  )
  return(out)
}

#' Check training predictor complexity
#'
#' @description
#' This internal function tests an existing set of covariates (usually the training data)
#' on the number of unique values within.
#' If fewer values than a given threshold (\code{'tr'}) is detected, then the predictor is removed, thus
#' reducing complexity.
#' @note Maybe in the future a more cleverer solution could be thought of, for instance using a singular value decomposition?
#' @param model A [`list`] of a model object containing the various predictors and biodiversity occurrence information.
#' @param tr A [`numeric`] value describing a threshold of minimum unique values to be retained.
#' @returns A [`vector`] with the variables that full fill the threshold.
#' @keywords internal
#' @noRd
rm_insufficient_covs <- function(model, tr = 5){
  assertthat::assert_that(
    is.list(model),
    is.numeric(tr),
    tr > 0
  )

  # Check that biodiversity information is present
  assertthat::assert_that(
    assertthat::has_name(model, "observations"),
    assertthat::has_name(model, "predictors_names")
  )

  # Now get all continuous ones
  vars_num <- model$predictors_types$predictors[model$predictors_types$type=="numeric"]
  vars_fac <- model$predictors_types$predictors[model$predictors_types$type=="factor"]

  vars_uniques <- apply(model$predictors[,vars_num], 2, function(x) length(unique(x,na.rm = TRUE)) )

  # Get all variables smaller than the threshold and return the original data.frame without them
  sufficient <- which(vars_uniques >= tr)

  # Get the factor variables in it as well
  if(length(vars_fac)>0){
    vars_uniques <- apply(model$predictors[vars_fac], 2, function(x) length(unique(x,na.rm = TRUE)) )
    # Get all factor variables with at least 2 levels
    sufficient_fac <- which(vars_uniques >= 2)
    if(length(sufficient_fac)>0){
      sufficient <- c(names(sufficient), names(sufficient_fac))
    }
  } else {
    sufficient <- names(sufficient)
  }
  assertthat::assert_that(all(sufficient %in% model$predictors_names)) # This should return a character of covariate names
  if(length(sufficient)==0){
    return(NULL)
  } else {
    return(sufficient)
  }
}

#' Calculate weights for Point Process models
#'
#' @param df The [`data.frame`] for which weights are to be calculated.
#' @param presence A [`vector`] with the observed species. Has to be in range \code{0} to \code{Inf}.
#' @param bg A background [`raster`] layer.
#' @param use_area A [`logical`] on whether area is to be used instead of grid counts.
#' @param weight A [`numeric`] weight to be used in down-weighted regressions.
#' @param type Accepting either “Infinitely weighted logistic regression” \code{'IWLR'} for use with binomial
#' logistic regressions or “Down-weighted Poisson regression” \code{'DWPR'} (Default).
#' @references
#' * Renner, I.W., Elith, J., Baddeley, A., Fithian, W., Hastie, T., Phillips, S.J., Popovic, G. and Warton, D.I., 2015. Point process models for presence‐only analysis. Methods in Ecology and Evolution, 6(4), pp.366-379.
#' * Fithian, W. & Hastie, T. (2013) Finite-sample equivalence in statistical models for presence-only data. The Annals of Applied Statistics 7, 1917–1939
#' @return A vector with the weights
#' @keywords utils
#' @noRd
ppm_weights <- function(df, pa, bg, use_area = FALSE, weight = 1e-6, type = "DWPR"){
  assertthat::assert_that(
    is.data.frame(df),
    length(unique(pa)) > 1,
    nrow(df) == length(pa),
    is.logical(use_area),
    is.numeric(weight),
    is.character(type)
  )
  type <- match.arg(type, c("DWPR", "IWLR"),several.ok = FALSE)

  if(use_area){
    suppressWarnings( ar <- terra::cellSize(bg) )
    ar <- terra::mask(ar, bg)
    nc <- terra::global(ar, "sum", na.rm = TRUE)[,1]
  } else {
    # number of non-NA cells
    nc <- terra::global(!is.na(bg), "sum", na.rm = TRUE)[,1]
  }

  # Set output weight as default
  if(type == "DWPR"){
    w <- rep( weight, nrow(df) )
    w[which(pa == 0)] = nc / sum(pa == 0)
  } else {
    w = (10^6)^(1 - pa)
  }

  assertthat::assert_that(
    length(unique(w)) > 1,
    length(w) == nrow(df),
    is.vector(w) && is.numeric(w)
  )
  return(w)
}
