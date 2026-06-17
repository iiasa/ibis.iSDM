#' Built formula for STAN model
#'
#' @description This function built a formula for a `engine_stan()` model.
#'
#' @param model A [`list()`] object containing the prepared model data for a given
#' biodiversity dataset.
#' @param x A [`BiodiversityDistribution`] object.
#' @param id The id for the species formula.
#' @param settings A [`Settings`] object.
#'
#' @note Function is not meant to be run outside the train() call.
#'
#' @author Martin Jung
#'
#' @noRd
#'
#' @keywords internal
built_formula_stan <- function(model, id, x, settings){
  assertthat::assert_that(
    is.list(model),
    length(model) > 0,
    assertthat::has_name(model, "predictors"),
    inherits(x, "BiodiversityDistribution"),
    inherits(settings, 'Settings'),
    is.character(id) || is.Id(id),
    msg = "Error in model object. This function is not meant to be called outside of train()."
  )
  obj <- model$biodiversity[[id]]

  # Default equation found (e.g. no separate specification of effects). Stan
  # assembles offsets and dataset intercepts directly, so this formula is only
  # used to declare the ecological covariate terms.
  if(is.character(obj$equation) && obj$equation == '<Default>'){
    form <- stats::reformulate(obj$predictors_names,
                               response = "observed",
                               intercept = FALSE)
  } else {
    if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','yellow','Use custom model equation.')
    form <- to_formula(obj$equation)
    assertthat::assert_that(
      all(formula_terms(form) %in% model[['predictors_names']]),
      msg = "Predictors in custom Stan formula not found!"
    )
    form <- stats::reformulate(formula_terms(form),
                               response = "observed",
                               intercept = FALSE)
  }
  return(form)
}

#' Checks whether cmdstanr is available and otherwise tries to install it
#'
#' @param install A [`logical`] factor to indicate whether cmdstanr should be
#' directly installed (Default: \code{TRUE}).
#' @param ask [`logical`] whether the cmdstanr package is to be installed
#' (Default: \code{FALSE}).
#'
#' @keywords stan utils
#'
#' @noRd
#'
#' @keywords internal
stan_check_cmd <- function(install = TRUE, ask = FALSE){
  assertthat::assert_that(
    is.logical(install), is.logical(ask)
  )
  # Check if available
  if(!requireNamespace("cmdstanr", quietly = TRUE)){
    if(install){
      if(ask){ a <- utils::askYesNo("Install cmdstanr?") } else { a <- TRUE}
      if(a){
        utils::install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
        cmdstanr::check_cmdstan_toolchain()
        cmdstanr::install_cmdstan(cores = 2)
      }
    } else {
      check_package("cmdstanr")
    }
  } else {
    invisible()
  }
}

#' Prepare Stan input for shared-beta multi-dataset SDMs
#'
#' @keywords internal
#' @noRd
stan_prepare_model <- function(model, settings, template){
  assertthat::assert_that(
    is.list(model),
    assertthat::has_name(model, "biodiversity"),
    assertthat::has_name(model, "predictors"),
    inherits(settings, "Settings"),
    is.Raster(template)
  )

  types <- vapply(model$biodiversity, function(z) z$type, character(1))
  families <- vapply(model$biodiversity, function(z) z$family, character(1))
  assertthat::assert_that(
    all(types %in% c("poipo", "poipa")),
    msg = "Stan supports point presence-only and point presence-absence datasets only."
  )
  assertthat::assert_that(
    all((types == "poipo" & families == "poisson") |
          (types == "poipa" & families == "binomial")),
    msg = "Stan supports poipo/poisson and poipa/binomial combinations only."
  )
  assertthat::assert_that(
    length(unique(vapply(model$biodiversity, function(z) z$use_intercept, logical(1)))) == 1,
    msg = "'separate_intercept' must be identical for all Stan datasets."
  )

  model <- explode_factor_predictors(model, expand_raster = TRUE, verbose = TRUE)
  assertthat::assert_that(!is.na(terra::global(template, "min", na.rm = TRUE)[, 1]))

  for(i in seq_along(model$biodiversity)){
    bd <- explode_factor_predictors(model$biodiversity[[i]], expand_raster = FALSE)

    if(bd$type == "poipo"){
      presabs <- add_pseudoabsence(df = bd$observations,
                                   field_occurrence = 'observed',
                                   template = template,
                                   settings = bd$pseudoabsence_settings)
      if(inherits(presabs, 'sf')) presabs <- sf::st_drop_geometry(presabs)

      abs <- subset(presabs, observed == 0)
      if(nrow(abs) > 0){
        envs <- get_rastervalue(coords = abs[, c('x', 'y')],
                                env = model$predictors_object$get_data(df = FALSE),
                                rm.na = FALSE)
        if(utils::hasName(bd$predictors, "Intercept")) envs$Intercept <- 1
      } else {
        envs <- bd$predictors[0, , drop = FALSE]
      }

      base_cols <- unique(c('x', 'y', 'Intercept', bd$predictors_names))
      base_cols <- base_cols[base_cols %in% names(bd$predictors)]
      pres <- bd$predictors[, base_cols, drop = FALSE]
      envs <- envs[, base_cols, drop = FALSE]
      df <- rbind(pres, envs)

      user_weight <- bd$expect
      if(length(user_weight) != nrow(presabs)){
        user_weight <- c(user_weight, rep(1, nrow(presabs) - length(user_weight)))
      }
      keep <- stats::complete.cases(df[, c('x', 'y', bd$predictors_names), drop = FALSE])
      df <- df[keep, , drop = FALSE]
      presabs <- presabs[keep, , drop = FALSE]
      user_weight <- user_weight[keep]
      assertthat::assert_that(nrow(df) == nrow(presabs), length(user_weight) == nrow(df))

      quad <- ppm_weights(df = df,
                          pa = presabs[['observed']],
                          bg = template,
                          weight = 1e-6,
                          type = "DWPR")
      quad[!is.finite(quad) | quad <= 0] <- .Machine$double.eps
      user_weight[!is.finite(user_weight) | user_weight <= 0] <- .Machine$double.eps

      spatial_offset <- rep(0, nrow(df))
      if(!is.Waiver(model$offset)){
        if(assertthat::has_name(model, "offset_object") &&
           !is.null(model$offset_object) &&
           !is.Waiver(model$offset_object)){
          ofs <- get_rastervalue(coords = df[, c("x", "y")],
                                 env = model$offset_object,
                                 rm.na = FALSE)
          names(ofs)[which(names(ofs) == names(model$offset_object))] <- "spatial_offset"
          if(utils::hasName(ofs, "spatial_offset")) spatial_offset <- ofs[["spatial_offset"]]
        } else if(!is.null(bd$offset) &&
                  is.data.frame(bd$offset) &&
                  utils::hasName(bd$offset, "spatial_offset") &&
                  nrow(bd$offset) == nrow(df)){
          spatial_offset <- bd$offset[["spatial_offset"]]
        }
      }
      spatial_offset[!is.finite(spatial_offset)] <- 0

      bd$observations <- presabs
      bd$predictors <- df
      bd$stan_offset <- log(quad) + as.numeric(spatial_offset)
      bd$stan_weight <- as.numeric(user_weight)
      bd$expect <- bd$stan_weight

      pres_points <- bd$observations[bd$observations$observed > 0, c("x", "y"), drop = FALSE]
      pres_raster <- terra::rasterize(x = guess_sf(pres_points),
                                      y = template, fun = 'count', background = 0)
      full_quad <- ppm_weights(df = model$predictors,
                               pa = as.vector(pres_raster[]),
                               bg = template,
                               weight = 1,
                               type = "DWPR")
      full_quad[!is.finite(full_quad) | full_quad <= 0] <- .Machine$double.eps
      bd$stan_projection_offset <- log(full_quad)
      bd$stan_default_offset <- stats::median(log(quad[presabs$observed == 0]), na.rm = TRUE)
      if(!is.finite(bd$stan_default_offset)){
        bd$stan_default_offset <- stats::median(log(quad), na.rm = TRUE)
      }
    } else {
      bd$observations[['observed']] <- ifelse(bd$observations[['observed']] >= 1, 1, 0)
      prNum <- sum(bd$observations[['observed']] == 1)
      bgNum <- sum(bd$observations[['observed']] == 0)
      class_weight <- ifelse(bd$observations[['observed']] == 1,
                             1,
                             ifelse(bgNum > 0, prNum / bgNum, 1))
      user_weight <- as.numeric(class_weight * bd$expect)
      user_weight[!is.finite(user_weight) | user_weight <= 0] <- .Machine$double.eps

      spatial_offset <- rep(0, nrow(bd$predictors))
      if(!is.Waiver(model$offset)){
        if(assertthat::has_name(model, "offset_object") &&
           !is.null(model$offset_object) &&
           !is.Waiver(model$offset_object)){
          ofs <- get_rastervalue(coords = bd$predictors[, c("x", "y")],
                                 env = model$offset_object,
                                 rm.na = FALSE)
          names(ofs)[which(names(ofs) == names(model$offset_object))] <- "spatial_offset"
          if(utils::hasName(ofs, "spatial_offset")) spatial_offset <- ofs[["spatial_offset"]]
        } else if(!is.null(bd$offset) &&
                  is.data.frame(bd$offset) &&
                  utils::hasName(bd$offset, "spatial_offset") &&
                  nrow(bd$offset) == nrow(bd$predictors)){
          spatial_offset <- bd$offset[["spatial_offset"]]
        }
      }
      spatial_offset[!is.finite(spatial_offset)] <- 0

      bd$stan_offset <- as.numeric(spatial_offset)
      bd$stan_weight <- user_weight
      bd$expect <- user_weight
    }

    model$biodiversity[[i]] <- bd
  }

  feature_sets <- lapply(model$biodiversity, function(z) z$predictors_names)
  if(length(unique(vapply(feature_sets, paste, character(1), collapse = "\r"))) > 1){
    myLog('[Estimation]', 'yellow',
          'Stan datasets use different formulas or predictor sets. Combining predictors into one shared beta set.')
  }
  feature_names <- unique(unlist(feature_sets, use.names = FALSE))
  missing_full_features <- feature_names[!feature_names %in% names(model$predictors)]
  assertthat::assert_that(
    length(missing_full_features) == 0,
    msg = paste0("Stan predictors missing from the full prediction grid: ",
                 paste(missing_full_features, collapse = ", "))
  )
  assertthat::assert_that(length(feature_names) > 0)

  for(i in seq_along(model$biodiversity)){
    bd <- model$biodiversity[[i]]
    missing_features <- feature_names[!feature_names %in% names(bd$predictors)]
    if(length(missing_features) > 0){
      assertthat::assert_that(
        assertthat::has_name(model, "predictors_object"),
        !is.null(model$predictors_object),
        !is.Waiver(model$predictors_object),
        msg = paste0("Cannot recover missing Stan predictors: ",
                     paste(missing_features, collapse = ", "))
      )
      env <- get_rastervalue(coords = bd$predictors[, c("x", "y")],
                             env = model$predictors_object$get_data(df = FALSE),
                             rm.na = FALSE)
      for(name in missing_features){
        assertthat::assert_that(name %in% names(env),
                                msg = paste0("Missing Stan predictor not found in raster data: ", name))
        bd$predictors[[name]] <- env[[name]]
      }
    }
    keep <- stats::complete.cases(bd$predictors[, feature_names, drop = FALSE])
    if(any(!keep)){
      bd$predictors <- bd$predictors[keep, , drop = FALSE]
      bd$observations <- bd$observations[keep, , drop = FALSE]
      bd$stan_offset <- bd$stan_offset[keep]
      bd$stan_weight <- bd$stan_weight[keep]
      bd$expect <- bd$expect[keep]
    }
    bd$predictors_names <- feature_names
    model$biodiversity[[i]] <- bd
  }

  model$predictors_names <- feature_names
  pt <- model$predictors_types[match(feature_names, model$predictors_types$predictors), , drop = FALSE]
  missing_pt <- is.na(pt$predictors)
  if(any(missing_pt)){
    pt[missing_pt, ] <- data.frame(predictors = feature_names[missing_pt], type = "numeric")
  }
  model$predictors_types <- pt

  model$stan <- list(
    components = stan_component_table(model),
    feature_names = feature_names,
    use_dataset_intercepts = unique(vapply(model$biodiversity, function(z) z$use_intercept, logical(1)))
  )
  model$stan$data <- stan_model_data(model, settings)
  model$stan$code <- stan_model_code(model, settings)
  model
}

#' Stan component table
#'
#' @keywords internal
#' @noRd
stan_component_table <- function(model){
  assertthat::assert_that(
    is.list(model),
    assertthat::has_name(model, "biodiversity"),
    length(model$biodiversity) > 0
  )

  ids <- names(model$biodiversity)
  if(is.null(ids)) ids <- as.character(seq_along(model$biodiversity))
  types <- vapply(model$biodiversity, function(x) x$type, character(1))
  families <- vapply(model$biodiversity, function(x) x$family, character(1))

  assertthat::assert_that(
    all(types %in% c("poipo", "poipa")),
    msg = "Stan supports point presence-only and point presence-absence datasets only."
  )
  assertthat::assert_that(
    all((types == "poipo" & families == "poisson") |
          (types == "poipa" & families == "binomial")),
    msg = "Stan supports poipo/poisson and poipa/binomial combinations only."
  )

  mixed_po_pa <- any(types == "poipo") && any(types == "poipa")
  links <- vapply(seq_along(model$biodiversity), function(i){
    bd <- model$biodiversity[[i]]
    if(bd$type == "poipo") return("log")
    li <- bd$link
    if(is.Waiver(li) || is.null(li) || length(li) == 0){
      if(mixed_po_pa) "cloglog" else "logit"
    } else {
      match.arg(li, c("logit", "cloglog"), several.ok = FALSE)
    }
  }, character(1))

  data.frame(
    index = seq_along(model$biodiversity),
    id = as.character(ids),
    name = vapply(model$biodiversity, function(x) x$name, character(1)),
    type = types,
    family = families,
    link = links,
    family_code = ifelse(types == "poipo", 1L, 2L),
    link_code = match(links, c("log", "logit", "cloglog")),
    stringsAsFactors = FALSE
  )
}

#' Resolve a Stan prediction component
#'
#' @keywords internal
#' @noRd
stan_component <- function(model, component = NULL){
  components <- model$stan$components
  assertthat::assert_that(is.data.frame(components), nrow(components) > 0)

  if(is.null(component)){
    po <- which(components$type == "poipo")
    return(if(length(po) > 0) po[1] else 1L)
  }

  if(is.numeric(component)){
    assertthat::assert_that(
      length(component) == 1,
      component >= 1,
      component <= nrow(components),
      msg = "Stan component index out of range."
    )
    return(as.integer(component))
  }

  component <- as.character(component)
  wh <- which(components$id %in% component | components$name %in% component)
  assertthat::assert_that(
    length(wh) == 1,
    msg = paste0(
      "Could not resolve Stan component '", component,
      "'. Use a dataset name, dataset id, or integer index."
    )
  )
  wh
}

#' Formula used for Stan posterior prediction
#'
#' @keywords internal
#' @noRd
stan_prediction_formula <- function(model){
  assertthat::assert_that(
    assertthat::has_name(model, "stan"),
    length(model$stan$feature_names) > 0
  )
  stats::reformulate(model$stan$feature_names, response = "observed", intercept = FALSE)
}

#' Align projection data to Stan feature names
#'
#' @keywords internal
#' @noRd
stan_align_newdata <- function(newdata, model){
  assertthat::assert_that(is.data.frame(newdata))
  feature_names <- model$stan$feature_names

  missing_features <- feature_names[!feature_names %in% names(newdata)]
  if(length(missing_features) > 0){
    for(col in names(newdata)){
      if(!is.factor(newdata[[col]]) && !is.character(newdata[[col]])) next
      z <- explode_factor(as.factor(newdata[[col]]), name = col)
      if(any(colnames(z) %in% missing_features)){
        newdata[[col]] <- NULL
        newdata <- cbind(newdata, z)
      }
    }
  }

  missing_features <- feature_names[!feature_names %in% names(newdata)]
  assertthat::assert_that(
    length(missing_features) == 0,
    msg = paste0(
      "Feature name mismatch between Stan model and newdata. Missing: ",
      paste(missing_features, collapse = ", ")
    )
  )
  newdata
}

#' Build link-scale offsets for a selected Stan component
#'
#' @keywords internal
#' @noRd
stan_prediction_offset <- function(model, newdata, component = NULL, offset = NULL){
  component <- stan_component(model, component)
  components <- model$stan$components
  out <- rep(0, nrow(newdata))

  if(components$type[component] == "poipo"){
    bd <- model$biodiversity[[component]]
    if(!is.null(bd$stan_projection_offset) &&
       length(bd$stan_projection_offset) == nrow(newdata) &&
       nrow(newdata) == nrow(model$predictors)){
      out <- bd$stan_projection_offset
    } else {
      default_offset <- bd$stan_default_offset
      if(is.null(default_offset) || !is.finite(default_offset)) default_offset <- 0
      out <- rep(default_offset, nrow(newdata))
    }
  }

  if(!is.null(offset)){
    assertthat::assert_that(length(offset) == nrow(newdata))
    offset[!is.finite(offset)] <- 0
    out <- out + as.numeric(offset)
  } else if(!is.Waiver(model$offset)){
    spatial_offset <- rep(0, nrow(newdata))
    if(is.data.frame(model$offset) &&
       utils::hasName(model$offset, "spatial_offset") &&
       nrow(model$offset) == nrow(newdata)){
      spatial_offset <- model$offset[["spatial_offset"]]
    } else if(assertthat::has_name(model, "offset_object") &&
              !is.null(model$offset_object) &&
              !is.Waiver(model$offset_object) &&
              all(c("x", "y") %in% names(newdata))){
      ofs <- get_rastervalue(coords = newdata[, c("x", "y")],
                             env = model$offset_object,
                             rm.na = FALSE)
      if(ncol(ofs) > 0){
        names(ofs)[which(names(ofs) == names(model$offset_object))] <- "spatial_offset"
        if(utils::hasName(ofs, "spatial_offset")) spatial_offset <- ofs[["spatial_offset"]]
      }
    } else if(is.data.frame(model$offset) &&
              utils::hasName(model$offset, "spatial_offset")){
      spatial_offset <- rep(mean(model$offset[["spatial_offset"]], na.rm = TRUE), nrow(newdata))
    }
    spatial_offset[!is.finite(spatial_offset)] <- 0
    out <- out + as.numeric(spatial_offset)
  }

  out[!is.finite(out)] <- 0
  as.numeric(out)
}

#' Predict from a prepared Stan SDM
#'
#' @keywords internal
#' @noRd
stan_predict <- function(fit, model, settings, newdata, component = NULL,
                         offset = NULL, type = NULL){
  if(is.null(type)) type <- settings$get("type")
  if(type == "link") type <- "predictor"
  type <- match.arg(type, c("predictor", "response"), several.ok = FALSE)

  component <- stan_component(model, component)
  components <- model$stan$components
  newdata <- stan_align_newdata(newdata, model)

  if(isTRUE(settings$get("clamp"))){
    newdata <- clamp_predictions(model, newdata)
  }
  if(!is.Waiver(settings$get('bias_variable'))){
    for(i in seq_along(settings$get('bias_variable'))){
      if(!(settings$get('bias_variable')[i] %in% names(newdata))) next
      newdata[[settings$get('bias_variable')[i]]] <- settings$get('bias_value')[i]
    }
  }

  pred_offset <- stan_prediction_offset(model, newdata, component, offset)
  intercept_index <- if(model$stan$use_dataset_intercepts) component else 1L
  posterior_predict_stanfit(
    obj = fit,
    form = stan_prediction_formula(model),
    newdata = newdata,
    offset = pred_offset,
    family = components$family[component],
    link = components$link[component],
    intercept = paste0("Intercept[", intercept_index, "]"),
    feature_names = model$stan$feature_names,
    type = type
  )
}

#' Build Stan model code for shared-beta multi-dataset SDMs
#'
#' @keywords internal
#' @noRd
stan_model_code <- function(model, settings){
  assertthat::assert_that(
    assertthat::has_name(model, "stan"),
    inherits(settings, "Settings")
  )
  feature_names <- model$stan$feature_names
  components <- model$stan$components
  assertthat::assert_that(length(feature_names) > 0)

  optim_hyperparam <- isTRUE(settings$get("optim_hyperparam")) && is.Waiver(model$priors)
  template <- if(nrow(components) == 1 && components$type[1] == "poipo"){
    "model_single_poipo.stan"
  } else if(nrow(components) == 1 && components$type[1] == "poipa"){
    "model_single_poipa.stan"
  } else {
    "model_multi_integrated.stan"
  }
  stan_file <- system.file(file.path("stanfiles", template),
                           package = "ibis.iSDM",
                           mustWork = TRUE)
  data_file <- system.file("stanfiles/data_parameters.stan",
                           package = "ibis.iSDM",
                           mustWork = TRUE)
  other_functions_file <- system.file("stanfiles/other_functions.stan",
                                      package = "ibis.iSDM",
                                      mustWork = TRUE)
  code <- paste(readLines(stan_file, warn = FALSE), collapse = "\n")
  data_parameters <- paste(readLines(data_file, warn = FALSE), collapse = "\n")
  other_functions <- paste(readLines(other_functions_file, warn = FALSE), collapse = "\n")

  if(optim_hyperparam){
    prior_file <- system.file("stanfiles/prior_functions.stan",
                              package = "ibis.iSDM",
                              mustWork = TRUE)
    functions_extra <- paste(readLines(prior_file, warn = FALSE), collapse = "\n")
    data_extra <- paste(
      "real<lower=0> hs_df;",
      "real<lower=0> hs_df_global;",
      "real<lower=0> hs_df_slab;",
      "real<lower=0> hs_scale_global;",
      "real<lower=0> hs_scale_slab;",
      sep = "\n"
    )
    parameters <- paste(
      "vector[K] zb;",
      "vector<lower=0>[K] hs_local;",
      "real<lower=0> hs_global;",
      "real<lower=0> hs_slab;",
      "vector[J_intercept] Intercept;",
      sep = "\n"
    )
    transformed_parameters <- paste(
      "vector[K] beta;",
      "beta = horseshoe(zb, hs_local, hs_global, hs_scale_slab^2 * hs_slab);",
      sep = "\n"
    )
    coefficient_priors <- paste(
      "target += std_normal_lpdf(zb);",
      "target += student_t_lpdf(hs_local | hs_df, 0, 1) - rows(hs_local) * log(0.5);",
      "target += student_t_lpdf(hs_global | hs_df_global, 0, hs_scale_global) - log(0.5);",
      "target += inv_gamma_lpdf(hs_slab | 0.5 * hs_df_slab, 0.5 * hs_df_slab);",
      sep = "\n"
    )
  } else {
    functions_extra <- ""
    data_extra <- ""
    parameters <- paste(
      "vector[K] beta;",
      "vector[J_intercept] Intercept;",
      sep = "\n"
    )
    transformed_parameters <- ""
    coefficient_priors <- character(length(feature_names))
    for(i in seq_along(feature_names)){
      if(!is.Waiver(model$priors) && feature_names[i] %in% model$priors$varnames()){
        pp <- model$priors$get(feature_names[i])
      } else {
        pp <- c(0, 2)
      }
      assertthat::assert_that(
        is.numeric(pp), length(pp) == 2, is.finite(pp[1]),
        is.finite(pp[2]), pp[2] > 0,
        msg = paste0("Invalid STANPrior for feature ", feature_names[i], ".")
      )
      coefficient_priors[i] <- paste0("target += normal_lpdf(beta[", i, "] | ", pp[1], ", ", pp[2], ");")
    }
    coefficient_priors <- paste(coefficient_priors, collapse = "\n")
  }

  replacements <- list(
    "{{data_parameters}}" = data_parameters,
    "{{other_functions}}" = other_functions,
    "{{functions_extra}}" = functions_extra,
    "{{data_extra}}" = data_extra,
    "{{parameters}}" = parameters,
    "{{transformed_parameters}}" = transformed_parameters,
    "{{coefficient_priors}}" = coefficient_priors
  )
  for(token in names(replacements)){
    code <- gsub(token, replacements[[token]], code, fixed = TRUE)
  }
  code
}

#' Build Stan data for a prepared Stan model
#'
#' @keywords internal
#' @noRd
stan_model_data <- function(model, settings){
  assertthat::assert_that(
    assertthat::has_name(model, "stan"),
    inherits(settings, "Settings")
  )
  feature_names <- model$stan$feature_names
  components <- model$stan$components
  use_dataset_intercepts <- model$stan$use_dataset_intercepts

  X <- list()
  observed <- offsets <- weights <- dataset <- intercept_id <- list()
  for(i in seq_along(model$biodiversity)){
    bd <- model$biodiversity[[i]]
    X[[i]] <- as.matrix(bd$predictors[, feature_names, drop = FALSE])
    storage.mode(X[[i]]) <- "double"
    observed[[i]] <- as.integer(bd$observations[["observed"]])
    if(bd$type == "poipa") observed[[i]] <- ifelse(observed[[i]] >= 1L, 1L, 0L)
    offsets[[i]] <- as.numeric(bd$stan_offset)
    weights[[i]] <- as.numeric(bd$stan_weight)
    weights[[i]][!is.finite(weights[[i]]) | weights[[i]] <= 0] <- .Machine$double.eps
    offsets[[i]][!is.finite(offsets[[i]])] <- 0
    dataset[[i]] <- rep(i, nrow(X[[i]]))
    intercept_id[[i]] <- rep(if(use_dataset_intercepts) i else 1L, nrow(X[[i]]))
  }

  dl <- list(
    N = sum(vapply(X, nrow, integer(1))),
    K = length(feature_names),
    X = do.call(rbind, X),
    observed = unlist(observed, use.names = FALSE),
    offsets = unlist(offsets, use.names = FALSE),
    weights = unlist(weights, use.names = FALSE),
    J = nrow(components),
    J_intercept = if(use_dataset_intercepts) nrow(components) else 1L,
    dataset = unlist(dataset, use.names = FALSE),
    likelihood = as.integer(components$family_code),
    link = as.integer(components$link_code),
    intercept_id = unlist(intercept_id, use.names = FALSE)
  )

  if(isTRUE(settings$get("optim_hyperparam")) && is.Waiver(model$priors)){
    dl$hs_df <- 1
    dl$hs_df_global <- 1
    dl$hs_df_slab <- 4
    dl$hs_scale_global <- 1
    dl$hs_scale_slab <- 2
  }
  dl
}

#' Write a cmdstanr model output to a specific file
#'
#' @description Write a cmdstanr model output to a specific destination
#'
#' @param mod A supplied cmdstanr model
#' @param dir The model directory where the model chould be written. Should be a
#' character / existing dir.
#'
#' @keywords stan utils
#'
#' @noRd
#'
#' @keywords internal
write_stanmodel <- function( mod, dir = tempdir() ) {
  assertthat::assert_that(
    dir.exists(dir)
  )
  fname <- file.path( dir , paste0("rt_cmdstanr_", digest::digest(mod,"md5")) )
  file_stan <- paste0( fname, ".stan" )
  fileConn <- file( file_stan )
  writeLines( mod , fileConn )
  close(fileConn)
  return(file_stan)
}

#' Fit a cmdstanr model
#'
#' @description This function fits a stan model using the light-weight interface
#' provided by cmdstanr. The code was adapted from McElreath rethinking package.
#'
#' @param model_code A [`character`] pointing to the stan modelling code.
#' @param data A [`list`] with all the parameters required to run the model_code
#' in stan.
#' @param algorithm A [`character`] giving the algorithm to use. Either \code{'sampling'}
#' (Default), \code{'optimize'} or \code{'variational'} for penalized likelihood estimation.
#' @param chains A [`numeric`] indicating the number of chains to use for estimation.
#' @param cores Number of threads for sampling. Default set to \code{'getOption("ibis.nthread")'}.
#' See [ibis_options()].
#' @param threads [`numeric`] giving the number of threads to be run per chain.
#' Has to be specified in accordance with cores.
#' @param iter A [`numeric`] value giving the number of MCMC samples to generate.
#' @param warmup [`numeric`] for the number of warm-up samples for MCMC. Default
#' set to 1/2 of iter.
#' @param control A [`list`] with further control options for stan.
#' @param cpp_options A [`list`] with options for the Cpp compiling.
#' @param force [`logical`] indication whether to force recompile the model
#' (Default: \code{FALSE}).
#' @param path [`character`] indicating a path to be made available to the stan
#' compiler.
#' @param save_warmup A [`logical`] flag whether to save the warmup samples.
#' @param return_stanfit A [`logical`] flag whether to convert sampling output
#' to an [`rstan`] stanfit object. Defaults to \code{FALSE}; the native
#' cmdstanr CmdStanFit object is used otherwise.
#' @param ... Other non-specified parameters.
#'
#' @returns A cmdstanr object by default, or a rstan object when requested and
#' conversion succeeds.
#'
#' @seealso rethinking R package
#'
#'
#' @keywords misc stan
#'
#' @examples
#' \dontrun{
#' stan_file <- tempfile(fileext = ".stan")
#' writeLines(
#'   "parameters { real y; } model { y ~ normal(0, 1); }",
#'   stan_file
#' )
#' fit <- run_stan(stan_file, data = list(), chains = 1, cores = 1, iter = 100)
#' }
#'
#' @export
run_stan <- function( model_code, data = list(),
                      algorithm = "sampling",
                      chains = 4, cores = getOption("ibis.nthread"),
                      threads = 1,
                      iter = 1000, warmup = floor(iter / 2),
                      control = list(adapt_delta = 0.95),
                      cpp_options = list(),
                      force = FALSE,
                      path = base::getwd(),
                      save_warmup = TRUE,
                      return_stanfit = FALSE, ... ) {
  assertthat::assert_that(
    is.numeric(chains), is.numeric(cores),
    is.numeric(iter), is.numeric(warmup),
    is.numeric(threads),
    threads <= cores,
    is.list(data),
    is.list(control), is.list(cpp_options),
    is.logical(save_warmup),
    is.logical(return_stanfit),
    is.logical(force)
  )
  # Check that cmdstanr is available
  check_package("cmdstanr")
  cmdstanr::check_cmdstan_toolchain(quiet = TRUE)

  # Match the algorithm to be used
  algorithm <- match.arg(algorithm, c("sampling", "optimize", "variational"), several.ok = FALSE)

  if( threads > 1 ) cpp_options[['stan_threads']] <- TRUE

  # Check extension
  assertthat::assert_that(
    is.character(model_code),
    assertthat::has_extension(model_code, "stan")
  )

  # Now compile the model
  mod <- cmdstanr::cmdstan_model( model_code,
                                  compile = TRUE,
                                  force_recompile = force,
                                  cpp_options = cpp_options,
                                  include_paths = path
                                  # stanc_options = list("O1") # Can result in substantial speedups!
                                  )

  # Final parameters for sampling
  samp <- iter - warmup
  warm <- warmup

  # pull out any control arguments
  carg_adapt_delta <- 0.95
  if ( !is.null( control[['adapt_delta']] ) )
    carg_adapt_delta <- as.numeric(control[['adapt_delta']])
  carg_max_treedepth <- 11
  if ( !is.null( control[['max_treedepth']] ) )
    carg_max_treedepth <- as.numeric(control[['max_treedepth']])

  if(algorithm == "sampling"){
    # Sample
    if ( threads > 1 ) {
      cmdstanfit <- mod$sample( data = data,
                                chains = chains,
                                parallel_chains = cores,
                                iter_sampling = samp, iter_warmup = warm,
                                adapt_delta = carg_adapt_delta,
                                max_treedepth = carg_max_treedepth,
                                threads_per_chain = threads,
                                save_warmup = save_warmup,
                                ... )
    } else {
      cmdstanfit <- mod$sample( data = data,
                                chains = chains,
                                parallel_chains = cores,
                                iter_sampling = samp , iter_warmup = warm,
                                adapt_delta = carg_adapt_delta,
                                max_treedepth = carg_max_treedepth,
                                save_warmup = save_warmup,
                                ... )
    }
    if(isTRUE(return_stanfit)){
      # Conversion is optional because rstan can lag behind CmdStan CSV metadata
      # changes even when the CmdStan fit itself is valid.
      stanfit <- try(rstan::read_stan_csv( cmdstanfit$output_files() ), silent = TRUE)
      if(inherits(stanfit, "try-error")){
        cli::cli_alert_warning("Could not convert CmdStan output to a {.cls stanfit}; keeping the {.cls CmdStanFit} object.")
        stanfit <- cmdstanfit
      }
    } else {
      stanfit <- cmdstanfit
    }

  } else if(algorithm == "optimize"){
    # Optimize for getting point estimates
    stanfit <- mod$optimize(data = data,
                               #seed = seed, # This could be passed on
                               threads = threads
                               )
  } else if(algorithm == "variational") {
    # Variational for approximating the posterior
    stanfit <- mod$variational(data = data,
                                  # seed = seed,
                                  threads = threads
    )
  }

  return(stanfit)
}

#' Summarize Stan fit objects from rstan or cmdstanr
#'
#' @keywords internal
#' @noRd
stan_fit_summary <- function(obj){
  assertthat::assert_that(
    inherits(obj, "stanfit") || inherits(obj, "CmdStanFit"),
    msg = "Expected a stanfit or CmdStanFit object."
  )
  if(inherits(obj, "stanfit")){
    rstan::summary(obj)$summary |>
      as.data.frame() |>
      tibble::rownames_to_column(var = "parameter") |>
      as.data.frame()
  } else {
    out <- posterior::summarise_draws(obj$draws()) |> as.data.frame()
    names(out)[names(out) == "variable"] <- "parameter"
    out
  }
}

#' Create a posterior prediction from a Stan fit object
#'
#' @description This function does simulates from the posterior of a created
#' stan model, therefore providing a fast and efficient way to project coefficients
#' obtained from Bayesian models to new/novel contexts.
#'
#' @param obj A \code{"stanfit"} object from rstan or a \code{"CmdStanFit"}
#' object from cmdstanr.
#' @param form A [`formula`] object created for the [ibis.iSDM::DistributionModel].
#' @param newdata A [data.frame] with new data to be used for prediction.
#' @param type A [`character`] of whether the linear `predictor` or the `response`
#' is to be summarized.
#' @param family A [`character`] giving the family for simulating linear response
#' values (Default: \code{NULL})
#' @param offset A [vector] with an optionally specified offset.
#' @param draws [numeric] indicating whether a specific number of draws should be taken.
#' @param intercept Optional intercept parameter name or fixed numeric intercept.
#' @param link Optional inverse-link name. Supports \code{"log"}, \code{"logit"},
#' \code{"cloglog"}, and \code{"identity"}.
#' @param feature_names Optional ordered feature names matching the Stan beta
#' vector.
#'
#' @references
#' * [https://medium.com/@alex.pavlakis/making-predictions-from-stan-models-in-r-3e349dfac1ed](https://medium.com/@alex.pavlakis/making-predictions-from-stan-models-in-r-3e349dfac1ed).
#' * The brms R-package.
#'
#' @return A [`data.frame`] of posterior predictions summarized over draws.
#'
#' @examples
#' \dontrun{
#' newdata <- data.frame(temperature = c(0.1, 0.4, 0.8))
#' posterior_predict_stanfit(
#'   obj = stan_fit,
#'   form = observed ~ temperature,
#'   newdata = newdata,
#'   family = "binomial",
#'   type = "response"
#' )
#' }
#'
#' @export
posterior_predict_stanfit <- function(obj, form, newdata, type = "predictor",
                                      family = NULL, offset = NULL, draws = NULL,
                                      intercept = NULL, link = NULL,
                                      feature_names = NULL){
  assertthat::assert_that(
    inherits(obj, "stanfit") || inherits(obj, "CmdStanFit"),
    is.formula(form),
    is.data.frame(newdata),
    is.character(type),
    is.null(family) || is.character(family),
    is.null(link) || is.character(link),
    is.null(draws) || is.numeric(draws),
    is.null(offset) || (length(offset) == nrow(newdata)),
    is.null(intercept) || is.character(intercept) || is.numeric(intercept),
    is.null(feature_names) || is.character(feature_names)
  )
  if(type == "link") type <- "predictor"
  type <- match.arg(type, c("predictor", "response"), several.ok = FALSE)
  # Build model matrix
  # Note: This removes all NA cells from matrix
  A <- stats::model.matrix(object = stats::delete.response(stats::terms(form)),
                           data = newdata)
  assertthat::assert_that(nrow(A)>0, inherits(A, "matrix") || inherits(A, "dgCMatrix"))
  if("(Intercept)" %in% colnames(A)){
    A <- A[, colnames(A) != "(Intercept)", drop = FALSE]
  }
  if(!is.null(feature_names)){
    assertthat::assert_that(
      all(feature_names %in% colnames(A)),
      msg = paste0("Missing Stan prediction features: ",
                   paste(feature_names[!feature_names %in% colnames(A)], collapse = ", "))
    )
    A <- A[, feature_names, drop = FALSE]
  }

  # Draw from the posterior
  if(inherits(obj, "stanfit")) {
    pp <- posterior::as_draws_df(obj)
  } else {
    pp <- posterior::as_draws_df(obj$draws())
  }
  pp <- as.data.frame(pp)
  # Create a subset?
  if (!is.null(draws)) {
    pp <- pp[sample.int(nrow(pp), draws),]
  }

  # Extract beta coefficients in numeric vector order.
  beta_cols <- grep("^beta\\[", colnames(pp), value = TRUE)
  if(length(beta_cols) > 0){
    beta_index <- as.integer(sub("^beta\\[([0-9]+)\\]$", "\\1", beta_cols))
    beta_cols <- beta_cols[order(beta_index)]
  } else if("beta" %in% colnames(pp)) {
    beta_cols <- "beta"
  }
  assertthat::assert_that(
    length(beta_cols) == ncol(A),
    msg = paste0("Stan posterior beta count (", length(beta_cols),
                 ") does not match prediction matrix columns (", ncol(A), ").")
  )
  beta <- as.matrix(pp[, beta_cols, drop = FALSE])

  # Extract the requested intercept. Numeric values are treated as fixed offsets.
  if(is.null(intercept)){
    if("Intercept" %in% colnames(pp)){
      intercept_draws <- pp[["Intercept"]]
    } else if("Intercept[1]" %in% colnames(pp)){
      intercept_draws <- pp[["Intercept[1]"]]
    } else {
      intercept_draws <- rep(0, nrow(pp))
    }
  } else if(is.numeric(intercept)){
    assertthat::assert_that(length(intercept) == 1)
    intercept_draws <- rep(intercept, nrow(pp))
  } else {
    assertthat::assert_that(
      length(intercept) == 1,
      intercept %in% colnames(pp),
      msg = paste0("Requested Stan intercept parameter not found: ", intercept)
    )
    intercept_draws <- pp[[intercept]]
  }

  # Prepare offset if set
  if(!is.null(offset)) {
    # Get only the rows in the A matrix (minus NA)
    offset <- offset[as.numeric(row.names(A))]
  } else { offset <- rep(0, nrow(A) ) }

  # Security checks
  assertthat::assert_that(
    nrow(A)>0, nrow(beta) > 0,
    ncol(beta) == ncol(A),
    is.numeric(offset)
  )

  a <- beta %*% t(A)
  a <- sweep(a, 1, intercept_draws, "+")
  a <- sweep(a, 2, offset, "+")

  if(is.null(link)){
    link <- switch(family,
                   "poisson" = "log",
                   "binomial" = "logit",
                   "log")
  }
  link <- match.arg(link, c("log", "logit", "cloglog", "identity"), several.ok = FALSE)

  if(type == "response"){
    a <- switch(link,
                "log" = exp(a),
                "logit" = ilink(a, link = "logit"),
                "cloglog" = ilink(a, link = "cloglog"),
                "identity" = a)
  }

  # Finally summarize
  preds <- cbind(
    matrixStats::colMeans2(a, na.rm = TRUE),
    matrixStats::colQuantiles(a, probs = c(.05,.5,.95), na.rm = TRUE),
    matrixStats::colSds(a, na.rm = TRUE)
  )

  # ---- #
  # Create output with cellid
  out <- tibble::rowid_to_column(newdata, var = "cellid")["cellid"] |> as.data.frame()
  out$cv <- out$q95 <- out$q50 <- out$q05 <- out$sd <- out$mean <- NA
  out$mean[as.numeric(row.names(A))] <- preds[,1]
  out$sd[as.numeric(row.names(A))] <- preds[,5]
  out$q05[as.numeric(row.names(A))] <- preds[,2]
  out$q50[as.numeric(row.names(A))] <- preds[,3]
  out$q95[as.numeric(row.names(A))] <- preds[,4]
  out$cv[as.numeric(row.names(A))] <- preds[,5] / preds[,1]
  out$cellid <- NULL

  return(out)
}

#' Show the stan code from a trained model
#'
#' @description This helper function shows the code from a trained
#' [DistributionModel] using the [`engine_stan`]. This function is emulated
#' after a similar functionality in the brms R-package.
#' **It only works with models inferred with stan!**
#'
#' @param obj Any prepared object.
#' @param ... not used.
#'
#' @return None.
#'
#' @seealso rstan, cmdstanr, brms
#' @keywords engine
#'
#' @examples
#' \dontrun{
#' stancode(fitted_stan_model)
#' }
#'
#' @name stancode
NULL

#' @rdname stancode
#' @export
methods::setGeneric("stancode",
                    signature = methods::signature("obj"),
                    function(obj, ...) standardGeneric("stancode"))

#' @rdname stancode
#' @export
stancode.DistributionModel <- function(obj, ...) obj$stancode()
