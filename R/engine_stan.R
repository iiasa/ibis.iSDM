#' @include class-engine.R class-distributionmodel.R
NULL

#' Use Stan as engine
#'
#' @description Stan is probabilistic programming language that can be used to
#' specify linear Bayesian species distribution models with one or more point
#' presence-only and presence-absence biodiversity datasets. Multiple datasets
#' are fitted jointly with shared covariate effects and dataset-specific
#' intercepts when requested via the biodiversity dataset settings.
#' **Requires the \code{"cmdstanr"} package to be installed!**
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param chains A positive [`integer`] specifying the number of Markov chains
#'   (Default: \code{4} chains).
#' @param iter A positive [`integer`] specifying the number of iterations for
#'   each chain (including warmup). (Default: \code{2000}).
#' @param warmup A positive [`integer`] specifying the number of warmup
#'   iterations per chain. The default is \code{iter/2}.
#' @param cores If set to NULL take values from specified ibis option
#'   \code{getOption('ibis.nthread')}.
#' @param init Initial values for parameters (Default: \code{'random'}).
#' @param algorithm Mode used to sample from the posterior. Available options are
#'   \code{"sampling"}, \code{"optimize"}, or \code{"variational"}.
#' @param control See \code{"cmdstanr"} for more details on sampler controls.
#' @param type The mode used for creating posterior predictions. Either
#'   \code{"predictor"} or \code{"response"} (Default: \code{"response"}).
#' @param ... Other variables.
#'
#' @returns An [Engine].
#'
#' @seealso rstan, cmdstanr
#' @family engine
#'
#' @examples
#' \dontrun{
#' x <- distribution(background) |> engine_stan(iter = 1000)
#' }
#'
#' @name engine_stan
NULL

#' @rdname engine_stan
#' @export
engine_stan <- function(x,
                        chains = 4,
                        iter = 2000,
                        warmup = floor(iter/2),
                        init = "random",
                        cores = getOption("ibis.nthread"),
                        algorithm = 'sampling',
                        control = list(adapt_delta = 0.95),
                        type = "response",
                        ...) {
  check_package('rstan')
  if(!isNamespaceLoaded("rstan")) {
    attachNamespace("rstan")
    requireNamespace('rstan')
  }
  stan_check_cmd(install = TRUE)
  check_package("cmdstanr")
  assertthat::assert_that(cmdstanr::cmdstan_version() > "2.26.0")

  assertthat::assert_that(
    inherits(x, "BiodiversityDistribution"),
    inherits(x$background, 'sf'),
    is.numeric(chains), is.numeric(iter), is.numeric(warmup),
    is.null(cores) || is.numeric(cores),
    is.character(init) || is.list(init),
    is.null(control) || is.list(control),
    is.character(algorithm),
    msg = 'Input parameters wrongly specified!'
  )
  algorithm <- match.arg(algorithm, c("sampling", "optimize", "variational"), several.ok = FALSE)
  type <- match.arg(type, c("response", "predictor"), several.ok = FALSE)
  if(is.null(cores)) cores <- getOption('ibis.nthread')

  template <- create_background(x)
  if(!is.Waiver(x$predictors) && x$predictors$is_spatial()){
    template <- terra::mask(template, sum(x$predictors$get_data(), na.rm = TRUE))
  }

  if(!is.Waiver(x$engine)) myLog('[Setup]', 'yellow', 'Replacing currently selected engine.')

  eg <- Engine

  eg$set("public", "set_control", function(chains = 4,
                                           iter = 2000,
                                           warmup = floor(iter/2),
                                           init = "random",
                                           cores = NULL,
                                           control = NULL,
                                           type = NULL){
    if(is.null(type)) type <- self$stan_param$type
    self$stan_param <- list(
      chains = chains, iter = iter, warmup = warmup, init = init,
      cores = cores, algorithm = algorithm, control = control, type = type
    )
    invisible()
  }, overwrite = TRUE)

  eg$set("public", "get_equation_latent_spatial", function(){ NULL }, overwrite = TRUE)

  eg$set("public", "setup", function(model, settings = NULL, ...){
    assertthat::assert_that(
      assertthat::has_name(model, 'background'),
      assertthat::has_name(model, 'biodiversity'),
      inherits(settings, 'Settings') || is.null(settings),
      nrow(model$predictors) == terra::ncell(self$get_data('template'))
    )
    assertthat::assert_that(
      is.numeric(self$stan_param$chains),
      is.numeric(self$stan_param$iter),
      is.numeric(self$stan_param$warmup)
    )
    options(mc.cores = self$stan_param$cores)

    if(getOption('ibis.setupmessages', default = TRUE)) {
      myLog('[Estimation]', 'green', 'Building Stan code.')
    }
    model <- stan_prepare_model(model, settings, self$get_data("template"))
    self$set_data("stancode", model$stan$code)

    if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]', 'green', 'Engine setup.')
    model
  }, overwrite = TRUE)

  eg$set("public", "train", function(model, settings, ...){
    if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]', 'green', 'Starting fitting...')

    settings$set('algorithm', self$stan_param$algorithm)
    settings$set('cores', self$stan_param$cores)
    settings$set('chains', self$stan_param$chains)
    settings$set('iter', self$stan_param$iter)
    settings$set('warmup', self$stan_param$warmup)
    settings$set('type', self$stan_param$type)

    dl <- model$stan$data

    fpath_code <- write_stanmodel(self$get_data("stancode"))
    fit_stan <- run_stan(
      model_code = fpath_code,
      data = dl,
      algorithm = settings$get('algorithm'),
      cores = self$stan_param$cores,
      chains = self$stan_param$chains,
      iter = self$stan_param$iter,
      warmup = self$stan_param$warmup,
      control = self$stan_param$control,
      path = getwd(),
      force = TRUE
    )

    if(!settings$get('inference_only')){
      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]', 'green', 'Starting prediction...')
      full <- model$predictors[, unique(c("x", "y", model$stan$feature_names)), drop = FALSE]
      out <- stan_predict(fit_stan, model, settings, full, component = NULL,
                          type = self$stan_param$type)
      prediction <- self$get_data('template')
      prediction <- fill_rasters(post = out, background = prediction)
      prediction <- terra::mask(prediction, model$background)
      try({ rm(out) }, silent = TRUE)
    } else {
      prediction <- NULL
    }

    settings$set('end.time', Sys.time())
    settings$set('stan_components', model$stan$components)

    obj <- DistributionModel

    obj$set("public", "project", function(newdata, offset = NULL,
                                          type = NULL, layer = "mean",
                                          component = NULL){
      assertthat::assert_that(
        nrow(newdata) > 0,
        all(c("x", "y") %in% names(newdata)),
        is.null(offset) || is.numeric(offset),
        is.character(type) || is.null(type)
      )
      fit <- self$get_data("fit_best")
      model <- self$model
      settings <- self$settings
      assertthat::assert_that(inherits(fit, "stanfit") || inherits(fit, "CmdStanFit"))

      newdata_copy <- newdata
      pred_stan <- stan_predict(fit, model, settings, newdata,
                                component = component,
                                offset = offset,
                                type = type)

      if(nrow(newdata_copy) == nrow(model$predictors)){
        prediction <- try({ model_to_background(model) }, silent = TRUE)
      } else {
        prediction <- try({
          terra::rast(newdata_copy[, c("x", "y")],
                      crs = terra::crs(model$background),
                      type = "xyz") |>
            emptyraster()
        }, silent = TRUE)
      }
      prediction <- fill_rasters(pred_stan, prediction)
      if(!is.null(layer) && terra::nlyr(prediction) > 1){
        assertthat::assert_that(
          layer %in% names(prediction),
          msg = paste0("Requested projection layer '", layer,
                       "' not found. Available layers: ",
                       paste(names(prediction), collapse = ", "))
        )
        prediction <- prediction[[layer]]
      }
      return(prediction)
    }, overwrite = TRUE)

    obj$set("public", "partial", function(x.var = NULL, constant = NULL,
                                          variable_length = 100,
                                          values = NULL, newdata = NULL,
                                          plot = FALSE, type = "predictor",
                                          component = NULL){
      fit <- self$get_data('fit_best')
      model <- self$model
      settings <- self$settings
      if(is.null(type)) type <- settings$get("type")
      assertthat::assert_that(
        inherits(fit, 'stanfit') || inherits(fit, "CmdStanFit"),
        is.character(x.var) || is.null(x.var),
        is.numeric(variable_length) && variable_length > 1,
        is.null(newdata) || is.data.frame(newdata),
        is.null(constant) || is.numeric(constant)
      )

      variables <- model$stan$feature_names
      if(is.null(x.var)){
        x.var <- variables
      } else {
        x.var <- match.arg(x.var, variables, several.ok = TRUE)
      }

      if(is.null(newdata)){
        rr <- sapply(model$predictors[, variables, drop = FALSE],
                     function(x) range(x, na.rm = TRUE)) |> as.data.frame()
        if(!is.null(values)) variable_length <- length(values)
        df_partial <- list()
        if(is.null(constant)){
          for(n in names(rr)) df_partial[[n]] <- rep(mean(model$predictors[[n]], na.rm = TRUE), variable_length)
        } else {
          for(n in names(rr)) df_partial[[n]] <- rep(constant, variable_length)
        }
        df_partial <- do.call(cbind, df_partial) |> as.data.frame()
      } else {
        df_partial <- stan_align_newdata(newdata, model)
        df_partial <- dplyr::select(df_partial, dplyr::any_of(variables))
      }

      o <- vector(mode = "list", length = length(x.var))
      names(o) <- x.var
      for(v in x.var){
        df_temp <- df_partial
        if(!is.null(values)){
          df_temp[, v] <- values
        } else {
          df_temp[, v] <- seq(rr[1, v], rr[2, v], length.out = variable_length)
        }

        pred_part <- stan_predict(fit, model, settings, df_temp,
                                  component = component,
                                  type = type)
        pred_part <- cbind("variable" = v,
                           "partial_effect" = df_temp[, v],
                           pred_part)
        o[[v]] <- pred_part
      }

      o <- do.call(what = rbind, args = c(o, make.row.names = FALSE))
      if(plot){
        pm <- ggplot2::ggplot(data = o, ggplot2::aes(x = partial_effect)) +
          ggplot2::theme_classic() +
          ggplot2::geom_ribbon(ggplot2::aes(ymin = mean - sd, ymax = mean + sd), fill = "grey85") +
          ggplot2::geom_line(ggplot2::aes(y = mean)) +
          ggplot2::facet_wrap(. ~ variable, scales = "free") +
          ggplot2::labs(x = "Variable", y = "Partial effect")
        print(pm)
      }
      o
    }, overwrite = TRUE)

    obj$set("public", "spartial", function(x.var, constant = NULL,
                                           newdata = NULL, plot = TRUE,
                                           type = "predictor",
                                           component = NULL, ...){
      fit <- self$get_data('fit_best')
      model <- self$model
      settings <- self$settings
      assertthat::assert_that(
        inherits(fit, 'stanfit') || inherits(fit, "CmdStanFit"),
        is.character(x.var),
        is.null(constant) || is.numeric(constant)
      )
      x.var <- match.arg(x.var, model$stan$feature_names, several.ok = FALSE)

      if(is.null(newdata)){
        df_partial <- model$predictors[, unique(c("x", "y", model$stan$feature_names)), drop = FALSE]
      } else {
        df_partial <- newdata
      }
      df_partial <- stan_align_newdata(df_partial, model)
      if(is.null(constant)){
        for(n in model$stan$feature_names){
          if(n != x.var) df_partial[[n]] <- mean(model$predictors[[n]], na.rm = TRUE)
        }
      } else {
        for(n in model$stan$feature_names){
          if(n != x.var) df_partial[[n]] <- constant
        }
      }

      pred_part <- stan_predict(fit, model, settings, df_partial,
                                component = component,
                                type = type)
      template <- if(nrow(df_partial) == nrow(model$predictors)){
        model_to_background(model)
      } else {
        terra::rast(df_partial[, c("x", "y")],
                    crs = terra::crs(model$background),
                    type = "xyz") |>
          emptyraster()
      }
      template <- fill_rasters(pred_part, template)
      if(plot){
        terra::plot(template[[c("mean", "sd")]], col = ibis_colours$ohsu_palette)
      }
      template
    }, overwrite = TRUE)

    obj$set("public", "has_converged", function(){
      fit <- self$get_data("fit_best")
      if(is.Waiver(fit)) return(FALSE)
      TRUE
    }, overwrite = TRUE)

    obj$set("public", "get_residuals", function(){
      message("Not yet implemented.. :-( ")
      new_waiver()
    }, overwrite = TRUE)

    obj$set("public", "get_coefficients", function(){
      cofs <- self$summary()
      if(nrow(cofs) == 0) return(NULL)
      cofs <- subset(cofs, select = c("parameter", "mean", "sd"))
      names(cofs) <- c("Feature", "Beta", "Sigma")
      int <- grep("Intercept", cofs$Feature, ignore.case = TRUE)
      if(length(int) > 0) cofs <- cofs[-int, ]
      cofs
    }, overwrite = TRUE)

    obj$set("public", "plot_spatial", function(out, plot = TRUE){ NULL }, overwrite = TRUE)

    obj$set("public", "stancode", function(){
      message(self$get_data("sm_code"))
    }, overwrite = TRUE)

    out <- obj$new(name = "STAN-Model")
    out$id <- model$id
    out$model <- model
    out$settings <- settings
    out$fits <- list(
      "fit_best" = fit_stan,
      "prediction" = prediction,
      "sm_code" = self$get_data("stancode")
    )
    out
  }, overwrite = TRUE)

  eg <- eg$new(engine = "STAN-Engine", name = "<STAN>")
  eg$data <- list('template' = template)
  eg$stan_param <- list(
    chains = chains,
    iter = iter,
    warmup = warmup,
    init = init,
    cores = cores,
    algorithm = algorithm,
    control = control,
    type = type
  )

  y <- x$clone(deep = TRUE)
  y$set_engine(eg)
}
