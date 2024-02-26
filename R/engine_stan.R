#' @include class-engine.R class-distributionmodel.R
NULL
#' Use Stan as engine
#'
#' @description Stan is probabilistic programming language that can be used to
#' specify most types of statistical linear and non-linear regression models.
#' Stan provides full Bayesian inference for continuous-variable models
#' through Markov chain Monte Carlo methods such as the No-U-Turn sampler, an
#' adaptive form of Hamiltonian Monte Carlo sampling. Stan code has to be
#' written separately and this function acts as compiler to build the
#' stan-model.
#' **Requires the \code{"cmdstanr"} package to be installed!**
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param chains A positive [`integer`] specifying the number of Markov chains
#' (Default: \code{4} chains).
#' @param iter A positive [`integer`] specifying the number of iterations for
#' each chain (including warmup). (Default: \code{2000}).
#' @param warmup A positive [`integer`] specifying the number of warmup (aka
#' burnin) iterations per chain. If step-size adaptation is on (Default: \code{TRUE}),
#' this also controls the number of iterations for which adaptation is run (and
#' hence these warmup samples should not be used for inference). The number of
#' warmup iterations should be smaller than \code{iter} and the default is \code{iter/2}.
#' @param cores If set to NULL take values from specified ibis option \code{getOption('ibis.nthread')}.
#' @param init Initial values for parameters (Default: \code{'random'}). Can also
#' be specified as [list] (see: \code{"rstan::stan"})
#' @param algorithm Mode used to sample from the posterior. Available options are
#' \code{"sampling"}, \code{"optimize"}, or \code{"variational"}. See \code{"cmdstanr"}
#' package for more details. (Default: \code{"sampling"}).
#' @param control See \code{"rstan::stan"} for more details on specifying the controls.
#' @param type The mode used for creating posterior predictions. Either summarizing
#' the linear \code{"predictor"} or \code{"response"} (Default: \code{"response"}).
#' @param ... Other variables
#'
#' @details By default the posterior is obtained through sampling, however stan
#' also supports approximate inference forms through penalized maximum
#' likelihood estimation (see Carpenter et al. 2017).
#'
#' @note The function \code{obj$stancode()} can be used to print out the
#' stancode of the model.
#'
#' @returns An [Engine].
#'
#' @references
#' * Jonah Gabry and Rok Češnovar (2021). cmdstanr: R Interface to 'CmdStan'.
#' https://mc-stan.org/cmdstanr, https://discourse.mc-stan.org.
#' * Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B., Betancourt, M.,
#' ... & Riddell, A. (2017). Stan: A probabilistic programming language. Journal of
#' statistical software, 76(1), 1-32.
#' * Piironen, J., & Vehtari, A. (2017). Sparsity information and regularization
#' in the horseshoe and other shrinkage priors. Electronic Journal of Statistics, 11(2), 5018-5051.
#'
#' @seealso rstan, cmdstanr
#' @family engine
#'
#' @examples
#' \dontrun{
#' # Add Stan as an engine
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
  # Check whether INLA package is available
  check_package('rstan')
  if(!isNamespaceLoaded("rstan")) { attachNamespace("rstan");requireNamespace('rstan') }
  stan_check_cmd(install = TRUE)
  check_package("cmdstanr")
  assertthat::assert_that( cmdstanr::cmdstan_version()>"2.26.0")

  # assert that arguments are valid
  assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                          inherits(x$background,'sf')
  )
  # Other checks of parameters
  assertthat::assert_that(
    is.numeric(chains), is.numeric(iter), is.numeric(warmup),
    is.null(cores) || is.numeric(cores),
    is.character(init) || is.list(init),
    is.null(control) || is.list(control),
    is.character(algorithm),
    msg = 'Input parameters wrongly specified!'
  )
  # Match algorithm and posterior prediction type
  algorithm <- match.arg(algorithm, c("sampling", "optimize", "variational"),several.ok = FALSE)
  type <- match.arg(type, c("response", "predictor"), several.ok = FALSE)
  if(is.null(cores)) cores <- getOption('ibis.nthread')

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

  # Define new engine object of class
  eg <- Engine

  # Function to respecify the control parameters
  eg$set("public", "set_control", function(chains = 4,
                            iter = 2000,
                            warmup = 500,
                            init = "random",
                            cores = NULL,
                            control = NULL){

    # Overwrite existing
    self$stan_param <- list(
      chains = chains, iter = iter,
      warmup = warmup, init = init,
      cores = cores, algorithm = algorithm,
      control = control
    )

  },overwrite = TRUE)
  # Spatial latent effect
  eg$set("public", "get_equation_latent_spatial", function(){ return( NULL )},overwrite = TRUE)

  # Setup a model
  eg$set("public", "setup", function(model, settings = NULL, ...){
    # Simple security checks
    assertthat::assert_that(
      assertthat::has_name(model, 'background'),
      assertthat::has_name(model, 'biodiversity'),
      inherits(settings,'Settings') || is.null(settings),
      nrow(model$predictors) == ncell(self$get_data('template'))
    )
    # Check that all stan parameters are appropriately set
    assertthat::assert_that(
      is.numeric(self$stan_param$chains),
      is.numeric(self$stan_param$iter),
      is.numeric(self$stan_param$warmup)
    )
    # Set cores
    options(mc.cores = self$stan_param$cores)

    # FIXME: Stan should handle factors directly. For now outsourced to split up
    if(any(model$predictors_types$type=="factor")){
      vf <- model$predictors_types$predictors[model$predictors_types$type=="factor"]
      for(k in vf){
        o <- explode_factor(model$predictors[[k]],name = k)
        model$predictors <- cbind(model$predictors, o)
        model$predictors_names <- c(model$predictors_names, colnames(o))
        model$predictors_types <- rbind(model$predictors_types,
                                        data.frame(predictors = colnames(o), type = "numeric") )
        # Finally remove the original column from the predictor object
        model$predictors[[k]] <- NULL
        model$predictors_names <- model$predictors_names[-which( model$predictors_names == k )]
        model$predictors_types <- subset(model$predictors_types, subset = predictors != k)
        # Explode the columns in the raster object
        model$predictors_object$data <- c(
          model$predictors_object$data,
          explode_factorized_raster(model$predictors_object$data[[k]])
        )
        model$predictors_object$data <- terra::subset(model$predictors_object$data, k, negate = TRUE)
      }
    }

    # Stan procedure - First add integration points to all poipo datasets
    # FIXME: Possibly outsoure this across methods
    for(i in 1:length(model$biodiversity)){

      # If there any factor variables split them per type and explode them
      if(any(model$biodiversity[[i]]$predictors_types$type=="factor")){
        vf <- model$biodiversity[[i]]$predictors_types$predictors[model$biodiversity[[i]]$predictors_types$type=="factor"]
        fv <- model$biodiversity[[i]]$predictors[vf]
        for(k in 1:ncol(fv)){
          o <- explode_factor(fv[,k],name = colnames(fv)[k])
          # Add
          model$biodiversity[[i]]$predictors <- cbind(model$biodiversity[[i]]$predictors, o)
          model$biodiversity[[i]]$predictors_names <- c(model$biodiversity[[i]]$predictors_names, colnames(o))
          model$biodiversity[[i]]$predictors_types <- rbind(model$biodiversity[[i]]$predictors_types,
                                                            data.frame(predictors = colnames(o), type = "numeric") )
          # Finally remove the original column from the predictor object
          model$biodiversity[[i]]$predictors[[colnames(fv)[k]]] <- NULL
          model$biodiversity[[i]]$predictors_names <- model$biodiversity[[i]]$predictors_names[-which( model$biodiversity[[i]]$predictors_names == colnames(fv)[k] )]
          model$biodiversity[[i]]$predictors_types <- subset(model$biodiversity[[i]]$predictors_types, subset = predictors != colnames(fv)[k])
        }
      }

      # Add pseudo-absence points if necessary, by including nearest predictor values for each
      if('poipo' == model$biodiversity[[i]]$type){

        # Get background layer
        bg <- self$get_data("template") # model$engine$get_data('template')
        assertthat::assert_that(!is.na(terra::global(bg, "min", na.rm = TRUE)[,1] ))

        # Add pseudo-absence points
        presabs <- add_pseudoabsence(df = model$biodiversity[[i]]$observations,
                                     field_occurrence = 'observed',
                                     template = bg,
                                     settings = model$biodiversity[[i]]$pseudoabsence_settings)
        if(inherits(presabs, 'sf')) presabs <- presabs |> sf::st_drop_geometry()
        # Sample environmental points for absence only points
        abs <- subset(presabs, observed == 0)
        # Re-extract environmental information for absence points
        envs <- get_rastervalue(coords = abs[,c('x','y')],
                                env = model$predictors_object$get_data(df = FALSE),
                                rm.na = FALSE)
        if(assertthat::has_name(model$biodiversity[[i]]$predictors, "Intercept")){ envs$Intercept <- 1}

        # Format out
        df <- rbind(model$biodiversity[[i]]$predictors[,c('x','y','Intercept', model$biodiversity[[i]]$predictors_names)],
                    envs[,c('x','y','Intercept', model$biodiversity[[i]]$predictors_names)] )
        any_missing <- which(apply(df, 1, function(x) any(is.na(x))))
        if(length(any_missing)>0){
          presabs <- presabs[-any_missing,] # This works as they are in the same order
          model$biodiversity[[i]]$expect <- model$biodiversity[[i]]$expect[-any_missing]
        }
        df <- subset(df, stats::complete.cases(df))
        assertthat::assert_that(nrow(presabs) == nrow(df))

        # Check that expect matches
        if(length(model$biodiversity[[1]]$expect)!=nrow(df)){
          # Fill the absences with 1 as multiplier. This works since absences follow the presences
          model$biodiversity[[1]]$expect <- c( model$biodiversity[[1]]$expect,
                                               rep(1, nrow(presabs)-length(model$biodiversity[[1]]$expect) ))
        }

        # Overwrite observation data
        model$biodiversity[[i]]$observations <- presabs

        # Preprocessing security checks
        assertthat::assert_that( all( model$biodiversity[[i]]$observations[['observed']] >= 0 ),
                                 any(!is.na(presabs[['observed']])),
                                 length(model$biodiversity[[i]]$expect)==nrow(model$biodiversity[[i]]$observations),
                                 nrow(df) == nrow(model$biodiversity[[i]]$observations)
        )

        # Add offset if existent
        if(!is.Waiver(model$offset)) {
          # Respecify offset if not set
          of <- model$offset; of[, "spatial_offset" ] <- ifelse(is.na(of[, "spatial_offset" ]), 1, of[, "spatial_offset"])
          of1 <- get_ngbvalue(coords = model$biodiversity[[i]]$observations[,c("x","y")],
                              env =  of,
                              longlat = terra::is.lonlat(bg),
                              field_space = c('x','y')
          )
          df[["spatial_offset"]] <- of1
        }

        # Define expectation as very small vector following Renner et al.
        w <- ppm_weights(df = df,
                         pa = model$biodiversity[[i]]$observations[['observed']],
                         bg = bg,
                         weight = 1e-6
        )
        df$w <- w * (1/model$biodiversity[[i]]$expect) # Also add as column

        model$biodiversity[[i]]$predictors <- df
        model$biodiversity[[i]]$expect <- df$w
      } else {
        # calculating the case weights (equal weights)
        # the order of weights should be the same as presences and backgrounds in the training data
        prNum <- as.numeric(table(model$biodiversity[[i]]$observations[['observed']])["1"]) # number of presences
        bgNum <- as.numeric(table(model$biodiversity[[i]]$observations[['observed']])["0"]) # number of backgrounds
        w <- ifelse(model$biodiversity[[i]]$observations[['observed']] == 1, 1, prNum / bgNum)
        model$biodiversity[[i]]$expect <- w * model$biodiversity[[i]]$expect # Multiply with provided weights
      }
    }
    # --- #
    # Messenger
    if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','green','Building stan code.')
    sm_code <- vector("list", 7)
    names(sm_code) <- c("functions","data","transformed_data","parameters","transformed_parameters",
                        "model","generated_quantities")

    # Has intercept?
    has_intercept <- attr(stats::terms.formula(model$biodiversity[[1]]$equation), "intercept") == 1
    # Family and link function
    fam <- model$biodiversity[[1]]$family
    li <- model$biodiversity[[1]]$link

    # Any spatial or other functions needed?
    if(!is.null(self$get_equation_latent_spatial())){
      # Stan functions for CAR and GP models
      ir <- readLines( system.file("inst/stanfiles/spatial_functions.stan",package = "ibis.iSDM",mustWork = TRUE) )
      assertthat::assert_that(length(ir)>0)
      for(i in ir) sm_code$functions <- append(sm_code$functions, i)
    }

    # Load all the data parameters
    ir <- readLines( system.file("stanfiles/data_parameters.stan",package = "ibis.iSDM",mustWork = TRUE) )
    assertthat::assert_that(length(ir)>0)
    for(i in ir) sm_code$data <- append(sm_code$data, i)

    # Append prior to transformed parameters
    sm_code$transformed_parameters <- append(sm_code$transformed_parameters,"
                            // Prior contribution to log posterior
                            real lprior = 0;")

    # Equation has overall intercept
    if(has_intercept){
      # Add data
      sm_code$transformed_data <- append(sm_code$transformed_data,"
                                              int Kc = K - 1;
                                              matrix[N, Kc] Xc;  // centered version of X without an intercept
                                              vector[Kc] means_X;  // column means of X before centering
                                              for (i in 2:K) {
                                                means_X[i - 1] = mean(X[, i]);
                                                Xc[, i - 1] = X[, i] - means_X[i - 1];
                                              }
                                             ")
      # Add population level effects for rest of coefficients
      sm_code$parameters <- append(sm_code$parameters,"
                                          vector[Kc] b;  // population-level effects
                                          real Intercept;  // temporary intercept for centered predictors
                                       ")
      # add a prior on the intercept
      sm_code$transformed_parameters <- append(sm_code$transformed_parameters,
                                               paste0("
                                    // priors including constants
                                    lprior += student_t_lpdf(Intercept | 3, ",
                                    ifelse(fam == "poisson", -2, 0), # Adapted student prior for poisson
                                    ", 2.5);
                                  ")
      )
      # Generate actual population-level intercept
      sm_code$generated_quantities <- append(sm_code$generated_quantities,"
                                                // actual population-level intercept
                                                real b_Intercept = Intercept - dot_product(means_X, b);
                                                 ")
    }

    # Transformed parameters
    # Add (gaussian) priors to model likelihood if set
    if((!is.Waiver(model$priors) || settings$get(what='optim_hyperparam') == FALSE)){
      # If no intercept is specified, add beta
      if(has_intercept){
        # Parameters
        sm_code$parameters <- append(sm_code$parameters, "
                                         vector[Kc] beta;")
      } else {
        sm_code$parameters <- append(sm_code$parameters, "vector[K] beta;")
      }

      # Add priors for each variable for which it is set to the model
      sm_code$transformed_parameters <- append(sm_code$transformed_parameters, "// beta priors including constants")
      # Now add for each one a normal effect
      for(i in 1:length(model$predictors_names)){
        if(!is.Waiver(model$priors)){
          if(model$predictors_names[i] %in% model$priors$varnames()) {
            # Get prior estimats
            pp <- model$priors$get(model$predictors_names[i])
            sm_code$transformed_parameters <- append(sm_code$transformed_parameters, paste0(
              "lprior += normal_lpdf(beta[",i,"] | ",pp[1],", ",pp[2],");"
            ))
          } else {
            # Default gaussian prior
            sm_code$transformed_parameters <- append(sm_code$transformed_parameters,
                                                     paste0("lprior += normal_lpdf(beta[",i,"] | 0, 2);"))
          }
        } else {
          # Default gaussian prior
          sm_code$transformed_parameters <- append(sm_code$transformed_parameters,
                                                   paste0("lprior += normal_lpdf(beta[",i,"] | 0, 2);"))
        }
      }
    } else
      if( settings$get(what='optim_hyperparam') == TRUE ){
        if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','green','Adding regularized Bayesian priors.')
        # Add regularized horseshoe prior
        # See brms::horseshoe
        ir <- readLines( system.file("stanfiles/prior_functions.stan",package = "ibis.iSDM",mustWork = TRUE) )
        assertthat::assert_that(length(ir)>0)
        for(i in ir) sm_code$functions <- append(sm_code$functions, i)

        sm_code$data <- append(sm_code$data,"
                                    // data for the horseshoe prior
                                    real<lower=0> hs_df;  // local degrees of freedom
                                    real<lower=0> hs_df_global;  // global degrees of freedom
                                    real<lower=0> hs_df_slab;  // slab degrees of freedom
                                    real<lower=0> hs_scale_global;  // global prior scale
                                    real<lower=0> hs_scale_slab;  // slab prior scale"
        )

        # Parameters and transformed parameters
        sm_code$parameters <- append(sm_code$parameters,"
            // local parameters for horseshoe prior
            vector[K] zb;
            vector<lower=0>[K] hs_local;
            // horseshoe shrinkage parameters
            real<lower=0> hs_global;  // global shrinkage parameters
            real<lower=0> hs_slab;  // slab regularization parameter
                                         ")
        sm_code$transformed_parameters <- append(sm_code$transformed_parameters,"
              vector[K] beta;  // population-level effects
              // compute actual regression coefficients
              beta = horseshoe(zb, hs_local, hs_global, hs_scale_slab^2 * hs_slab);
            ")

        # Finally add priors to model
        sm_code$model <- append(sm_code$model, "
            // priors including constants
            target += std_normal_lpdf(zb);
            target += student_t_lpdf(hs_local | hs_df, 0, 1)
            - rows(hs_local) * log(0.5);
            target += student_t_lpdf(hs_global | hs_df_global, 0, hs_scale_global)
            - 1 * log(0.5);
            target += inv_gamma_lpdf(hs_slab | 0.5 * hs_df_slab, 0.5 * hs_df_slab);
            ")
      }

    # Now add the model depending on the type
    if(length(model$biodiversity)>1){
      # For integrated model
      stop("TBD")

    } else if(model$biodiversity[[1]]$type == "poipo" && model$biodiversity[[1]]$family == "poisson"){
      # For poisson process model add likelihood
      if(has_intercept){
        ir <- readLines( system.file("stanfiles/poipo_ll_poisson_intercept.stan",package = "ibis.iSDM",mustWork = TRUE))
      } else {
        ir <- readLines( system.file("stanfiles/poipo_ll_poisson.stan",package = "ibis.iSDM",mustWork = TRUE) )
      }
      assertthat::assert_that(length(ir)>0)
      for(i in ir) sm_code$model <- append(sm_code$model, i)

    } else if(model$biodiversity[[1]]$type == "poipa" && model$biodiversity[[1]]$family == "binomial"){
      # For logistic regression
      if(has_intercept){
        ir <- readLines( system.file("stanfiles/poipa_ll_bernoulli_intercept.stan",package = "ibis.iSDM",mustWork = TRUE) )
      } else {
        ir <- readLines( system.file("stanfiles/poipa_ll_bernoulli.stan",package = "ibis.iSDM",mustWork = TRUE) )
      }
      assertthat::assert_that(length(ir)>0)
      for(i in ir) sm_code$model <- append(sm_code$model, i)
    } else {
      # Else
      stop("Model as of now not implemented for Stan!")
    }
    # Append prior contributions to model
    sm_code$model <- append(sm_code$model, "
                                // Prior contributions
                                target += lprior;")

    # Wrap list entries in model code and save in model object
    self$set_data("stancode", wrap_stanmodel(sm_code))

    # --- #
    # Messenger
    if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','green','Engine setup.')

    # Return modified model object
    return(model)
  },overwrite = TRUE)

  eg$set("public", "train", function(model, settings, ...){
    # Messenger
    if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','green','Starting fitting...')

    # Define an algorithm for MCMC sampling
    # an be "sampling" for MCMC (the default), "optimize" for optimization,
    # "variational" for variational inference with independent normal distributions,
    settings$set('algorithm', self$stan_param$algorithm)
    settings$set('cores', self$stan_param$cores)
    settings$set('chains', self$stan_param$chains)
    settings$set('iter', self$stan_param$iter)
    settings$set('warmup', self$stan_param$warmup)
    settings$set('type', self$stan_param$type)

    # --- #
    # Collect data for stan modelling
    if(length(model$biodiversity)>1){
      stop("done")
    } else {
      has_intercept <- attr(stats::terms(model$biodiversity[[1]]$equation), "intercept")
      # Format data list
      if(has_intercept == 1){
        pn <- c("Intercept",model$biodiversity[[1]]$predictors_names)
      } else { pn <- model$biodiversity[[1]]$predictors_names }

      dl <- list(
        N = nrow( model$biodiversity[[1]]$observations),
        observed = model$biodiversity[[1]]$observations[["observed"]],
        X = as.matrix( model$biodiversity[[1]]$predictors[, pn] ),
        K = length( pn ),
        offsets = log(model$biodiversity[[1]]$expect), # Notice that exposure is log-transformed here!
        has_intercept = attr(stats::terms(model$biodiversity[[1]]$equation), "intercept"),
        has_spatial = ifelse(is.null(self$get_equation_latent_spatial()), 0, 1),
        # Horseshoe prior default parameters
        # FIXME: Allow passing this one via a parameter
        hs_df = 1,
        hs_df_global = 1, hs_df_slab = 4,
        hs_scale_global = 1, hs_scale_slab = 2
      )
      # If any additional offset is set, simply to the existing one in sum
      # This works as log(2) + log(5) == log(2*5)
      if(!is.Waiver(model$offset)) dl$offsets <- dl$offsets + model$biodiversity[[1]]$offset[,"spatial_offset"]
    }

    # Model estimation
    # ---- #
    # Fitting
    fpath_code <- write_stanmodel( self$get_data("stancode") )
    fit_stan <- run_stan(
      model_code = fpath_code,
      data = dl,
      algorithm = settings$get('algorithm'),
      cores = self$stan_param$cores,
      chains = self$stan_param$chains,
      iter = self$stan_param$iter,
      warmup = self$stan_param$warmup,
      path = getwd(),
      force = TRUE # Force recompile
    )

    # Prediction
    if(!settings$get('inference_only')){
      # Messager
      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','green','Starting prediction...')

      # Prepare prediction dataset
      prediction <- self$get_data('template') # Get output raster and new data
      # Full data for prediction
      full <- subset(model$predictors, select = c('x','y',model$predictors_names))

      # Clamp?
      if( settings$get("clamp") ) full <- clamp_predictions(model, full)

      if(has_intercept==1) full$Intercept <- 1

      # If poipo, add w to prediction container
      bd_poipo <- sapply(model$biodiversity, function(x) x$type) == "poipo"
      if(any(bd_poipo)){
        # FIXME: Bit hackish. See if works for other projections
        full$w <- unique(model$biodiversity[[which(bd_poipo)]]$expect)[2] # Absence location being second unique value
      }
      # Add offset if set
      if(!is.Waiver(model$offset)) {
        # Offsets are simply added linearly (albeit transformed)
        if(utils::hasName(full,"w")) full$w <- full$w + model$offset[,"spatial_offset"] else full$w <- model$offset[,"spatial_offset"]
      }
      suppressWarnings(
        full <- sp::SpatialPointsDataFrame(coords = full[,c("x","y")],
                                           data = full,
                                           proj4string = sp::CRS(sp::proj4string(methods::as(model$background, "Spatial")))
        )
      )
      full <- methods::as(full, 'SpatialPixelsDataFrame')

      # Set target variables to bias_value for prediction if specified
      if(!is.Waiver(settings$get('bias_variable'))){
        for(i in 1:length(settings$get('bias_variable'))){
          if(settings$get('bias_variable')[i] %notin% names(full)) next()
          full[[settings$get('bias_variable')[i]]] <- settings$get('bias_value')[i]
        }
      }

      # For Integrated model, follow poisson
      fam <- ifelse(length(model$biodiversity)>1, "poisson", model$biodiversity[[1]]$family)

      # Do the prediction by sampling from the posterior
      out <- posterior_predict_stanfit(obj = fit_stan,
                                       form = to_formula(paste0("observed ~ ",
                                                                ifelse(has_intercept==1, "Intercept + ", ""),
                                                                paste(model$biodiversity[[1]]$predictors_names,collapse = " + "))),
                                       newdata = full@data,
                                       offset = (full$w),
                                       family = fam, # Family
                                       mode = self$stan_param$type # Type
      )

      # Convert full to raster
      prediction <- terra::rast(full)
      # Fill output
      prediction <- fill_rasters(post = out, background = prediction)
      prediction <- terra::mask(prediction, model$background) # Mask with background
      # plot(prediction$mean, col = ibis.iSDM:::ibis_colours$sdm_colour)
      try({ rm(out) })
    } else {
      prediction <- NULL
    }

    # Compute end of computation time
    settings$set('end.time', Sys.time())

    # Definition of STAN Model object ----
    obj <- DistributionModel # Make a copy to set new functions

    # Project function
    obj$set("public", "project", function(newdata, offset = NULL,
                                          type = NULL, layer = "mean"){
      assertthat::assert_that(
        nrow(newdata) > 0,
        all( c("x", "y") %in% names(newdata) ),
        is.null(offset) || is.numeric(offset),
        is.character(type) || is.null(type)
      )
      # Check that fitted model exists
      obj <- self$get_data("fit_best")
      model <- self$model
      settings <- self$settings
      if(is.null(type)) type <- settings$get("type")
      assertthat::assert_that(inherits(obj, "stanfit"),
                              all(model$predictors_names %in% colnames(newdata)))

      # Clamp?
      if( settings$get("clamp") ) newdata <- clamp_predictions(model, newdata)

      # Set target variables to bias_value for prediction if specified
      if(!is.Waiver(settings$get('bias_variable'))){
        for(i in 1:length(settings$get('bias_variable'))){
          if(settings$get('bias_variable')[i] %notin% names(newdata)) next()
          newdata[[settings$get('bias_variable')[i]]] <- settings$get('bias_value')[i]
        }
      }

      # For Integrated model, follow poisson
      fam <- ifelse(length(model$biodiversity)>1, "poisson", model$biodiversity[[1]]$family)

      # Build prediction stack
      full <- subset(newdata, select = c('x','y', model$predictors_names))

      # If poipo, add w to prediction container
      bd_poipo <- sapply(model$biodiversity, function(x) x$type) == "poipo"
      if(any(bd_poipo)){
        # FIXME: Bit hackish. See if works for other projections
        full$w <- unique(model$biodiversity[[which(bd_poipo)]]$expect)[2] # Absence location being second unique value
      }
      # Add offset if set
      if(!is.null(offset)) {
        # Offsets are simply added linearly (albeit transformed)
        if(utils::hasName(full,"w")) full$w <- full$w + offset else full$w <- offset
      }
      suppressWarnings(
        full <- sp::SpatialPointsDataFrame(coords = full[,c("x","y")],
                                           data = full,
                                           proj4string = sp::CRS(sp::proj4string(methods::as(model$background, "Spatial")))
        )
      )
      full <- methods::as(full, 'SpatialPixelsDataFrame')

      # Do the prediction by sampling from the posterior
      pred_stan <- posterior_predict_stanfit(obj = obj,
                                             form = to_formula(paste0("observed ~ ", paste(model$predictors_names,collapse = " + "))),
                                             newdata = full@data,
                                             offset = (full$w),
                                             family = fam,
                                             mode = type # Linear predictor
      )

      # Fill output with summaries of the posterior
      prediction <- try({emptyraster( self$model$predictors_object$get_data()[[1]] )},silent = TRUE) # Background
      if(inherits(prediction, "try-error")){
        prediction <- terra::rast(self$model$predictors[,c("x", "y")], crs = terra::crs(model$background),type = "xyz") |>
          emptyraster()
      }
      prediction <- fill_rasters(pred_stan, prediction)

      return(prediction)

    },overwrite = TRUE)
    # Partial effect
    obj$set("public", "partial", function(x.var = NULL, constant = NULL, variable_length = 100,
                           values = NULL, newdata = NULL, plot = FALSE, type = "predictor"){
      # Get model and intercept if present
      mod <- self$get_data('fit_best')
      model <- self$model
      has_intercept <- attr(stats::terms(model$biodiversity[[1]]$equation), "intercept")
      if(is.null(type)) type <- self$settings$get("type")
      assertthat::assert_that(inherits(mod,'stanfit'),
                              is.character(x.var) || is.null(x.var),
                              is.numeric(variable_length) && variable_length > 1,
                              is.null(newdata) || is.data.frame(newdata),
                              is.null(constant) || is.numeric(constant)
      )

      # get variable names
      variables <- model$predictors_names

      # Match x.var to argument
      if(is.null(x.var)) {
        x.var <- variables
      } else {
        x.var <- match.arg(x.var, variables, several.ok = TRUE)
      }

      if(is.null(newdata)){

        # get range of data
        rr <- sapply(model$predictors[, variables], function(x) range(x, na.rm = TRUE)) |> as.data.frame()

        df_partial <- list()

        if(!is.null(values)){ variable_length <- length(values) }

        # Add all others as constant
        if(is.null(constant)){
          for(n in names(rr)) df_partial[[n]] <- rep( mean(model$predictors[[n]], na.rm = TRUE), variable_length )
        } else {
          for(n in names(rr)) df_partial[[n]] <- rep( constant, variable_length )
        }
        df_partial <- do.call(cbind, df_partial) |> as.data.frame()
      } else {
        df_partial <- dplyr::select(newdata, dplyr::any_of(model$predictors_names))
      }

      # For Integrated model, follow poisson
      fam <- ifelse(length(model$biodiversity)>1, "poisson", model$biodiversity[[1]]$family)

      # Add intercept if present
      if(has_intercept == 1) df_partial$Intercept <- 1
      # If poipo, add w to prediction container
      bd_poipo <- sapply(model$biodiversity, function(x) x$type) == "poipo"
      if(any(bd_poipo)){
        # FIXME: Bit hackish. See if works for other projections
        df_partial$w <- unique(model$biodiversity[[which(bd_poipo)]]$expect)[2] # Absence location being second unique value
      }
      # Add offset if set
      if(!is.Waiver(model$offset)) {
        # Offsets are simply added linearly (albeit transformed).
        # FIXME: Taken here as average. To be re-evaluated as use case develops!
        if(utils::hasName(df_partial,"w")) df_partial$w <- df_partial$w + mean(model$offset,na.rm = TRUE) else df_partial$w <- mean(model$offset,na.rm = TRUE)
      }

      # create list to store results
      o <- vector(mode = "list", length = length(x.var))
      names(o) <- x.var

      # loop through x.var
      for(v in x.var) {

        df_temp <- df_partial

        if(!is.null(values)){
          df_temp[,v ] <- values
        } else {
          df_temp[, v] <- seq(rr[1, v], rr[2, v], length.out = variable_length)
        }

        # Simulate from the posterior
        pred_part <- posterior_predict_stanfit(obj = mod,
                                               form = to_formula(paste0("observed ~ ",
                                                                        ifelse(has_intercept==1, "Intercept +", ""),
                                                                        paste(model$predictors_names,collapse = " + "))),
                                               newdata = df_temp,
                                               offset = df_temp$w,
                                               family = fam,
                                               mode = type) # Linear predictor

        # FIXME: Something wrong here I guess
        # Also attach the partial variable
        pred_part <- cbind("variable" = v, "partial_effect" = df_temp[, v],
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
      return(o) # Return the partial data
    },overwrite = TRUE)

    # Spatial partial effect plots
    obj$set("public", "spartial", function(x.var, constant = NULL, newdata = NULL,
                            plot = TRUE,type = "predictor", ...){
      # Get model object and check that everything is in order
      mod <- self$get_data('fit_best')
      model <- self$model
      has_intercept <- attr(stats::terms(model$biodiversity[[1]]$equation), "intercept")
      assertthat::assert_that(inherits(mod,'stanfit'),
                              'model' %in% names(self),
                              is.character(x.var),
                              is.null(constant) || is.numeric(constant)
      )

      # Match variable name
      x.var <- match.arg(x.var, model$predictors_names, several.ok = FALSE)

      # Make spatial container for prediction
      suppressWarnings(
        df_partial <- sp::SpatialPointsDataFrame(coords = model$predictors[,c('x', 'y')],
                                                 data = model$predictors[, names(model$predictors) %notin% c('x','y')],
                                                 proj4string = sp::CRS(sp::proj4string(methods::as(model$background, "Spatial")))
        )
      )
      df_partial <- methods::as(df_partial, 'SpatialPixelsDataFrame')

      # Add all others as constant
      if(is.null(constant)){
        for(n in names(df_partial)) if(n != x.var) df_partial[[n]] <- mean(model$predictors[[n]], na.rm = TRUE)
      } else {
        for(n in names(df_partial)) if(n != x.var) df_partial[[n]] <- constant
      }
      bd_poipo <- sapply(model$biodiversity, function(x) x$type) == "poipo"
      if(any(bd_poipo)){
        # FIXME: Bit hackish. See if works for other projections
        df_partial$w <- unique(model$biodiversity[[which(bd_poipo)]]$expect)[2] # Absence location being second unique value
      }

      # For Integrated model, follow poisson
      fam <- ifelse(length(model$biodiversity)>1, "poisson", model$biodiversity[[1]]$family)

      # Check intercept
      if(has_intercept==1) df_partial$Intercept <- 1
      # If poipo, add w to prediction container
      bd_poipo <- sapply(model$biodiversity, function(x) x$type) == "poipo"
      if(any(bd_poipo)){
        # FIXME: Bit hackish. See if works for other projections
        df_partial$w <- unique(model$biodiversity[[which(bd_poipo)]]$expect)[2] # Absence location being second unique value
      } else { df_partial$w <- NULL }

      # Simulate from the posterior
      pred_part <- posterior_predict_stanfit(obj = mod,
                                             form = to_formula(paste0("observed ~ Intercept + ",
                                                                      paste(model$predictors_names,collapse = " + "))),
                                             newdata = df_partial@data,
                                             offset = df_partial$w,
                                             family = fam,
                                             mode = type # Linear predictor
      )

      # Get container
      template <- model_to_background(model)
      template <- fill_rasters(pred_part, template)

      # Do plot and return result
      if(plot){
        terra::plot(template[[c("mean","sd")]],
                    col = ibis_colours$ohsu_palette)
      }
      return(template)
    },overwrite = TRUE)

    # Model convergence check
    obj$set("public", "has_converged", function(){
      fit <- self$get_data("fit_best")
      if(is.Waiver(fit)) return(FALSE)
      return(TRUE)
    },overwrite = TRUE)

    # Residual function
    obj$set("public", "get_residuals", function(){
      # Get best object
      message("Not yet implemented.. :-( ")
      new_waiver()
    },overwrite = TRUE)

    # Get coefficients
    obj$set("public", "get_coefficients", function(){
      # Returns a vector of the coefficients with direction/importance
      cofs <- self$summary()
      if(nrow(cofs)==0) return(NULL)
      cofs <- subset(cofs, select = c("parameter", "mean", "sd"))
      names(cofs) <- c("Feature", "Beta", "Sigma")
      # Remove intercept(s)
      int <- grep("Intercept",cofs$Feature,ignore.case = TRUE)
      if(length(int)>0) cofs <- cofs[-int,]
      return(cofs)
    },overwrite = TRUE)

    # Spatial latent effect
    obj$set("public", "plot_spatial", function(out, plot = TRUE){return(NULL)},overwrite = TRUE)

    # Custom function to show stan code
    obj$set("public", "stancode", function(){
      message(
        self$get_data("sm_code")
      )
    },overwrite = TRUE)

    out <- obj$new(name = "STAN-Model")
    out$id = model$id
    out$model = model
    out$settings = settings
    out$fits = list(
      "fit_best" = fit_stan,
      "prediction" = prediction,
      "sm_code" = self$get_data("stancode")
    )
    # Return
    return(out)
  },overwrite = TRUE) # End of Train

  # Define engine object
  eg <- eg$new(engine = "STAN-Engine", name = "<STAN>")
  eg$data <- list(
    'template' = template
  )
  # Stan options
  eg$stan_param = list(
    chains = chains, iter = iter,
    warmup = warmup, init = init,
    cores = cores, algorithm = algorithm,
    control = control,
    type = type
  )

  # Set engine in distribution object
  y <- x$clone(deep = TRUE)
  return( y$set_engine(eg) )
}
