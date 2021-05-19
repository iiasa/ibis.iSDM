#' @include utils.R bdproto-biodiversitydistribution.R utils-spatial.R
NULL

#' Train the model from a given engine
#'
#' Train a [distribution()] model with the specified engine.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object).
#' @param runname A [`character`] name of the trained run
#' @param rm_corPred Remove highly correlated predictors. Default is True
#' @param varsel Perform a variable selection on the set of predictors
#' @param inference_only Fit model only without spatial prediction (Default: FALSE)
#' @param ... further arguments passed on.
#'
#' @details
#'
#' @return A distribution prediction object
#' @name train
#' @exportMethod train
#' @aliases train, train-method
#' @export
NULL

#' @name train
#' @rdname train
#' @exportMethod train
#' @export
methods::setGeneric(
  "train",
  signature = methods::signature("x", "runname","rm_corPred","varsel"),
  function(x, runname, rm_corPred = FALSE, varsel = FALSE, inference_only = FALSE,...) standardGeneric("train"))

#' @name train
#' @rdname train
#' @usage \S4method{train}{BiodiversityDistribution}(x)
methods::setMethod(
  "train",
  methods::signature(x = "BiodiversityDistribution", runname = "character"),
  function(x, runname, rm_corPred = FALSE, varsel = FALSE, inference_only = FALSE, ...) {
    # Make load checks
    assertthat::assert_that(
      inherits(x, "BiodiversityDistribution"),
      is.character(runname),
      is.logical(rm_corPred),
      is.logical(inference_only)
    )
    # Now make checks on completeness of the object
    assertthat::assert_that(!is.Waiver(x$engine),
                            msg = 'No engine set for training the distribution model.')
    assertthat::assert_that( x$show_biodiversity_length() > 0,
                             msg = 'No biodiversity data specified.')
    assertthat::assert_that('observed' %notin% x$get_predictor_names(), msg = 'observed is not an allowed predictor name.' )

    # Check whether prior objects match the used engine, otherwise raise warning
    # if( unique(x$priors$classes()) == 'GDBPrior' && x$engine$name == '<INLA>' ) warning('Priors for wrong models specified, which will be ignored in the inference.')

    # --- #
    #### Defining model objects ----
    # Set model object for fitting
    model <- list()

    # Set model name
    model[['runname']] <- runname

    # Specify a unique id for the run
    model[['id']] <- new_id()

    # Save the background
    model[['background']] <- x$background

    # Get overall Predictor data
    if(is.Waiver(x$get_predictor_names())) {
      # Dummy covariate of background raster
      dummy <- as.data.frame( raster::raster(extent(x$background),nrow=100,ncol=100,val=1), xy = TRUE );names(dummy)[3] <- 'dummy'
      model[['predictors']] <- dummy
      model[['predictors_names']] <- 'dummy'
      model[['predictors_types']] <- data.frame(predictors = 'dummy', type = 'numeric')
    } else {
      # Convert Predictors to data.frame
      model[['predictors']] <- x$predictors$get_data(df = TRUE, na.rm = FALSE)
      # Also set predictor names
      model[['predictors_names']] <- x$get_predictor_names()
      # Try to guess predictor types
      # FIXME: Allow option to determine this directly during add_predictors
      lu <- apply(model[['predictors']][model[['predictors_names']]],
                  2, function(x) length(unique(x[])))
      model[['predictors_types']] <- data.frame(predictors = names(lu), type = ifelse(lu < 50,'factor','numeric') )
      rm(lu)
    }

    # Get biodiversity data
    model[['biodiversity']] <- list()
    # Specify list of ids
    biodiversity_ids <- as.character( x$biodiversity$get_ids() )
    for(id in biodiversity_ids) {
      model[['biodiversity']][[id]][['name']]         <- x$biodiversity$data[[id]]$name # Name of the species
      model[['biodiversity']][[id]][['observations']] <- x$biodiversity$get_data(id) # Observational data
      model[['biodiversity']][[id]][['type']]         <- x$biodiversity$get_types(short = TRUE)[[id]] # Type
      model[['biodiversity']][[id]][['family']]       <- x$biodiversity$get_families()[[id]] # Family
      model[['biodiversity']][[id]][['equation']]     <- x$biodiversity$get_equations()[[id]]
      # --- #
      # Rename observation column to 'observed'. Needs to be consistent for INLA
      # FIXME: try and not use dplyr as dependency (although it is probably loaded already)
      model$biodiversity[[id]]$observations <- model$biodiversity[[id]]$observations %>% dplyr::rename('observed' = x$biodiversity$get_columns_occ()[[id]])
      names(model$biodiversity[[id]]$observations) <- tolower(names(model$biodiversity[[id]]$observations)) # Also generally transfer everything to lower case
      # FIXME: For polygons this won't work. Ideally switch to WKT as default in future
      model$biodiversity[[id]]$observations <- model$biodiversity[[id]]$observations[,c('observed','x','y')] # Get only observed column and coordinates

      # Now extract coordinates and extract estimates
      env <- get_ngbvalue(
        coords = x$biodiversity$get_coordinates(id),
        env = model[['predictors']],
        field_space = c('x','y'),
        longlat = raster::isLonLat(x$background)
      )
      # Remove missing values as several engines can't deal with those easily
      # TODO: Possibly add a log entry for this
      miss <- complete.cases(env)
      model[['biodiversity']][[id]][['observations']] <- model[['biodiversity']][[id]][['observations']][miss,]
      env <- subset(env, miss)
      # Add intercept
      env$intercept <- 1

      # Add offset if specified
      if(!is.Waiver(x$offset)){
        # Extract offset for each observed point
        ofs <- get_ngbvalue(
          coords = x$biodiversity$get_coordinates(id),
          env = raster::as.data.frame(x$offset, xy = TRUE),
          field_space = c('x','y'),
          longlat = raster::isLonLat(x$background)
        )
        ofs <- subset(ofs, miss)
        assertthat::assert_that(nrow(ofs) == nrow( model$biodiversity[[id]]$observations ))
        model[['biodiversity']][[id]][['offset']] <- ofs
      }

      # Security check
      assertthat::assert_that(
        nrow(env) == nrow( model[['biodiversity']][[id]][['observations']] ),
        'observed' %in% names( model[['biodiversity']][[id]][['observations']] ),
        all( apply(env, 1, function(x) all(!is.na(x) )) ),msg = 'Missing values in extracted environmental predictors.'
      )

      # Check whether predictors should be refined and do so
      if(rm_corPred && model[['predictors_names']] != 'dummy'){
        message('Removing highly correlated variables...')
        test <- env;test$x <- NULL;test$y <- NULL;test$intercept <- NULL

        # Ignore variables for which we have priors
        if(!is.Waiver(x$priors)) keep <- as.character(x$priors$varnames()) else keep <- NULL

        co <- find_correlated_predictors(env = test,
                                         keep = keep,
                                         cutoff = 0.9,
                                         method = 'pear')
        if(length(co)>0){
          env %>% dplyr::select(-dplyr::all_of(co)) -> env
        }
      } else { co <- NULL }

      # Save predictors extracted for biodiversity extraction
      model[['biodiversity']][[id]][['predictors']] <- env
      model[['biodiversity']][[id]][['predictors_names']] <- model[['predictors_names']][which( model[['predictors_names']] %notin% co )]
      model[['biodiversity']][[id]][['predictors_types']] <- model[['predictors_types']][model[['predictors_types']]$predictors %notin% co,]
    }

    # Get and assign Priors
    # FIXME: Type-specific priors?
    if(!is.Waiver(x$priors)){
      # First clean and remove all priors that are not relevant to the engine
      spec_priors <- switch(
        x$engine$name,
        "<GDB>" = x$priors$classes() == 'GDBPrior',
        "<INLA>" = x$priors$classes() == 'INLAPrior'
      )
      spec_priors <- x$priors$collect( names(which(spec_priors)) )
    } else { spec_priors <- new_waiver() }
    model[['priors']] <- spec_priors

    # Calculate latent variables
    if(!is.Waiver(x$latentfactors)){
      # Calculate latent spatial factor (saved in engine data)
      if( length( grep('Spatial',x$get_latent() ) ) > 0 ){
        # Calculate the spatial model
        x$engine$calc_latent_spatial(type = attr(x$get_latent(),'spatial_model'), priors = model[['priors']] )
      }
    }

    # Set offset if existing
    # FIXME: Type-specific offset?
    if(!is.Waiver(x$offset)){
      # Check that they align
      if(!is_comparable_raster(x$offset, x$predictors$data) ){
        new <- raster::resample(x$offset,x$predictors$data)
        assertthat::assert_that(compareRaster(new,x$predictors$data))
        suppressWarnings( x <- x$set_offset(new) )
      }
      # Save overall offset
      model[['offset']] <- as.data.frame(x$offset, xy = TRUE)
    }

    # ----------------- #
    # Number of dataset types, families and ids
    types <- as.character( sapply(model$biodiversity, function(x) x$type) )
    fams <- as.character( sapply( model$biodiversity, function(z) z$family ) )
    ids <- names(model$biodiversity)
    # Engine specific preparations
    #### INLA Engine ####
    if( inherits(x$engine,'INLA-Engine') ){

      # Process per supplied dataset
      for(id in ids) {

        # Default equation found
        if(model$biodiversity[[id]]$equation=='<Default>'){
          # Check potential for rw2 fits
          # MJ: Can lead to to convergence issues with rw2. Removed for now
          # var_rw2 <- apply(poipo_env[model[['predictors_names']]], 2, function(x) length(unique(x)))
          # var_rw2 <- names(which(var_rw2 > 100)) # Get only those greater than 100 unique values
          var_lin <- model$biodiversity[[id]][['predictors_names']] #[which( model[['predictors_names']] %notin% var_rw2 )]

          # Construct formula with all variables
          form <- paste('observed', '~', 0, '+ intercept + ',
                        ifelse(length(types)==1, # Check whether a single intercept model is to be constructed
                               '',
                               paste(paste0('intercept_',sapply( model$biodiversity, function(x) x$type ),collapse = ' + '),' + ')
                               )
                        )
          # Check whether priors have been specified and if yes, use those
          if(!is.Waiver(model$priors)){
            # Loop through all provided INLA priors
            supplied_priors <- as.vector(model$priors$varnames())

            for(v in supplied_priors){
              # FIXME: This currently only work with normal, e.g. the type of the prior is ignored
              # pobj <- model$priors$priors[[model$priors$exists(v)]]
              # First add linear effects
              form <- paste(form, paste0('f(', v, ', model = \'linear\', ',
                                         'mean.linear = ', model$priors$get(v)[1],', ',
                                         'prec.linear = ', model$priors$get(v)[2],')',
                                         collapse = ' + ' ), ' + ' )
            }
            # Add linear for those missed ones
            miss <- var_lin[var_lin %notin% supplied_priors]
            if(length(miss)>0){
              # Add linear predictors without priors
              form <- paste(form,
                            paste('f(', miss,', model = \'linear\')', collapse = ' + ' )
              )
            }
          } else {
            # No priors specified, simply add variables with default
            # Linear for those with few observations
            if(length(var_lin)>0){
              form <- paste(form, ' + ',paste('f(', var_lin,', model = \'linear\')', collapse = ' + ' ) )
            }
            # Random walk 2 (spline) where feasible
            # if(length(var_rw2)>0){
            #   f <- paste(f, ' + ', paste('f(INLA::inla.group(', var_rw2,', n = 10, method = \'quantile\'), model = \'rw2\')', collapse = ' + ' ) )
            # }
          }
          form <- to_formula(form) # Convert to formula
          # Add offset if specified
          if(!is.Waiver(x$offset)){ form <- update.formula(form, paste0('~ . + offset(log(',x$get_offset(),'))') ) }
          if( length( grep('Spatial',x$get_latent() ) ) > 0 ){
            # Update with spatial term
            form <- update.formula(form, paste0(" ~ . + ",
                                          x$engine$get_equation_latent_spatial(
                                            spatial_model = attr(x$get_latent(),'spatial_model'))
                                          )
            )
          }
        } else{
          # If custom supplied formula, check that variable names match supplied predictors
          assertthat::assert_that(
            all( all.vars(form) %in% c('observed',
                                       paste0('intercept_',sapply( model$biodiversity, function(x) x$type )),
                                       model$biodiversity[[id]]$predictors_names) )
          )
          # FIXME: check that
          form <- to_formula( model$biodiversity[[id]]$equation )
        }
        # Update model formula in the model container
        model$biodiversity[[id]]$equation <- form
        rm(form)

        # For each type include expected data
        # expectation vector (area for integration points/nodes and 0 for presences)
        if(model$biodiversity[[id]]$family == 'poisson') model$biodiversity[[id]][['expect']] <- rep(0, nrow(model$biodiversity[[id]]$predictors) )
        if(model$biodiversity[[id]]$family == 'binomial') model$biodiversity[[id]][['expect']] <- rep(1, nrow(model$biodiversity[[id]]$predictors) )
      }

      # Run the engine setup script
      x$engine$setup(model)

      # Now train the model and create a predicted distribution model
      out <- x$engine$train(model,varsel = varsel, inference_only, ...)

      # ----------------------------------------------------------- #
      #### GDB Engine ####
    } else if( inherits(x$engine,"GDB-Engine") ){
      # GDB does not support joint likelihood models
      # Therefore the strategy is to either run the model on its own
      # or in sequence and informing the model following Miller et al. 2019

      # Single model type specified. Define formula
      if(length(types)==1){

        # Default equation found
        if(model$biodiversity[[1]]$equation=='<Default>'){
          # Construct formula with all variables
          form <- paste( 'observed' ,ifelse(model$biodiversity[[1]]$family=='poisson', '/w',''), '~ ')
          if(!is.Waiver(model$priors)){
            # Loop through all provided GDB priors
            supplied_priors <- as.vector(model$priors$varnames())
            for(v in supplied_priors){
              # First add linear effects
              form <- paste(form, paste0('bmono(', v,
                                         ', constraint = \'', model$priors$get(v) ,'\'',
                                         ')', collapse = ' + ' ), ' + ' )
            }
            # Add linear and smooth effects for all missing ones
            miss <- model$biodiversity[[1]]$predictors_names[model$biodiversity[[1]]$predictors_names %notin% supplied_priors]
            if(length(miss)>0){
              # Add linear predictors
              form <- paste(form, paste0('bols(',miss,')',collapse = ' + '), ' + ')
              # And smooth effects
              form <- paste(form, paste0('bbs(', miss,', knots = 4)',
                                         collapse = ' + '
              ))
            }
          } else {
            # Add linear predictors
            form <- paste(form, paste0('bols(',model$biodiversity[[1]]$predictors_names,')',collapse = ' + '), ' + ')
            # And smooth effects
            form <- paste(form, paste0('bbs(',
                                       model$biodiversity[[1]]$predictors_types$predictors[which(model$biodiversity[[1]]$predictors_types$type == 'numeric')],', knots = 4)',
                                       collapse = ' + '
            ))
          }
          # Convert to formula
          form <- to_formula(form)
          # Add offset if specified
          if(!is.Waiver(x$offset)){ form <- update.formula(form, paste0('~ . + offset(log(',x$get_offset(),'))') ) }
          if( length( grep('Spatial',x$get_latent() ) ) > 0 ){
            # Update with spatial term
            form <- update.formula(form, paste0(" ~ . + ",
                                                x$engine$get_equation_latent_spatial() )
            )
          }
        } else{
          # FIXME: Also make checks for correctness in supplied formula, e.g. if variable is contained within object
          form <- to_formula(model$biodiversity[[1]]$equation)
          assertthat::assert_that(
            all( all.vars(form) %in% c('observed', model[['predictors_names']]) )
          )
        }
        model$biodiversity[[1]]$equation <- form
        rm(form)

        # Add pseudo-absence points if necessary
        # Include nearest predictor values for each
        if('poipo' %in% types) {

          # Get background layer
          bg <- x$engine$get_data('template')
          assertthat::assert_that(!is.na(cellStats(bg,min)))

          abs <- create_pseudoabsence(
            env = model$predictors,
            presence = model$biodiversity[[1]]$observations,
            template = bg,
            npoints = 1000,
            replace = FALSE
          )
          abs$intercept <- 1 # Add dummy intercept
          # Combine absence and presence and save
          abs_observations <- abs[,c('x','y')]; abs_observations[['observed']] <- 0
          # Furthermore rasterize observed presences
          pres <- raster::rasterize(model$biodiversity[[1]]$predictors[,c('x','y')], bg, fun = 'count', background = 0)
          obs <- cbind( data.frame(observed = raster::extract(pres, model$biodiversity[[1]]$observations[,c('x','y')])),
                        model$biodiversity[[1]]$observations[,c('x','y')] )
          model$biodiversity[[1]]$observations <- rbind(obs, abs_observations)

          # Format out
          df <- rbind(model$biodiversity[[1]]$predictors, abs) %>% subset(., complete.cases(.) )

          # Preprocessing security checks
          assertthat::assert_that( all( names(abs) %in% names(model$biodiversity[[1]]$predictors) ),
                                   all( model$biodiversity[[1]]$observations[['observed']] >= 0 ),
                                   any(!is.na(rbind(obs, abs_observations)[['observed']] )),
                                   nrow(df) == nrow(model$biodiversity[[1]]$observations)
          )
          # Add offset if existent
          if(!is.Waiver(x$offset)) df[[x$get_offset()]] <- raster::extract(x$offset, df[,c('x','y')])

          # Define expectation as very small vector following Renner et al.
          w <- ppm_weights(df = df,
                           pa = model$biodiversity[[1]]$observations[['observed']],
                           bg = bg,
                           weight = 1e-4
                           )
          df$w <- w # Also add as column

          model$biodiversity[[1]]$predictors <- df
          model$biodiversity[[1]]$expect <- w
        }

      } else {
        # Check whether binomial is present. Construct a binomial model and use
        # its predictions to inform the PPM
        if(length(unique(fams)) == 1 || length(fams)>2) stop('Currently only two types of datasets implemented.')

        # Process binomial
        id = which(fams=='binomial')
        # Default equation found
        if(model$biodiversity[[id]]$equation=='<Default>'){
          # Construct formula with all variables
          form <- paste( 'observed' ,'~ ')
          if(!is.Waiver(model$priors)){
            # Loop through all provided GDB priors
            supplied_priors <- as.vector(model$priors$varnames())
            for(v in supplied_priors){
              # First add linear effects
              form <- paste(form, paste0('bmono(', v,
                                         ', constraint = \'', model$priors$get(v) ,'\'',
                                         ')', collapse = ' + ' ), ' + ' )
            }
            # Add linear and smooth effects for all missing ones
            miss <- model$biodiversity[[id]]$predictors_names[model$biodiversity[[id]]$predictors_names %notin% supplied_priors]
            if(length(miss)>0){
              # Add linear predictors
              form <- paste(form, paste0('bols(',miss,')',collapse = ' + '), ' + ')
              # And smooth effects
              form <- paste(form, paste0('bbs(', miss,', knots = 5)',
                                         collapse = ' + '
              ))
            }
          } else {
            # Add linear predictors
            form <- paste(form, paste0('bols(',model$biodiversity[[id]]$predictors_names,')',collapse = ' + '), ' + ')
            # And smooth effects
            form <- paste(form, paste0('bbs(',
                                       model$biodiversity[[id]]$predictors_types$predictors[which(model$biodiversity[[id]]$predictors_types$type == 'numeric')],', knots = 5)',
                                       collapse = ' + '
            ))
          }
          # Convert to formula
          form <- to_formula(form)
          # Add offset if specified
          # TODO: Think of whether the offset should be set here or rather below
          if(!is.Waiver(x$offset)){ form <- update.formula(form, paste0('~ . + offset(log(',x$get_offset(),'))') ) }
          if( length( grep('Spatial',x$get_latent() ) ) > 0 ){
            # Update with spatial term
            form <- update.formula(form, paste0(" ~ . + ",
                                                x$engine$get_equation_latent_spatial() )
            )
          }
        } else{
          # FIXME: Also make checks for correctness in supplied formula, e.g. if variable is contained within object
          form <- to_formula(model$biodiversity[[id]]$equation)
          assertthat::assert_that(
            all( all.vars(form) %in% c('observed', model[['predictors_names']]) )
          )
        }
        model$biodiversity[[id]]$equation <- form
        rm(form)

        # Now setup and model for poipa only
        model2 <- model
        model2$biodiversity[[ which(fams!='binomial') ]] <- NULL
        # Run the engine setup script
        x$engine$setup(model2)

        # Now train the model and create a predicted distribution model using the first biodiversity dataset
        out <- x$engine$train(model2, inference_only = FALSE)
        # Also create a class prediction of the binomial outcome
        new <- out$fits$prediction
        # Binary layer caused trouble (cholesky errors, thus removing)
                     # predict_gdbclass(fit = out$fits$fit_best,
                     #                  nd = model2$predictors,
                     #                  template = out$fits$prediction)
                     # )
        names(new) <- 'poipa_mean'
        rm(model2) # Clean
        # ---- #
        # FIXME: The code below can possibly be summarized in a clever way
        id = which(fams=='poisson')

        # Now extract and use as prediction
        env <- as.data.frame( raster::extract(new, model$biodiversity[[id]]$observations[,c('x','y')]) )
        names(env) <- 'poipa_mean'
        # Join in to existing predictors
        assertthat::assert_that(nrow(env) == nrow(model$biodiversity[[id]]$predictors),
                                nrow(model$predictors) == nrow(as.data.frame(new)))

        model$biodiversity[[id]]$predictors <- cbind(model$biodiversity[[id]]$predictors, env)
        model$biodiversity[[id]]$predictors_names <- c(model$biodiversity[[id]]$predictors_names, names(env))
        model$biodiversity[[id]]$predictors_types <- rbind(model$biodiversity[[id]]$predictors_types, data.frame(predictors = names(new), type = c('numeric')))
        model$predictors <- cbind(model$predictors, data.frame(poipa_mean = new[]) )
        model$predictors_names <- c(model$predictors_names, names(new))
        model$predictors_types <- rbind(model$predictors_types, data.frame(predictors = names(new), type = c('numeric')))

        # Set monotonic priors
        if(is.Waiver(model$priors)){
          model$priors <- priors(GDBPrior(names(new)[1],hyper = 'increasing'))
        } else {
          model$priors <- model$priors$combine(priors(GDBPrior(names(new)[1],hyper = 'increasing')))
        }
        # Now set up the PPM model
        # TODO: Duplication of above code. Functionize further to reduce code!
        if(model$biodiversity[[id]]$equation=='<Default>'){
          # Construct formula with all variables
          form <- paste( 'observed' ,ifelse(model$biodiversity[[id]]$family=='poisson', '/w',''), '~ ')
          if(!is.Waiver(model$priors)){
            # Loop through all provided GDB priors
            supplied_priors <- as.vector(model$priors$varnames())
            for(v in supplied_priors){
              # First add linear effects
              form <- paste(form, paste0('bmono(', v,
                                         ', constraint = \'', model$priors$get(v) ,'\'',
                                         ')', collapse = ' + ' ), ' + ' )
            }
            # Add linear and smooth effects for all missing ones
            miss <- model$biodiversity[[id]]$predictors_names[model$biodiversity[[id]]$predictors_names %notin% supplied_priors]
            if(length(miss)>0){
              # Add linear predictors
              form <- paste(form, paste0('bols(',miss,')',collapse = ' + '), ' + ')
              # And smooth effects
              form <- paste(form, paste0('bbs(', miss,', knots = 5)',
                                         collapse = ' + '
              ))
            }
          }
          # Convert to formula
          form <- to_formula(form)
          # Add offset if specified
          if(!is.Waiver(x$offset)){ form <- update.formula(form, paste0('~ . + offset(log(',x$get_offset(),'))') ) }
          if( length( grep('Spatial',x$get_latent() ) ) > 0 ){
            # Update with spatial term
            form <- update.formula(form, paste0(" ~ . + ",
                                                x$engine$get_equation_latent_spatial() )
            )
          }
        }
        model$biodiversity[[id]]$equation <- form
        rm(form)

        # --- #
        # Add pseudo-absence points
        bg <- x$engine$get_data('template')
        assertthat::assert_that(!is.na(cellStats(bg,min)))

        # FIXME: Ideally don't sample in regions where there presences of the other dataset
        abs <- create_pseudoabsence(
          env = model$predictors,
          presence = model$biodiversity[[id]]$observations,
          template = x$engine$get_data('template'),
          npoints = 1000,
          replace = FALSE
        )
        abs$intercept <- 1 # Add dummy intercept
        # Combine absence and presence and save
        abs_observations <- abs[,c('x','y')]; abs_observations[['observed']] <- 0
        # Furthermore rasterize observed presences
        pres <- raster::rasterize(model$biodiversity[[id]]$predictors[,c('x','y')], bg, fun = 'count', background = 0)
        obs <- cbind( data.frame(observed = raster::extract(pres, model$biodiversity[[id]]$observations[,c('x','y')])),
                      model$biodiversity[[id]]$observations[,c('x','y')] )
        model$biodiversity[[id]]$observations <- rbind(obs, abs_observations)

        assertthat::assert_that( all( names(abs) %in% names(model$biodiversity[[id]]$predictors) ) )
        # Format out
        df <- rbind(model$biodiversity[[id]]$predictors, abs) %>% subset(., complete.cases(.) )

        # Add offset if existent
        if(!is.Waiver(x$offset)) df[[x$get_offset()]] <- raster::extract(x$offset, df[,c('x','y')])

        # Define expectation as very small vector following Renner et al.
        w <- ppm_weights(df = df,
                         pa = model$biodiversity[[id]]$observations[['observed']],
                         bg = bg,
                         weight = 1e-6
        )
        df$w <- w # Also add as column

        model$biodiversity[[id]]$predictors <- df
        model$biodiversity[[id]]$expect <- w

        # Finally remove the first biodiversity dataset from the model
        model$biodiversity[[which(fams!='poisson')]] <- NULL
      } # Multiple likelihoods

      # Run the engine setup script
      x$engine$setup(model)

      # Now train the model and create a predicted distribution model
      out <- x$engine$train(model, inference_only, ...)

    } else { stop('Specified Engine not implemented yet.') }

    # return output object
    return(out)
  }
)
