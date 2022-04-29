#' @include utils.R bdproto-biodiversitydistribution.R utils-spatial.R
NULL

#' Train the model from a given engine
#'
#' @description
#' This function trains a [distribution()] model with the specified engine and
#' furthermore has some generic parameters that apply to all engines (regardless of type).
#'
#' Users are advised to check the help files for individual [engine]s for advice on how
#' the estimation is being done.
#' @details
#' The resulting object contains both a [`fit_best`] object of the estimated model and, if \code{inference_only} is \code{FALSE}
#' a [RasterLayer] object named [`prediction`] that contains the spatial prediction of the model.
#' These objects can be requested via \code{object$get_data("fit_best")}.
#'
#' Note that several engines do not support fully integrated likelihood models (such [engine_gdb] or
#' [engine_xgboost]) for example. In this case the estimation strategy for multiple added biodiversity datasets
#' (e.g. more than 1) is to fit recurrent models in sequence, where
#' the \code{nth'} spatial projection is provided as input to the \code{nth'+1} dataset.
#' See also Miller et al. (2019) in the references for more details on this strategy. Of course,
#' if users want more control about this aspect, another option is to fit separate models
#' and make use of the [add_offset], [add_offset_range] and [ensemble] functionalities.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object).
#' @param runname A [`character`] name of the trained run.
#' @param rm_corPred Remove highly correlated predictors (Default: \code{FALSE}). This option
#' removes - based on pairwise comparisons - those covariates that are highly collinear (Pearson's r = \code{>0.7}).
#' @param varsel Perform a variable selection on the set of predictors either prior to building the model
#' or via variable selection / regularization of the model. Available options are:
#' * [`none`] for no or default priors and no extensive hyperparameter search.
#' * [`reg`] Model selection either through DIC or regularization / hyperparameter tuning depending on the
#' engine (Default).
#' * [`abess`] A-priori adaptive best subset selection of covariates via the [abess] package (see References).
#' Note that this effectively fits a separate generalized linear model to reduce the number of covariates.
#' Can be helpful for engines that don't directly support efficient variable regularization and when \code{N>100}.
#' @param inference_only By default the [engine] is used to create
#' a spatial prediction of the suitability surface, which can take time. If only inferences of
#' the strength of relationship between covariates and observations are required, this parameter
#' can be set to \code{TRUE} to ignore any spatial projection (Default: \code{FALSE}).
#' @param only_linear Fit model only on linear baselearners and functions. Depending
#' on the [engine] setting this option to \code{FALSE} will result in non-linear relationships
#' between observations and covariates, often increasing processing time (Default: \code{TRUE}).
#' How non-linearity is captured depends on the used [engine].
#' @param bias_variable A [`vector`] with names of variables to be set to *bias_value* (Default: \code{NULL}).
#' This option can for instance be used to 'partial' out certain biases after predictions have been made.
#' See Examples.
#' @param bias_value A [`vector`] with values to be set to *bias_variable* (Default: \code{NULL}).
#' Specifying a [`numeric`] value here sets \code{bias_variable} to the target value.
#' @param ... further arguments passed on.
#' @references
#' * Miller, D.A.W., Pacifici, K., Sanderlin, J.S., Reich, B.J., 2019. The recent past and promising future for data integration methods to estimate species’ distributions. Methods Ecol. Evol. 10, 22–37. https://doi.org/10.1111/2041-210X.13110
#' * Zhu, J., Wen, C., Zhu, J., Zhang, H., & Wang, X. (2020). A polynomial algorithm for best-subset selection problem. Proceedings of the National Academy of Sciences, 117(52), 33117-33123.
#' @seealso [engine_gdb], [engine_xgboost], [engine_bart], [engine_inla], [engine_inlabru], [engine_breg]
#' @returns A [DistributionModel] object.
#' @examples
#' \dontrun{
#'  # Fit a linear penalized logistic regression model via stan
#'  x <- distribution(background) %>%
#'         # Presence-absence data
#'         add_biodiversity_poipa(surveydata) %>%
#'         # Add predictors and scale them
#'         add_predictors(env = predictors, transform = "scale", derivates = "none") %>%
#'         # Use stan for estimation
#'         engine_stan(chains = 2, iter = 1000, warmup = 500)
#'  # Train the model
#'  mod <- train(x, only_linear = TRUE, varsel = 'reg')
#'  mod
#' }
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
  function(x, runname, rm_corPred = FALSE, varsel = "none", inference_only = FALSE,
           only_linear = TRUE,
           bias_variable = NULL, bias_value = NULL, verbose = FALSE,...) standardGeneric("train"))

#' @name train
#' @rdname train
#' @usage \S4method{train}{BiodiversityDistribution}(x)
methods::setMethod(
  "train",
  methods::signature(x = "BiodiversityDistribution", runname = "character"),
  function(x, runname, rm_corPred = FALSE, varsel = "none", inference_only = FALSE,
           only_linear = TRUE,
           bias_variable = NULL, bias_value = NULL, verbose = FALSE,...) {
    # Make load checks
    assertthat::assert_that(
      inherits(x, "BiodiversityDistribution"),
      is.character(runname),
      is.logical(rm_corPred),
      is.logical(inference_only),
      is.null(bias_variable) || is.character(bias_variable),
      is.null(bias_value) || is.numeric(bias_value),
      is.logical(only_linear),
      is.logical(verbose)
    )
    # Now make checks on completeness of the object
    assertthat::assert_that(!is.Waiver(x$engine),
                            msg = 'No engine set for training the distribution model.')
    assertthat::assert_that( x$show_biodiversity_length() > 0,
                             msg = 'No biodiversity data specified.')
    assertthat::assert_that('observed' %notin% x$get_predictor_names(), msg = 'observed is not an allowed predictor name.' )
    if(!is.null(bias_variable)) assertthat::assert_that(bias_variable %in% x$get_predictor_names(),length(bias_variable) == length(bias_value)) else {
      bias_variable <- new_waiver(); bias_value <- new_waiver()
    }
    # Messenger
    if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Collecting input parameters.')
    # --- #
    #rm_corPred = TRUE; varsel = "none"; inference_only = FALSE; verbose = TRUE;only_linear=TRUE;bias_variable = new_waiver();bias_value = new_waiver()
    # Match variable selection
    if(is.logical(varsel)) varsel <- ifelse(varsel, "reg", "none")
    varsel <- match.arg(varsel, c("none", "reg", "abess"), several.ok = FALSE)
    # Define settings object for any other information
    settings <- bdproto(NULL, Settings)
    settings$set('rm_corPred', rm_corPred)
    settings$set('varsel', varsel)
    settings$set('only_linear',only_linear)
    settings$set('inference_only', inference_only)
    settings$set('verbose', verbose)
    settings$set('bias_variable', bias_variable)
    settings$set('bias_value',bias_value)
    # Other settings
    mc <- match.call(expand.dots = FALSE)
    settings$data <- c( settings$data, mc$... )
    # Start time
    settings$set('start.time', Sys.time())

    # Set up logging if specified
    if(!is.Waiver(x$log)) x$log$open()

    # --- #
    #### Defining model objects ----
    # Set model object for fitting
    model <- list()

    # Set model name
    model[['runname']] <- runname

    # Specify a unique id for the run
    model[['id']] <- new_id()
    settings$modelid <- model[['id']]

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
      # Get predictor types
      lu <- sapply(model[['predictors']][model[['predictors_names']]], is.factor)
      model[['predictors_types']] <- data.frame(predictors = names(lu), type = ifelse(lu,'factor', 'numeric') )
      # Assign attribute to predictors to store the name of object
      model[['predictors_object']] <- x$predictors
      rm(lu)
    }

    # Calculate latent variables if set
    if(!is.Waiver(x$latentfactors)){
      # Get the method and check whether it is supported by the engine
      m <- attr(x$get_latent(),'method')
      if(x$get_engine() %notin% c("<INLA>", "<INLABRU>") & m == 'spde'){
        if(getOption('ibis.setupmessages')) myLog('[Setup]','yellow',paste0(m, ' terms are not supported for engine. Switching to poly...'))
        x$set_latent(type = '<Spatial>', 'poly')
      }
      if(x$get_engine()=="<GDB>" & m == 'poly'){
        if(getOption('ibis.setupmessages')) myLog('[Setup]','yellow','Replacing polynominal with P-splines for GDB.')
      }
      if(x$get_engine()=="<BART>" & m == 'car'){
        if(getOption('ibis.setupmessages')) myLog('[Setup]','yellow',paste0(m, ' terms are not supported for engine. Switching to poly...'))
        x$set_latent(type = '<Spatial>', 'poly')
      }
      # Calculate latent spatial terms (saved in engine data)
      if( length( grep('Spatial',x$get_latent() ) ) > 0 ){
        # If model is polynominal, get coordinates of first entry for names of transformation
        if(m == 'poly' & x$get_engine()!="<GDB>"){
          # And the full predictor container
          coords_poly <- polynominal_transform(model$predictors[,c('x','y')], degree = 2)
          model$predictors <- cbind(model$predictors, coords_poly)
          model$predictors_names <- c(model$predictors_names, names(coords_poly))
          model$predictors_types <- rbind(model$predictors_types,
                                          data.frame(predictors = names(coords_poly), type = "numeric"))
          # Also add to predictor object
          pred <- model$predictors_object$get_data(df = FALSE)
          new <-  fill_rasters(coords_poly, emptyraster(pred))
          for(val in names(new)){
            model$predictors_object$set_data(val, new[[val]] )
          }
          rm(pred, new)
        } else {
          # Calculate the spatial model
          x$engine$calc_latent_spatial(type = attr(x$get_latent(),'method'), priors = model[['priors']])
        }
      }
    }

    # Set offset if existing
    if(!is.Waiver(x$offset)){
      # Aggregate offset if necessary
      if(raster::nlayers(x$offset)>1){
        ras_of <- sum(x$offset, na.rm = TRUE)
        names(ras_of) <- "spatial_offset"
      } else {
        ras_of <- x$offset
        names(ras_of) <- "spatial_offset"
      }
      # Save overall offset
      ofs <- as.data.frame(ras_of, xy = TRUE)
      names(ofs)[which(names(ofs)==names(ras_of))] <- "spatial_offset"
      model[['offset']] <- ofs
      # Also add offset object for faster extraction
      model[['offset_object']] <- ras_of
    } else { model[['offset']] <- new_waiver() }

    # Get biodiversity data
    model[['biodiversity']] <- list()
    # Specify list of ids
    biodiversity_ids <- as.character( x$biodiversity$get_ids() )
    for(id in biodiversity_ids) {
      model[['biodiversity']][[id]][['name']]         <- x$biodiversity$data[[id]]$name # Name of the species
      model[['biodiversity']][[id]][['observations']] <- x$biodiversity$get_data(id) # Observational data
      model[['biodiversity']][[id]][['type']]         <- x$biodiversity$get_types(short = TRUE)[[id]] # Type
      model[['biodiversity']][[id]][['family']]       <- x$biodiversity$get_families()[[id]] # Family
      model[['biodiversity']][[id]][['link']]         <- x$biodiversity$get_links()[[id]]
      model[['biodiversity']][[id]][['equation']]     <- x$biodiversity$get_equations()[[id]]
      model[['biodiversity']][[id]][['use_intercept']]<- x$biodiversity$data[[id]]$use_intercept # Separate intercept?
      # --- #
      # Rename observation column to 'observed'. Needs to be consistent for INLA
      # FIXME: try and not use dplyr as dependency (although it is probably loaded already)
      model$biodiversity[[id]]$observations <- model$biodiversity[[id]]$observations %>% dplyr::rename('observed' = x$biodiversity$get_columns_occ()[[id]])
      names(model$biodiversity[[id]]$observations) <- tolower(names(model$biodiversity[[id]]$observations)) # Also generally transfer everything to lower case

      # If the type is polygon, convert to regular sampled points per covered grid cells
      if(any(sf::st_geometry_type(guess_sf(model$biodiversity[[id]]$observations)) %in% c("POLYGON", "MULTIPOLYGON"))){
        o <- polygon_to_points(
          poly = guess_sf(model$biodiversity[[id]]$observations),
          template = emptyraster(x$predictors$get_data(df = FALSE)),
          field_occurrence = "observed" # renamed above
                               )
        model[['biodiversity']][[id]][['observations']] <- o |> as.data.frame()
        model[['biodiversity']][[id]][['type']] <- ifelse(model[['biodiversity']][[id]][['type']] == 'polpo', 'poipo', 'poipa')
        rm(o)
      } else {
        # FIXME: For polygons this won't work. Ideally switch to WKT as default in future
        model$biodiversity[[id]]$observations <- as.data.frame(model$biodiversity[[id]]$observations) # Get only observed column and coordinates
      }

      # Get pseudo-absence information if set, otherwise default options
      if(model[['biodiversity']][[id]][['type']] == "poipo"){
        psa <- x$biodiversity$data[[id]][["pseudoabsence_settings"]]
        if(!is.null(psa)){
          model[['biodiversity']][[id]][['pseudoabsence_settings']] <- psa
        } else { model[['biodiversity']][[id]][['pseudoabsence_settings']] <- getOption("ibis.pseudoabsence")}
      }

      # Now extract coordinates and extract estimates, shifted to raster extraction by default to improve speed!
      env <- get_rastervalue(coords = guess_sf(model$biodiversity[[id]]$observations),
                             env = x$predictors$get_data(df = FALSE),
                             rm.na = FALSE)

      # Remove missing values as several engines can't deal with those easily
      miss <- complete.cases(env)
      model[['biodiversity']][[id]][['observations']] <- model[['biodiversity']][[id]][['observations']][miss,]
      env <- subset(env, miss)
      if(nrow(env)<=2) stop("Too many missing data points in covariates. Check out 'predictor_homogenize_na'.")
      # Add intercept
      env$Intercept <- 1

      # Add offset if specified and model is of poisson type
      if(!is.Waiver(x$offset) ){
        # Extract offset for each observed point
        ofs <- get_rastervalue(
          coords = sf::st_coordinates( guess_sf(model$biodiversity[[id]]$observations) ),
          env = model$offset_object,
          rm.na = FALSE
        )
        ofs <- subset(ofs, miss)
        assertthat::assert_that(nrow(ofs) == nrow( model$biodiversity[[id]]$observations ))
        # Rename
        names(ofs)[which(names(ofs)==names(model$offset_object))] <- "spatial_offset"
        model[['biodiversity']][[id]][['offset']] <- ofs
      }

      # Security check
      assertthat::assert_that(
        nrow(env) == nrow( model[['biodiversity']][[id]][['observations']] ),
        'observed' %in% names( model[['biodiversity']][[id]][['observations']] ),
        all( apply(env, 1, function(x) all(!is.na(x) )) ),msg = 'Missing values in extracted environmental predictors.'
      )

      # Check whether predictors should be refined and do so
      if(settings$get('rm_corPred') && model[['predictors_names']] != 'dummy'){
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','yellow','Removing highly correlated variables...')
        test <- env;test$x <- NULL;test$y <- NULL;test$Intercept <- NULL

        # Ignore variables for which we have priors
        if(!is.Waiver(x$priors)){
          keep <- unique( as.character(x$priors$varnames()) )
          if('spde'%in% keep) keep <- keep[which(keep!='spde')] # Remove SPDE where existing
          test <- test[,-which(names(test) %in% keep)]
          assert_that(!any(keep %in% names(test)))
        } else keep <- NULL

        co <- find_correlated_predictors(env = test,
                                         keep = keep,
                                         cutoff = getOption('ibis.corPred'), # Probably keep default, but maybe sth. to vary in the future
                                         method = 'pearson')

        # For all factor variables, remove those with only the minimal value (e.g. 0)
        fac_min <- apply(test[,model$predictors_types$predictors[which(model$predictors_types$type=='factor')]], 2, function(x) min(x,na.rm = TRUE))
        fac_mean <- apply(test[,model$predictors_types$predictors[which(model$predictors_types$type=='factor')]], 2, function(x) mean(x,na.rm = TRUE))
        co <- unique(co, names(which(fac_mean == fac_min)) ) # Now add to co all those variables where the mean equals the minimum, indicating only absences
        if(length(co)>0){
          env %>% dplyr::select(-dplyr::all_of(co)) -> env
        }
      } else { co <- NULL }

      # Make use of adaptive best subset selection
      if(settings$get("varsel") == "abess"){
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','yellow','Applying abess method to reduce predictors...')
        if(!is.Waiver(x$priors)){
          keep <- unique( as.character(x$priors$varnames()) )
          if('spde'%in% keep) keep <- keep[which(keep!='spde')] # Remove SPDE where existing
        } else keep <- NULL

        # If PPM, calculate points per grid cell first
        if(model[['biodiversity']][[id]]$family == "poisson"){
          bg <- x$engine$get_data("template")
          if(!is.Raster(bg)) bg <- emptyraster(x$predictors$get_data() )

          obs <- aggregate_observations2grid(df = model[['biodiversity']][[id]]$observations,
                                              template = bg,field_occurrence = "observed")

          envs <- get_rastervalue(
            coords = obs[,c('x','y')],
            env = x$predictors$get_data(df = FALSE)[[ model[['predictors_names']][which( model[['predictors_names']] %notin% co )] ]],
            rm.na = TRUE
          )
        } else {
          obs <- model[['biodiversity']][[id]]$observations$observed
          envs <- env[,model[['predictors_names']][which( model[['predictors_names']] %notin% co )]]
        }
        # Add abess here
        co2 <- find_subset_of_predictors(
          env = envs,
          observed = obs$observed,
          family = model[['biodiversity']][[id]]$family,
          tune.type = "cv",
          weight = NULL,
          keep = keep
          )
        co <- c(co, co2) |> unique()
      }

      # Save predictors extracted for biodiversity extraction
      model[['biodiversity']][[id]][['predictors']] <- env
      model[['biodiversity']][[id]][['predictors_names']] <- model[['predictors_names']][which( model[['predictors_names']] %notin% co )]
      model[['biodiversity']][[id]][['predictors_types']] <- model[['predictors_types']][model[['predictors_types']]$predictors %notin% co,]
    }

    # Get and assign Priors
    if(!is.Waiver(x$priors)){
      # First clean and remove all priors that are not relevant to the engine
      spec_priors <- switch(
        x$engine$name,
        "<GDB>" = x$priors$classes() == 'GDBPrior',
        "<XGBOOST>" = x$priors$classes() == 'XGBPrior',
        "<BART>" = x$priors$classes() == 'BARTPrior',
        "<INLA>" = x$priors$classes() == 'INLAPrior',
        "<INLABRU>" = x$priors$classes() == 'INLAPrior',
        "<STAN>" = x$priors$classes() == 'STANPrior',
        "<BREG>" = x$priors$classes() == 'BREGPrior'
      )
      spec_priors <- x$priors$collect( names(which(spec_priors)) )
      # Check whether prior objects match the used engine, otherwise raise warning
      if(spec_priors$length() != x$priors$length()) warning('Some specified priors do not match the engine...')
      # Check whether all priors variables do exist as predictors, otherwise remove
      if(any(spec_priors$varnames() %notin% c( model$predictors_names, 'spde' ))){
        vv <- spec_priors$varnames()[which(spec_priors$varnames() %notin% model$predictors_names)]
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','red',paste0('Some specified priors (',paste(vv, collapse = "|"),') do not match any variable names!') )
        spec_priors$rm( spec_priors$exists(vv) )
      }
    } else { spec_priors <- new_waiver() }
    model[['priors']] <- spec_priors

    # Applying prediction filter based on model input data if specified
    # TODO: Potentially outsource to a function in the future
    if(!is.Waiver(x$limits)){
      # Get biodiversity data
      coords <- do.call(rbind, lapply(model$biodiversity, function(z) z[['observations']][,c('x','y','observed')] ) )
      coords <- subset(coords, observed > 0)
      # Get zones from the limiting area, e.g. those intersecting with input
      suppressMessages(
        suppressWarnings(
          zones <- sf::st_intersection(sf::st_as_sf(coords, coords = c('x','y'),
                                                    crs = sf::st_crs(model$background)),
                                       x$limits)
        )
      )
      # Limit zones
      zones <- subset(x$limits, limit %in% unique(zones$limit) )

      # Now clip all predictors and background to this
      model$background <- suppressMessages(suppressWarnings( sf::st_union( sf::st_intersection(zones, model$background), by_feature = TRUE) )) %>%
        sf::st_cast("MULTIPOLYGON")

      # Extract predictors and offsets again if set
      if(!is.Waiver(model$predictors_object)){
        # Using the raster operations is generally faster than point in polygon tests
        pred_ov <- model$predictors_object$get_data()
        # Make a rasterized mask of the background
        pred_ov <- raster::mask( pred_ov, model$background )
        # Convert Predictors to data.frame
        model[['predictors']] <- raster::as.data.frame(pred_ov, xy = TRUE)
      } else {
        model$predictors[which( is.na(
          point_in_polygon(poly = model$background, points = model$predictors[,c('x','y')] )[['limit']]
        )),model$predictors_names] <- NA # Fill with NA
      }
      # The same with offset if specified
      if(!is.Waiver(x$offset)){
        model$offset[which( is.na(
          point_in_polygon(poly = zones, points = model$offset[,c('x','y')] )[['limit']]
        )), "spatial_offset" ] <- NA # Fill with NA
      }
    }
    # Messenger
    if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Adding engine-specific parameters.')

    # --------------------------------------------------------------------- #
    #### Engine specific code starts below                               ####
    # --------------------------------------------------------------------- #
    # Number of dataset types, families and ids
    types <- as.character( sapply( model$biodiversity, function(x) x$type ) )
    fams <- as.character( sapply( model$biodiversity, function(z) z$family ) )
    ids <- names(model$biodiversity)
    # Engine specific preparations
    #### INLA Engine ####
    if( inherits(x$engine,'INLA-Engine') ){

      # Process per supplied dataset
      for(id in ids) {

        # Default equation found
        if(model$biodiversity[[id]]$equation=='<Default>'){
          # Check potential for rw1 fits
          if(settings$get('only_linear') == FALSE){
            # Get Numeric variables
            vf <- model$biodiversity[[id]]$predictors_types$predictors[model$biodiversity[[id]]$predictors_types$type=="numeric"]
            var_rw1 <- vf
            # Set remaining variables to linear, especially if they are factors
            var_lin <- c()
            if(any(model$biodiversity[[id]]$predictors_types$type=="factor")){
              vf <- model$biodiversity[[id]]$predictors_types$predictors[model$biodiversity[[id]]$predictors_types$type=="factor"]
              var_lin <- c(var_lin, names(explode_factor(model$biodiversity[[id]]$predictors[[vf]], vf)) )
            } else {
              var_lin <- model$biodiversity[[id]][['predictors_names']][which( model$biodiversity[[id]][['predictors_names']] %notin% var_rw1 )]
            }
          } else {
            var_rw1 <- c()
            # Set remaining variables to linear
            var_lin <- model$biodiversity[[id]]$predictors_types$predictors[model$biodiversity[[id]]$predictors_types$type=="numeric"]
            # If any factors are present, split them and add too
            if(any(model$biodiversity[[id]]$predictors_types$type=="factor")){
              vf <- model$biodiversity[[id]]$predictors_types$predictors[model$biodiversity[[id]]$predictors_types$type=="factor"]
              var_lin <- c(var_lin, names(explode_factor(model$biodiversity[[id]]$predictors[[vf]], vf)) )
            }
          }

          # Construct formula with all variables
          form <- paste('observed', '~ Intercept',
                        ifelse(length(types)>1 && model$biodiversity[[id]]$use_intercept, # Check whether a single intercept model is to be constructed
                               paste(' + ',paste0('Intercept_',
                                                  make.names(tolower(model$biodiversity[[id]]$name)),'_',model$biodiversity[[id]]$type
                                                  # make.names(tolower(sapply( model$biodiversity, function(x) x$name ))),'_', # Make intercept from name
                                                  # sapply( model$biodiversity, function(x) x$type ),
                                                  )#collapse = ' + ')
                                     ),
                               ""
                        )
          )
          # Check whether priors have been specified and if yes, use those
          if(!is.Waiver(model$priors)){
            # Loop through all provided INLA priors
            supplied_priors <- model$priors$ids()

            for(v in supplied_priors){
              # Prior variable name
              vn <- as.character( model$priors$varnames()[v] )
              if(vn == 'spde') next()
              # Prior variable type
              vt <- as.character( model$priors$types()[v] )
              # FIXME: This currently only work with normal, e.g. the type of the prior is ignored
              if(vt == 'clinear'){
                # Constrained linear effect
                form <- paste(form, '+', paste0('f(', vn, ', model = \'clinear\', ',
                                                'range = c(',model$priors$get(vn)[1],',',model$priors$get(vn)[2],') )',
                                                collapse = ' + ' ) )
              } else if(vt == 'normal') {
                # Add linear effects
                form <- paste(form, '+', paste0('f(', vn, ', model = \'linear\', ',
                                                'mean.linear = ', model$priors$get(vn)[1],', ',
                                                'prec.linear = ', model$priors$get(vn)[2],')',
                                                collapse = ' + ' ) )
              } else if(vt == 'pc.prec' || vt == 'loggamma'){
                # Add RW effects with pc priors. PC priors is on the KL distance (difference between probability distributions), P(sigma >2)=0.05
                # Default is a loggamma prior with mu 1, 5e-05. Better would be 1, 0.5 following Caroll 2015
                form <- paste0(form, '+', paste0('f(INLA::inla.group(', vn, '), model = \'rw1\', ',
                                                 # 'scale.model = TRUE,',
                                                 'hyper = list(theta = list(prior = ',vt,', param = c(',model$priors$get(vn)[1],',',model$priors$get(vn)[2],')) )
                                                 )',collapse = ' + ')
                )
              }
            }
            # Add linear for those missed ones
            miss <- c(var_lin, var_rw1)[c(var_lin, var_rw1) %notin% model$priors$varnames()]
            if(length(miss)>0){
              if(any(miss %in% var_lin)){
                # Add linear predictors without priors
                form <- paste(form, ' + ',
                              paste('f(', miss[which(miss%in%var_lin)],', model = \'linear\')', collapse = ' + ')
                )
              }
              if(length(var_rw1)>0 & (any(miss %in% var_rw1))){
                # Random walk where feasible and not already included
                form <- paste(form, ifelse(length(var_lin) == 0,'+',''), paste('f(INLA::inla.group(', miss[which(miss%in%var_rw1)],'),',
                                                                                 # 'scale.model = TRUE, ',
                                                                                 # Add RW effects with pc priors. PC priors is on the KL distance (difference between probability distributions), P(sigma >2)=0.05
                                                                                 # Default is a loggamma prior with mu 1, 5e-05. Better would be 1, 0.5 following Caroll 2015
                                                                                 'hyper = list(theta = list(prior = \'loggamma\', param = c(1, 0.5))),',
                                                                                 'model = \'rw1\')', collapse = ' + ' ) )
              }
            }
          } else {
            # No priors specified, simply add variables with default
            # Linear for those with few observations
            if(length(var_lin)>0){
              form <- paste(form, ' + ',paste('f(', var_lin,', model = \'linear\')', collapse = ' + ' ) )
            }
            # Random walk where feasible
            if(length(var_rw1)>0){
              form <- paste(form, ifelse(length(var_lin) == 0,'+',''), paste('f(INLA::inla.group(', var_rw1,'),',
                                                                               # 'scale.model = TRUE,',
                                                                               # Add RW effects with pc priors. PC priors is on the KL distance (difference between probability distributions), P(sigma >2)=0.05
                                                                               # Default is a loggamma prior with mu 1, 5e-05. Better would be 1, 0.5 following Caroll 2015
                                                                               'hyper = list(theta = list(prior = \'loggamma\', param = c(1, 0.5))),',
                                                                               'model = \'rw1\')', collapse = ' + ' ) )
            }
          }
          form <- to_formula(form) # Convert to formula
          # Add offset if specified
          if(!is.Waiver(x$offset) ){ form <- update.formula(form, paste0('~ . + offset(spatial_offset)') ) }
          if( length( grep('Spatial',x$get_latent() ) ) > 0 ){
            if(attr(x$get_latent(), "method") != "poly"){
              # Update with spatial term
              form <- update.formula(form, paste0(" ~ . + ",
                                                  x$engine$get_equation_latent_spatial(
                                                    method = attr(x$get_latent(),'method'),
                                                    vars = which(ids == id),
                                                    separate_spde = attr(x$get_latent(),'separate_spde')
                                                  )
                )
              )
            }
          }
        } else{
          # If custom supplied formula, check that variable names match supplied predictors
          assertthat::assert_that(
            all( all.vars(form) %in% c('observed','Intercept',
                                       paste0('Intercept_',sapply( model$biodiversity, function(x) x$type )),
                                       model$biodiversity[[id]]$predictors_names) )
          )
          # FIXME: check that
          # TODO: Remove elements from predictors that are not used in the formula
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
      # FIXME: Do some checks on whether an observation falls into the mesh?
      model <- x$engine$setup(model, settings)

      # Now train the model and create a predicted distribution model
      out <- x$engine$train(model, settings)

      # ----------------------------------------------------------- #
      #### INLABRU Engine ####
    } else if( inherits(x$engine,'INLABRU-Engine') ){

      # Process per supplied dataset
      for(id in ids) {

        # Default equation found (e.g. no separate specification of effects)
        if(model$biodiversity[[id]]$equation=='<Default>'){
          # Check whether to use dataset specific intercepts
          if(length(types)>1 && model$biodiversity[[id]]$use_intercept){
            ii <- paste0('+ Intercept_',
                         make.names(tolower(model$biodiversity[[id]]$name)),'_',model$biodiversity[[id]]$type
                         # make.names(tolower(sapply( model$biodiversity, function(x) x$name ))),'_', # Make intercept from name
                         # sapply( model$biodiversity, function(x) x$type ),
            )
          } else ii <- ""
          # Go through each variable and build formula for likelihood
          form <- to_formula(paste("observed ~ ", "Intercept +",
                               paste(model$biodiversity[[id]]$predictors_names,collapse = " + "),
                               # Check whether a single dataset is provided, otherwise add other intercepts
                               ii,
                               # # If multiple datasets, remove intercept
                               ifelse(length(ids)>1,"-1", ""),
                               collapse = " ")
                          )

          # Add offset if specified
          # TODO: Not quite sure if this formulation works for inlabru predictor expressions
          if(!is.Waiver(x$offset) ){ form <- update.formula(form, paste0('~ . + offset(spatial_offset)') ) }
          if( length( grep('Spatial',x$get_latent() ) ) > 0 ){
            if(attr(x$get_latent(), "method") != "poly"){
              # Update with spatial term
              form <- update.formula(form, paste0(" ~ . + ",
                                                  # For SPDE components, simply add spatial.field
                                                  paste0("spatial.field",which(ids == id))
                                         )
              )
            }
          }
        } else {
          # If custom likelihood formula is provided, check that variable names match supplied predictors
          form <- model$biodiversity[[id]]$equation
          assertthat::assert_that(
            all( all.vars(form) %in% c('observed',
                                       model$biodiversity[[id]]$predictors_names) )
          )
          # Remove non-covered predictors from the predictor names objects
          model$biodiversity[[id]]$predictors_names <- model$biodiversity[[id]]$predictors_names[which(model$biodiversity[[id]]$predictors_names %in% all.vars(form))]
          model$biodiversity[[id]]$predictors_types <- model$biodiversity[[id]]$predictors_types[
            which( model$biodiversity[[id]]$predictors_types$predictors %in% model$biodiversity[[id]]$predictors_names )
          ,]

          # Convert to formula to be safe
          form <- to_formula( model$biodiversity[[id]]$equation )
          # Add generic Intercept if not set in formula
          if("Intercept" %notin% all.vars(form)) form <- update.formula(form, ". ~ . + Intercept")
          # If length of ids is larger than 1, add dataset specific intercept too
          # Check whether to use dataset specific intercepts
          if(length(types)>1 && model$biodiversity[[id]]$use_intercept){
            form <- update.formula(form,
                                   paste0(". ~ . + ",
                                          paste0('Intercept_',
                                                 make.names(tolower(sapply( model$biodiversity, function(x) x$name ))),'_', # Make intercept from name
                                                 sapply( model$biodiversity, function(x) x$type ),collapse = ' + '))
            )
          }

          if( length( grep('Spatial',x$get_latent() ) ) > 0 ){
            if(attr(x$get_latent(), "method") != "poly"){
              # Update with spatial term
              form <- update.formula(form, paste0(" ~ . + ",
                                                  # For SPDE components, simply add spatial.field
                                                  paste0("spatial.field",which(ids == id))
                                                 )
              )
            }
          }

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
      # FIXME: Do some checks on whether an observation falls into the mesh?
      x$engine$setup(model, settings)

      # Now train the model and create a predicted distribution model
      out <- x$engine$train(model, settings)

      # ----------------------------------------------------------- #
      #### GDB Engine ####
    } else if( inherits(x$engine,"GDB-Engine") ){

      # For each formula, process in sequence
      for(id in ids){
        # Default equation found
        if(model$biodiversity[[id]]$equation=='<Default>'){
          # Construct formula with all variables
          form <- "observed ~ "

          # Use only variables that have sufficient covariate range for training
          # Finally check that a minimum of unique numbers are present in the predictor range and if not, remove them
          covariates <- rm_insufficient_covs(model = model$biodiversity[[id]], tr = 5)

          if(!is.Waiver(model$priors)){
            # Loop through all provided GDB priors
            supplied_priors <- as.vector(model$priors$varnames())
            for(v in supplied_priors){
              if(v %notin% covariates) next() # In case the variable has been removed
              # First add linear effects
              form <- paste(form, paste0('bmono(', v,
                                         ', constraint = \'', model$priors$get(v) ,'\'',
                                         ')', collapse = ' + ' ), ' + ' )
            }
            # Add linear and smooth effects for all missing ones
            miss <- covariates[covariates %notin% supplied_priors]
            if(length(miss)>0){
              # Add linear predictors
              form <- paste(form, paste0('bols(',miss,')',collapse = ' + '))
              if(is.Waiver(settings$get('only_linear'))){
                # And smooth effects for all numeric data
                miss <- miss[ miss %in% model$predictors_types$predictors[which(model$predictors_types$type=="numeric")] ]
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
                                                covariates[which(covariates %in% model$biodiversity[[id]]$predictors_types$predictors[model$biodiversity[[id]]$predictors_types$type=="numeric"] )],', knots = 4)',
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
          if(!is.Waiver(x$offset) ){ form <- update.formula(form, paste0('~ . + offset(spatial_offset)') ) }
          if( length( grep('Spatial',x$get_latent() ) ) > 0 ){
            # Update with spatial term
            form <- update.formula(form, paste0(" ~ . + ",
                                                x$engine$get_equation_latent_spatial())
            )
          }
        } else{
          # FIXME: Also make checks for correctness in supplied formula, e.g. if variable is contained within object
          if(getOption('ibis.setupmessages')) myLog('[Estimation]','yellow','Use custom model equation')
          form <- to_formula(model$biodiversity[[id]]$equation)
          # Update formula to weights if forgotten
          if(model$biodiversity[[id]]$family=='poisson') form <- update.formula(form, 'observed ~ .')
          assertthat::assert_that(
            all( all.vars(form) %in% c('observed', model[['predictors_names']]) )
          )
        }
        model$biodiversity[[id]]$equation <- form
        rm(form)

        # Remove those not part of the modelling
        model2 <- model
        model2$biodiversity <- NULL; model2$biodiversity[[id]] <- model$biodiversity[[id]]

        # Run the engine setup script
        model2 <- x$engine$setup(model2, settings)

        # Now train the model and create a predicted distribution model
        settings2 <- settings
        if(id != ids[length(ids)]) settings2$set('inference_only', FALSE)
        out <- x$engine$train(model2, settings2)

        # Add Prediction of model to next object if multiple are supplied
        if(length(ids)>1 && id != ids[length(ids)]){
          # Add to predictors frame
          new <- out$get_data("prediction")
          pred_name <- paste0(model$biodiversity[[id]]$type, "_", make.names(model$biodiversity[[id]]$name),"_mean")
          names(new) <- pred_name

          # Add the object to the overall prediction object
          model$predictors_object$data <- raster::addLayer(model$predictors_object$get_data(), new)

          # Now for each biodiversity dataset and the overall predictors
          # extract and add as variable
          for(k in names(model$biodiversity)){
            env <- as.data.frame(
              raster::extract(new, model$biodiversity[[k]]$observations[,c('x','y')]) )
            # Rename to current id dataset
            names(env) <- pred_name
            # Add
            model$biodiversity[[k]]$predictors <- cbind(model$biodiversity[[k]]$predictors, env)
            model$biodiversity[[k]]$predictors_names <- c(model$biodiversity[[k]]$predictors_names,
                                                          names(env) )
            model$biodiversity[[k]]$predictors_types <- rbind(
              model$biodiversity[[k]]$predictors_types,
              data.frame(predictors = names(env), type = c('numeric'))
            )
          }
          # Add to overall predictors
          model$predictors <- cbind(model$predictors, as.data.frame(new))
          model$predictors_names <- c(model$predictors_names, names(new))
          model$predictors_types <- rbind(model$predictors_types,
                                          data.frame(predictors = names(new), type = c('numeric')))

          # Set monotonic priors
          if(is.Waiver(model$priors)){
            model$priors <- priors(GDBPrior(names(new)[1],hyper = 'increasing'))
          } else {
            model$priors <- model$priors$combine(priors(GDBPrior(names(new)[1],hyper = 'increasing')))
          }

        }

      }
      # ----------------------------------------------------------- #
      #### XGBoost Engine ####
    } else if( inherits(x$engine,"XGBOOST-Engine") ){
      # Create XGBboost regression and classification

      # Process per supplied dataset and in order of supplied data
      for(id in ids) {

        # Default equation found
        if(model$biodiversity[[id]]$equation=='<Default>'){
          # XGboost does not explicitly work with formulas, thus all supplied objects are assumed to be part
          # a covariate
          form <- new_waiver()
          # Note: Priors are added in the fitted distribution object through the model object
        } else{
          # If custom supplied formula, check that variable names match supplied predictors
          stop("Custom formulas not yet implemented")
          # TODO: Remove elements from predictors that are not used in the formula
          form <- new_waiver()
        }
        # Update model formula in the model container
        model$biodiversity[[id]]$equation <- form
        rm(form)

        # Remove those not part of the modelling
        model2 <- model
        model2$biodiversity <- NULL; model2$biodiversity[[id]] <- model$biodiversity[[id]]

        # Run the engine setup script
        model2 <- x$engine$setup(model2, settings)

        # Now train the model and create a predicted distribution model
        settings2 <- settings
        if(id != ids[length(ids)]) settings2$set('inference_only', FALSE)
        out <- x$engine$train(model2, settings2)

        # Add Prediction of model to next object if multiple are supplied
        if(length(ids)>1 && id != ids[length(ids)]){
          # Add to predictors frame
          new <- out$get_data("prediction")
          pred_name <- paste0(model$biodiversity[[id]]$type, "_", make.names(model$biodiversity[[id]]$name),"_mean")
          names(new) <- pred_name
          # Add the object to the overall prediction object
          model$predictors_object$data <- raster::addLayer(model$predictors_object$get_data(), new)
          # Now for each biodiversity dataset and the overall predictors
          # extract and add as variable
          for(k in names(model$biodiversity)){
            env <- as.data.frame(
              raster::extract(new, model$biodiversity[[k]]$observations[,c('x','y')]) )
            # Rename to current id dataset
            names(env) <- pred_name
            # Add
            model$biodiversity[[k]]$predictors <- cbind(model$biodiversity[[k]]$predictors, env)
            model$biodiversity[[k]]$predictors_names <- c(model$biodiversity[[k]]$predictors_names,
                                                          names(env) )
            model$biodiversity[[k]]$predictors_types <- rbind(
                model$biodiversity[[k]]$predictors_types,
                data.frame(predictors = names(env), type = c('numeric'))
              )
          }
          # Add to overall predictors
          model$predictors <- cbind(model$predictors, as.data.frame(new) )
          model$predictors_names <- c(model$predictors_names, names(new))
          model$predictors_types <- rbind(model$predictors_types,
                                          data.frame(predictors = names(new), type = c('numeric')))
        }
      }

      # ----------------------------------------------------------- #
      #### BART Engine ####
    } else if( inherits(x$engine,"BART-Engine") ){

      # Process each id
      for(id in ids){
        # Default equation found
        if(model$biodiversity[[id]]$equation=='<Default>'){
          # Construct formula with all variables
          form <- paste( 'observed' ,ifelse(model$biodiversity[[id]]$family=='poisson', '/w',''), '~ ',
                         paste(model$biodiversity[[id]]$predictors_names,collapse = " + "))
          # Convert to formula
          form <- to_formula(form)
        } else {
          # FIXME: Also make checks for correctness in supplied formula, e.g. if variable is contained within object
          form <- to_formula(model$biodiversity[[id]]$equation)
          assertthat::assert_that(
            all( all.vars(form) %in% c('observed','w', model[['predictors_names']]) )
          )
        }
        model$biodiversity[[id]]$equation <- form
        rm(form)

        # Model and estimate
        model2 <- model
        model2$biodiversity <- NULL; model2$biodiversity[[id]] <- model$biodiversity[[id]]

        # Run the engine setup script
        model2 <- x$engine$setup(model2, settings)

        # Now train the model and create a predicted distribution model
        settings2 <- settings
        if(id != ids[length(ids)]) settings2$set('inference_only', FALSE)
        out <- x$engine$train(model2, settings2)

        # Add Prediction of model to next object if multiple are supplied
        if(length(ids)>1 && id != ids[length(ids)]){
          # Add to predictors frame
          new <- out$get_data("prediction")[[1]] # Only take the mean!
          pred_name <- paste0(model$biodiversity[[id]]$type, "_", make.names(model$biodiversity[[id]]$name),"_mean")
          names(new) <- pred_name
          # Add the object to the overall prediction object
          model$predictors_object$data <- raster::addLayer(model$predictors_object$get_data(), new)

          # Now for each biodiversity dataset and the overall predictors
          # extract and add as variable
          for(k in names(model$biodiversity)){
            env <- as.data.frame(
              raster::extract(new, model$biodiversity[[k]]$observations[,c('x','y')]) )
            # Rename to current id dataset
            names(env) <- pred_name
            # Add
            model$biodiversity[[k]]$predictors <- cbind(model$biodiversity[[k]]$predictors, env)
            model$biodiversity[[k]]$predictors_names <- c(model$biodiversity[[k]]$predictors_names,
                                                          names(env) )
            model$biodiversity[[k]]$predictors_types <- rbind(
              model$biodiversity[[k]]$predictors_types,
              data.frame(predictors = names(env), type = c('numeric'))
            )
          }
          # Add to overall predictors
          model$predictors <- cbind(model$predictors, as.data.frame(new) )
          model$predictors_names <- c(model$predictors_names, names(new))
          model$predictors_types <- rbind(model$predictors_types,
                                          data.frame(predictors = names(new), type = c('numeric')))
        }
      } # End of id loop
    } else if( inherits(x$engine,"STAN-Engine") ){
      # ----------------------------------------------------------- #
      #### STAN Engine ####
      # For stan, the actual model is built sequentially per id
      # Process per supplied dataset
      for(id in ids) {
        # TODO
        if(length(model$biodiversity)>1) stop("Not yet implemented")

        # Default equation found (e.g. no separate specification of effects)
        if(model$biodiversity[[id]]$equation=='<Default>'){

          # Go through each variable and build formula for likelihood
          form <- to_formula(paste("observed",
                                          " ~ ", "0 + ",
                                   ifelse(model$biodiversity[[id]]$family=='poisson', " offset(log(w)) + ", ""), # Use log area as offset
                                   paste(model$biodiversity[[id]]$predictors_names,collapse = " + "),
                                   # Check whether a single dataset is provided, otherwise add other intercepts
                                   ifelse(length(types)==1,
                                          '',
                                          paste('+',paste0('Intercept_',
                                                           make.names(tolower(sapply( model$biodiversity, function(x) x$name ))),'_', # Make intercept from name
                                                           sapply( model$biodiversity, function(x) x$type ),collapse = ' + ')
                                          )
                                   ),
                                   # # If multiple datasets, don't use intercept
                                   # ifelse(length(ids)>1,"-1", ""),
                                   collapse = " ")
                            )

          # Add offset if specified
          if(!is.Waiver(x$offset)){ form <- update.formula(form, paste0('~ . + offset(spatial_offset)') ) }
          # if( length( grep('Spatial',x$get_latent() ) ) > 0 ) {} # Possible to be implemented for CAR models
        } else {
          if(getOption('ibis.setupmessages')) myLog('[Estimation]','yellow','Use custom model equation.')
          form <- to_formula(model$biodiversity[[1]]$equation)
          # Update formula to weights if forgotten
          if(model$biodiversity[[1]]$family=='poisson') form <- update.formula(form, 'observed / w ~ .')
          assertthat::assert_that(
            all( all.vars(form) %in% c('observed','w', model[['predictors_names']]) )
          )
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
      model <- x$engine$setup(model, settings)

      # Now train the model and create a predicted distribution model
      out <- x$engine$train(model, settings)


    } else if (inherits(x$engine,"BREG-Engine") ){
    # ----------------------------------------------------------- #
    #### BREG Engine ####
    # For each formula, process in sequence
    for(id in ids){
      # Default equation found
      if(model$biodiversity[[id]]$equation=='<Default>'){
        # Construct formula with all variables
        form <- paste( 'observed' , ' ~ ')
        # Add linear predictors
        form <- paste(form, paste0(model$biodiversity[[id]]$predictors_names,collapse = ' + '))
        # NOTE: Non-linearity will be specified during engine setup!
        # Convert to formula
        form <- to_formula(form)
      } else{
        # FIXME: Also make checks for correctness in supplied formula, e.g. if variable is contained within object
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','yellow','Use custom model equation')
        form <- to_formula(model$biodiversity[[id]]$equation)
        assertthat::assert_that(
          all( all.vars(form) %in% c('observed', model[['predictors_names']]) )
        )
      }
      model$biodiversity[[id]]$equation <- form
      rm(form)

      # Remove those not part of the modelling
      model2 <- model
      model2$biodiversity <- NULL; model2$biodiversity[[id]] <- model$biodiversity[[id]]

      # Run the engine setup script
      model2 <- x$engine$setup(model2, settings)

      # Now train the model and create a predicted distribution model
      settings2 <- settings
      if(id != ids[length(ids)]) settings2$set('inference_only', FALSE)
      out <- x$engine$train(model2, settings2)

      # Add Prediction of model to next object if multiple are supplied
      if(length(ids)>1 && id != ids[length(ids)]){
        # Add to predictors frame
        new <- out$get_data("prediction")
        pred_name <- paste0(model$biodiversity[[id]]$type, "_", make.names(model$biodiversity[[id]]$name),"_mean")
        names(new) <- pred_name
        # Add the object to the overall prediction object
        model$predictors_object$data <- raster::addLayer(model$predictors_object$get_data(), new)

        # Now for each biodiversity dataset and the overall predictors
        # extract and add as variable
        for(k in names(model$biodiversity)){
          env <- as.data.frame(
            raster::extract(new, model$biodiversity[[k]]$observations[,c('x','y')]) )
          # Rename to current id dataset
          names(env) <- pred_name
          # Add
          model$biodiversity[[k]]$predictors <- cbind(model$biodiversity[[k]]$predictors, env)
          model$biodiversity[[k]]$predictors_names <- c(model$biodiversity[[k]]$predictors_names,
                                                        names(env) )
          model$biodiversity[[k]]$predictors_types <- rbind(
            model$biodiversity[[k]]$predictors_types,
            data.frame(predictors = names(env), type = c('numeric'))
          )
        }
        # Add to overall predictors
        model$predictors <- cbind(model$predictors, as.data.frame(new))
        model$predictors_names <- c(model$predictors_names, names(new))
        model$predictors_types <- rbind(model$predictors_types,
                                        data.frame(predictors = names(new), type = c('numeric')))
       }

    }
  # End of BREG engine
  } else { stop('Specified Engine not implemented yet.') }

  if(getOption('ibis.setupmessages')) myLog('[Done]','green',paste0('Completed after ', round( as.numeric(out$settings$duration()), 2),' ',attr(out$settings$duration(),'units') ))

  # Stop logging if specified
  if(!is.Waiver(x$log)) x$log$close()

  # return output object
  return(out)
  }
)
