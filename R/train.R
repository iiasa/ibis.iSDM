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
  function(x, runname,rm_corPred,varsel,...) standardGeneric("train"))

#' @name train
#' @rdname train
#' @usage \S4method{train}{BiodiversityDistribution}(x)
methods::setMethod(
  "train",
  methods::signature(x = "BiodiversityDistribution", runname = "character"),
  function(x, runname, rm_corPred = TRUE, varsel = TRUE, ...) {
    # Make load checks
    assertthat::assert_that(
      inherits(x, "BiodiversityDistribution"),
      is.character(runname),
      is.logical(rm_corPred)
    )
    # Now make checks on completeness of the object
    assertthat::assert_that(!is.Waiver(x$engine),
                            msg = 'No engine set for training the distribution model.')
    assertthat::assert_that( x$show_biodiversity_length() > 0,
                             msg = 'No biodiversity data specified.')

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

    # Get biodiversity data
    model[['data']] <- list()
    types <- names( x$biodiversity$get_types() )
    for(ty in types) model[['data']][[ty]] <- x$biodiversity$get_data(ty)

    # Observed columns
    for(ty in types) model[['data']][[paste0(ty,'_','response')]] <-  x$biodiversity$get_columns_occ()[[ty]]

    # Get predictors
    if(is.Waiver(x$get_predictor_names())) {
      # Dummy covariate of background raster
      dummy <- as.data.frame( raster(extent(x$background),nrow=100,ncol=100,val=1), xy = TRUE );names(dummy)[3] <- 'dummy'
      model[['predictors']] <- dummy
      model[['predictors_names']] <- 'dummy'
    } else {
      # Convert Predictors to data.frame
      model[['predictors']] <- x$predictors$get_data(df = TRUE, na.rm = FALSE)
      # Also set predictor names
      model[['predictors_names']] <- x$get_predictor_names()
    }

    # Extract estimates for point records
    poipo_env <- get_ngbvalue(
      coords = x$biodiversity$get_coordinates('poipo'),
      env = model[['predictors']],
      field_space = c('x','y'),
      longlat = raster::isLonLat(x$background)
    )

    # Check whether predictors should be refined and do so
    if(rm_corPred && model[['predictors_names']] != 'dummy'){
      message('Removing highly correlated variables...')
      test <- subset(poipo_env, complete.cases(poipo_env));test$x <- NULL;test$y <- NULL
      co <- find_correlated_predictors(env = test,
                                      keep = NULL, cutoff = 0.9, method = 'pear')
      if(length(co)>0){
        poipo_env %>% dplyr::select(-all_of(co)) -> poipo_env
        model[['predictors_names']] <- model[['predictors_names']][which( model[['predictors_names']] %notin% co )]
        model[['predictors']] <- model[['predictors']][which( model[['predictors']] %notin% co )]
      }
    }

    # assign default priors
    if(is.Waiver( x$priors )){
      # TODO: Define prior objects. Also look at PC priors https://www.tandfonline.com/doi/abs/10.1080/01621459.2017.1415907?journalCode=uasa20
      message('TODO: Define prior objects')
      model[['priors']] <- NULL
      #x$set_prior( add_default_priors() )
    }
    # Get latent variables
    if(!is.Waiver(x$latentfactors)){
      # Calculate latent spatial factor (saved in engine data)
      if(x$get_latent()=="<Spatial>") x$engine$calc_latent_spatial(type = 'pc')
    }

    # Format formulas
    model[['equation']] <- list()
    types <- names( x$biodiversity$get_types() )
    for(ty in types) {
      # Default equation found
      if(x$biodiversity$get_equations()[[ty]]=='<Default>'){
        # Construct formula with all variables
        f <- formula(
                  paste( x$biodiversity$get_columns_occ()[[ty]], '~ ',
                         0, #ifelse(x$show_biodiversity_length()==1,1,0),
                         ' +',
                         paste( model[['predictors_names']], collapse = ' + ' )  )
                )
        if(x$get_latent()=="<Spatial>"){
          # Update with spatial term
          f <- update.formula(f, paste0(" ~ . + ",x$engine$get_equation_latent_spatial() ) )
        }
      } else{
        stop('TBD')
        # FIXME: Also make checks for correct formula, e.g. if variable is contained within object
        }
      model[['equation']][[ty]] <- f
      rm(f)
    }

    # Engine specific preparations
    if( inherits(x$engine,'INLA-Engine') ){
      # Include nearest predictor values for each
      types <- names( x$biodiversity$get_types() )
      if('poipo' %in% types) {
        model[['data']][['poipo_values']] <- poipo_env
        # Calculate the expectaction
        # Calculate e as relative area from the points
        # sfmesh <- mesh_as_sf( self$get_data('mesh'))
        # # Set those not intersecting with the background to 0
        # ind <- suppressMessages(
        #   sf::st_join(sfmesh, x$background, join = st_within)[,3] %>% st_drop_geometry()
        # )
        # sfmesh[which( is.na(ind[,1]) ),'relarea'] <- NA
        # e <- suppressMessages(
        #   point_in_polygon(
        #     poly = sfmesh,
        #     points = model$data$poipo
        #   )$relarea
        # )
        model[['data']][['poipo_expect']] <- rep(0, nrow(poipo_env) )
      }

      # Run the engine setup script
      x$engine$setup(model)

      # Now train the model and create a predicted distribution model
      out <- x$engine$train(model,varsel = varsel,...)

    } else { stop('Engines not implemented yet')}

    # return output object
    return(out)
  }
)
