#' @include utils.R bdproto-biodiversitydistribution.R utils-spatial.R
NULL

#' Train the model from a given engine
#'
#' @description
#' This function trains a [distribution()] model with the specified engine and
#' furthermore has some generic options that apply to all engines (regardless of type).
#' See Details with regards to such options.
#'
#' Users are advised to check the help files for individual engines for advice on how
#' the estimation is being done.
#' @details
#' This function acts as a generic training function that - based on the provided [`BiodiversityDistribution-class`]
#' object creates a new distribution model.
#' The resulting object contains both a \code{"fit_best"} object of the estimated model and, if \code{inference_only} is \code{FALSE}
#' a [SpatRaster] object named \code{"prediction"} that contains the spatial prediction of the model.
#' These objects can be requested via \code{object$get_data("fit_best")}.
#'
#' Other parameters in this function:
#'
#' * \code{"filter_predictors"} The parameter can be set to various options to remove highly correlated variables or those
#' with little additional information gain from the model prior to any estimation. Available options are \code{"none"} (Default) \code{"pearson"} for
#' applying a \code{0.7} correlation cutoff, \code{"abess"} for the regularization framework by Zhu et al. (2020), or \code{"RF"} or \code{"randomforest"}
#' for removing the least important variables according to a randomForest model. **Note**: This function is only applied on
#' predictors for which no prior has been provided (e.g. potentially non-informative ones).
#'
#' * \code{"optim_hyperparam"} This option allows to make use of hyper-parameter search for several models, which can improve
#' prediction accuracy although through the a substantial increase in computational cost.
#'
#' * \code{"method_integration"} Only relevant if more than one [`BiodiversityDataset`] is supplied and when
#' the engine does not support joint integration of likelihoods.
#' See also Miller et al. (2019) in the references for more details on different types of integration. Of course,
#' if users want more control about this aspect, another option is to fit separate models
#' and make use of the [add_offset], [add_offset_range] and [ensemble] functionalities.
#'
#' * \code{"clamp"} Boolean parameter to support a clamping of the projection predictors to the range of values observed
#' during model training.
#'
#' @note
#' There are no silver bullets in (correlative) species distribution modelling and for each model the analyst has to
#' understand the objective, workflow and parameters than can be used to modify the outcomes. Different predictions can
#' be obtained from the same data and parameters and not all necessarily make sense or are useful.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object).
#' @param runname A [`character`] name of the trained run.
#' @param filter_predictors A [`character`] defining if and how highly correlated predictors are to be removed
#' prior to any model estimation.
#' Available options are:
#' * \code{"none"} No prior variable removal is performed (Default).
#' * \code{"pearson"}, \code{"spearman"} or \code{"kendall"} Makes use of pairwise comparisons to identify and
#' remove highly collinear predictors (Pearson's \code{r >= 0.7}).
#' * \code{"abess"} A-priori adaptive best subset selection of covariates via the [`abess`] package (see References).
#' Note that this effectively fits a separate generalized linear model to
#' reduce the number of covariates.
#' * \code{"boruta"} Uses the [`Boruta`] package to identify non-informative features.
#'
#' @param optim_hyperparam Parameter to tune the model by iterating over input parameters or selection
#' of predictors included in each iteration. Can be set to \code{TRUE} if extra precision is
#' needed (Default: \code{FALSE}).
#'
#' @param inference_only By default the engine is used to create
#' a spatial prediction of the suitability surface, which can take time. If only inferences of
#' the strength of relationship between covariates and observations are required, this parameter
#' can be set to \code{TRUE} to ignore any spatial projection (Default: \code{FALSE}).
#' @param only_linear Fit model only on linear baselearners and functions. Depending
#' on the engine setting this option to \code{FALSE} will result in non-linear relationships
#' between observations and covariates, often increasing processing time (Default: \code{TRUE}).
#' How non-linearity is captured depends on the used engine.
#' @param method_integration A [`character`] with the type of integration that should be applied if more
#' than one [`BiodiversityDataset-class`] object is provided in \code{x}. Particular relevant for engines
#' that do not support the integration of more than one dataset. Integration methods are generally sensitive
#' to the order in which they have been added to the  [`BiodiversityDistribution`] object.
#'
#' Available options are:
#' * \code{"predictor"} The predicted output of the first (or previously fitted) models are
#' added to the predictor stack and thus are predictors for subsequent models (Default).
#' * \code{"offset"} The predicted output of the first (or previously fitted) models are
#' added as spatial offsets to subsequent models. Offsets are back-transformed depending
#' on the model family. This option might not be supported for every [`engine`].
#' * \code{"interaction"} Instead of fitting several separate models, the observations from each dataset
#' are combined and incorporated in the prediction as a factor interaction with the "weaker" data source being
#' partialed out during prediction. Here the first dataset added determines the reference level
#' (see Leung et al. 2019 for a description).
#' * \code{"prior"} In this option we only make use of the coefficients from a previous model to define priors to be used in the next model.
#' Might not work with any engine!
#' * \code{"weight"} This option only works for multiple biodiversity datasets with the same type (e.g. \code{"poipo"}).
#' Individual weight multipliers can be determined while setting up the model (**Note: Default is 1**). Datasets are then combined for estimation
#' and weighted respectively, thus giving for example presence-only records less weight than survey records.
#'
#' **Note that this parameter is ignored for engines that support joint likelihood estimation.**
#' @param aggregate_observations [`logical`] on whether observations covering the same grid cell should be aggregated (Default: \code{TRUE}).
#' @param clamp [`logical`] whether predictions should be clamped to the range of predictor values observed during model fitting (Default: \code{FALSE}).
#' @param verbose Setting this [`logical`] value to \code{TRUE} prints out further information during the model fitting (Default: \code{FALSE}).
#' @param ... further arguments passed on.
#' @references
#' * Miller, D.A.W., Pacifici, K., Sanderlin, J.S., Reich, B.J., 2019. The recent past and promising future for data integration methods to estimate species’ distributions. Methods Ecol. Evol. 10, 22–37. https://doi.org/10.1111/2041-210X.13110
#' * Zhu, J., Wen, C., Zhu, J., Zhang, H., & Wang, X. (2020). A polynomial algorithm for best-subset selection problem. Proceedings of the National Academy of Sciences, 117(52), 33117-33123.
#' * Leung, B., Hudgins, E. J., Potapova, A. & Ruiz‐Jaen, M. C. A new baseline for countrywide α‐diversity and species distributions: illustration using &gt;6,000 plant species in Panama. Ecol. Appl. 29, 1–13 (2019).
#' @seealso [engine_gdb], [engine_xgboost], [engine_bart], [engine_inla], [engine_inlabru], [engine_breg], [engine_stan]
#' @returns A [DistributionModel] object.
#' @examples
#' \dontrun{
#'  # Fit a linear penalized logistic regression model via stan
#'  x <- distribution(background) |>
#'         # Presence-absence data
#'         add_biodiversity_poipa(surveydata) |>
#'         # Add predictors and scale them
#'         add_predictors(env = predictors, transform = "scale", derivates = "none") |>
#'         # Use Stan for estimation
#'         engine_stan(chains = 2, iter = 1000, warmup = 500)
#'  # Train the model
#'  mod <- train(x, only_linear = TRUE, filter_predictors = 'pearson')
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
  signature = methods::signature("x"),
  function(x, runname, filter_predictors = "none", optim_hyperparam = FALSE, inference_only = FALSE,
           only_linear = TRUE, method_integration = "predictor",
           aggregate_observations = TRUE, clamp = FALSE, verbose = FALSE,...) standardGeneric("train"))

#' @name train
#' @rdname train
#' @usage \S4method{train}{BiodiversityDistribution, character, character, logical, logical, logical, character, logical, logical, logical}(x,runname,filter_predictors,optim_hyperparam,inference_only,only_linear,method_integration,aggregate_observations,clamp,verbose)
methods::setMethod(
  "train",
  methods::signature(x = "BiodiversityDistribution"),
  function(x, runname, filter_predictors = "none", optim_hyperparam = FALSE, inference_only = FALSE,
           only_linear = TRUE, method_integration = "predictor",
           aggregate_observations = TRUE, clamp = FALSE, verbose = FALSE,...) {
    if(missing(runname)) runname <- "Unnamed run"

    # Make load checks
    assertthat::assert_that(
      inherits(x, "BiodiversityDistribution"),
      is.character(runname),
      is.logical(optim_hyperparam),
      is.character(filter_predictors),
      is.logical(inference_only),
      is.logical(only_linear),
      is.character(method_integration),
      is.logical(clamp),
      is.logical(verbose)
    )
    # Now make checks on completeness of the object
    assertthat::assert_that(!is.Waiver(x$engine),
                            msg = 'No engine set for training the distribution model.')
    assertthat::assert_that( x$show_biodiversity_length() > 0,
                             msg = 'No biodiversity data specified.')
    assertthat::assert_that('observed' %notin% x$get_predictor_names(), msg = 'observed is not an allowed predictor name.' )
    # Messenger
    if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Collecting input parameters.')
    # --- #
    #filter_predictors = "none"; optim_hyperparam = FALSE; runname = "test";inference_only = FALSE; verbose = TRUE;only_linear=TRUE;method_integration="predictor";aggregate_observations = TRUE; clamp = FALSE
    # Match variable selection
    filter_predictors <- match.arg(filter_predictors, c("none", "pearson", "spearman", "kendall", "abess", "RF", "randomForest", "boruta"), several.ok = FALSE)
    method_integration <- match.arg(method_integration, c("predictor", "offset", "interaction", "prior", "weight"), several.ok = FALSE)
    # Define settings object for any other information
    settings <- bdproto(NULL, Settings)
    settings$set('filter_predictors', filter_predictors)
    settings$set('optim_hyperparam', optim_hyperparam)
    settings$set('only_linear',only_linear)
    settings$set('inference_only', inference_only)
    settings$set('clamp', clamp)
    settings$set('ibis.cleannames', getOption("ibis.cleannames"))
    settings$set('verbose', verbose)
    settings$set('seed', getOption("ibis.seed"))
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
      if(getOption('ibis.setupmessages')) myLog('[Setup]','yellow',paste0('No predictor terms found. Using dummy.'))
      # Dummy covariate of background raster
      # Check if the engine has a template and if so use that one
      if(is.Raster(x$engine$get_data("template"))){
        dummy <- emptyraster(x$engine$get_data("template"));names(dummy) <- "dummy"
        dummy[] <- 1 ; dummy <- terra::mask(dummy, x$background)
      } else {
        dummy <- terra::rast( terra::ext(x$background),
                              nrow=100, ncol=100, val=1,
                              crs = terra::crs(x$background));names(dummy) <- 'dummy'
      }
      model[['predictors']] <- terra::as.data.frame(dummy, xy = TRUE, na.rm = FALSE)
      model[['predictors_names']] <- 'dummy'
      model[['predictors_types']] <- data.frame(predictors = 'dummy', type = 'numeric')
      model[['predictors_object']] <- bdproto(NULL, PredictorDataset, id = new_id(), data = dummy)
    } else {
      # Convert Predictors to data.frame
      model[['predictors']] <- x$predictors$get_data(df = TRUE, na.rm = FALSE)

      # Check whether any of the variables are fully NA, if so exclude
      if( any( apply(model[['predictors']], 2, function(z) all(is.na(z))) )){
        chk <- which( apply(model[['predictors']], 2, function(z) all(is.na(z))) )
        if(getOption('ibis.setupmessages')) myLog('[Setup]','red',
                                                  paste0('The following variables are fully missing and are removed:\n',
                                                         paste(names(chk),collapse = " | "))
        )
        model[['predictors']] <- model[['predictors']][,-chk]
        x$predictors$rm_data(names(chk)) # Remove the variables
      }

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
      model[["latent"]] <- attr(x$latentfactors, "method") # Save type for the record
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
            model$predictors_object <- model$predictors_object$set_data(val, new[[val]] )
          }
          rm(pred, new)
        } else if(m == "kde") {
          # Bivariate kernel density estimation
          # First get all points
          biodiversity_ids <- as.character( x$biodiversity$get_ids() )
          poi <- data.frame()
          for(id in biodiversity_ids) {
            # Get presence points
            o <- guess_sf( x$biodiversity$get_data(id) )
            o <- o[ which( o[[x$biodiversity$get_columns_occ()[[id]]]] > 0 ), ]
            o <- subset(o, select = "geometry")
            poi <- rbind(poi, o)
          }
          # Ensure we have a backgroudn raster
          if(inherits(x$background, "sf")){
            bg <- terra::rasterize(x$background, model$predictors_object$get_data(), 1)
          } else {bg <- x$background }
          # Then calculate
          ras <- st_kde(points = poi, background = bg, bandwidth = 3)
          # Add to predictor objects, names, types and the object
          model[['predictors']] <- cbind.data.frame( model[['predictors']], terra::as.data.frame(ras, na.rm = FALSE) )
          model[['predictors_names']] <- c( model[['predictors_names']], names(ras) )
          model[['predictors_types']] <- rbind.data.frame(model[['predictors_types']],
                                                          data.frame(predictors = names(ras),
                                                                     type = "numeric" )
          )
          if( !all(names(ras) %in% model[['predictors_object']]$get_names()) ){
            model[['predictors_object']]$data <- c(model[['predictors_object']]$data, ras)
          }

        } else if(m == "nnd") {
          # Nearest neighbour
          biodiversity_ids <- as.character( x$biodiversity$get_ids() )
          cc <- terra::rast()
          for(id in biodiversity_ids) {
            # Get presence points
            o <- guess_sf( x$biodiversity$get_data(id) )
            o <- o[ which( o[[x$biodiversity$get_columns_occ()[[id]]]] > 0 ), ]
            # Calculate point distance
            ras <- terra::distance(x = emptyraster( model$predictors_object$get_data() ),
                                   y = o)
            ras <- terra::mask(ras, model$background)
            names(ras) <- paste0("nearestpoint_", which(biodiversity_ids == id))
            suppressWarnings( cc <- c(cc, ras) )
            rm(ras, o )
          }
          # Add to predictor objects, names, types and the object
          model[['predictors']] <- cbind.data.frame( model[['predictors']], terra::as.data.frame(cc, na.rm = FALSE) )
          model[['predictors_names']] <- c( model[['predictors_names']], names(cc) )
          model[['predictors_types']] <- rbind.data.frame(model[['predictors_types']],
                                                          data.frame(predictors = names(cc),
                                                                     type = "numeric" )
          )
          if( !all(names(cc) %in% model[['predictors_object']]$get_names()) ){
            model[['predictors_object']]$data <- c(model[['predictors_object']]$data, cc)
          }
          rm(cc, biodiversity_ids)

        }
      }
    } else { model[["latent"]] <- new_waiver() }# End of latent factor loop

    # Set offset if existing
    if(!is.Waiver(x$offset)){
      # Aggregate offset if necessary
      if(terra::nlyr(x$offset)>1){
        # As log(x) + log(y) == log( x * y )
        ras_of <- sum(x$offset, na.rm = TRUE)
        # Normalize the result
        ras_of <- predictor_transform(ras_of, option = "norm")
        names(ras_of) <- "spatial_offset"
      } else {
        ras_of <- x$offset
        names(ras_of) <- "spatial_offset"
      }
      # Save overall offset
      ofs <- terra::as.data.frame(ras_of, xy = TRUE, na.rm = FALSE)
      names(ofs)[which(names(ofs)==names(ras_of))] <- "spatial_offset"
      model[['offset']] <- ofs
      # Also add offset object for faster extraction
      model[['offset_object']] <- ras_of
    } else { model[['offset']] <- new_waiver() }

    # Setting up variable bias control if set
    if(!is.Waiver( x$get_biascontrol())){
      if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding bias variable for bias control.')
      bias <- x$bias
      if(bias$method == "partial"){
        settings$set("bias_variable", names(bias$layer) )
        settings$set("bias_value", bias$bias_value )
        # Check that variable is already in the predictors object
        if(!(names(bias$layer) %in% model$predictors_names)){
          model$predictors_object <- model$predictors_object$set_data(names(bias$layer), bias$layer)
          # Also set predictor names
          model[['predictors_names']] <- model$predictors_object$get_names()
          model[['predictors']] <- model$predictors_object$get_data(df = TRUE, na.rm = FALSE)
          # Get predictor types
          lu <- sapply(model[['predictors']][model[['predictors_names']]], is.factor)
          model[['predictors_types']] <- data.frame(predictors = names(lu), type = ifelse(lu, 'factor', 'numeric') )
        }
        assertthat::assert_that(nrow(model[['predictors']]) == terra::ncell(model$predictors_object$get_data()))
      }
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
      model[['biodiversity']][[id]][['link']]         <- x$biodiversity$get_links()[[id]]
      model[['biodiversity']][[id]][['equation']]     <- x$biodiversity$get_equations()[[id]]
      model[['biodiversity']][[id]][['use_intercept']]<- x$biodiversity$data[[id]]$use_intercept # Separate intercept?
      model[['biodiversity']][[id]][['expect']]       <- x$biodiversity$get_weights()[[id]] # Weights per dataset
      # --- #
      # Rename observation column to 'observed'. Needs to be consistent for INLA
      # FIXME: try and not use dplyr as dependency (although it is probably loaded already)
      model$biodiversity[[id]]$observations <- model$biodiversity[[id]]$observations |> dplyr::rename('observed' = x$biodiversity$get_columns_occ()[[id]])
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
        # Check and reset multiplication weights
        if(nrow(o) != length( model[['biodiversity']][[id]][['expect']] )){
          if(length(unique( model[['biodiversity']][[id]][['expect']] ))>1){
            myLog('[Setup]','red', 'First weight is taken from the observations due to type conversion!')
          }
          val <- unique( model[['biodiversity']][[id]][['expect']] )[1]
          model[['biodiversity']][[id]][['expect']] <- rep(val, nrow(o))
        }
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

      # Aggregate observations if poipo
      if(aggregate_observations && model[['biodiversity']][[id]][['type']] == "poipo"){
        model$biodiversity[[id]]$observations <- aggregate_observations2grid(
          df = model$biodiversity[[id]]$observations,
          template = emptyraster(x$predictors$get_data(df = FALSE)),
          field_occurrence = "observed")
        # Check and reset multiplication weights
        model[['biodiversity']][[id]][['expect']] <- rep(unique( model[['biodiversity']][[id]][['expect']] )[1],
                                                         nrow(model$biodiversity[[id]]$observations))
      }

      # Now extract coordinates and extract estimates, shifted to raster extraction by default to improve speed!
      env <- get_rastervalue(coords = guess_sf(model$biodiversity[[id]]$observations),
                             env = model$predictors_object$get_data(df = FALSE),
                             rm.na = FALSE)

      # Remove missing values as several engines can't deal with those easily
      miss <- stats::complete.cases(env)
      if(sum( !miss )>0 && getOption('ibis.setupmessages')) {
        myLog('[Setup]','yellow', 'Excluded ', sum( !miss ), ' observations owing to missing values in covariates!' )
      }
      model[['biodiversity']][[id]][['observations']] <- model[['biodiversity']][[id]][['observations']][miss,]
      model[['biodiversity']][[id]][['expect']] <- model[['biodiversity']][[id]][['expect']][miss]
      env <- subset(env, miss)
      if(nrow(env)<=2) stop("Too many missing data points in covariates. Check out 'predictor_homogenize_na' and projections.")
      if( all( model[['biodiversity']][[id]][['observations']]$observed == 0) ) stop("All presence records fall outside the modelling background.")
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

      # Biodiversity dataset specific predictor refinement if the option is set
      if(settings$get("filter_predictors")!= "none"){
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','yellow', paste0('Filtering predictors via ',
                                                                                  settings$get("filter_predictors"),'...'))
        # Make backups
        test <- env;test$x <- NULL;test$y <- NULL;test$Intercept <- NULL

        # Ignore variables for which we have priors
        if(!is.Waiver(x$priors)){
          keep <- unique( as.character(x$priors$varnames()) )
          if('spde'%in% keep) keep <- keep[which(keep!='spde')] # Remove SPDE where existing
          test <- test[,-which(names(test) %in% keep)]
          assertthat::assert_that(!any(keep %in% names(test)))
        } else {keep <- NULL}
        # Add bias variable to keep as we risk filtering it out otherwise
        if(!is.Waiver(settings$get("bias_variable"))) keep <- c(keep, settings$get("bias_variable") )

        # Filter the predictors
        # Depending on the option this function returns the variables to be removed.
        co <- predictor_filter(env = test,
                               keep = keep,
                               cutoff = getOption('ibis.corPred'), # Probably keep default, but maybe sth. to vary in the future
                               method = settings$get("filter_predictors"),
                               observed = model[['biodiversity']][[id]]$observations[['observed']],
                               family = model[['biodiversity']][[id]]$family,
                               tune.type = "gic",
                               weight = NULL,
                               verbose = getOption('ibis.setupmessages')
                               )

        # For all factor variables, remove those with only the minimal value (e.g. 0)
        fac_min <- apply(test[,model$predictors_types$predictors[which(model$predictors_types$type=='factor')]], 2, function(x) min(x,na.rm = TRUE))
        fac_mean <- apply(test[,model$predictors_types$predictors[which(model$predictors_types$type=='factor')]], 2, function(x) mean(x,na.rm = TRUE))
        co <- unique(co, names(which(fac_mean == fac_min)) ) # Now add to co all those variables where the mean equals the minimum, indicating only absences
        # Remove variables if found
        if(length(co)>0){
          env |> dplyr::select(-dplyr::all_of(co)) -> env
        }

      } else { co <- NULL }

      # Save predictors extracted for biodiversity extraction
      model[['biodiversity']][[id]][['predictors']] <- env
      model[['biodiversity']][[id]][['predictors_names']] <- model[['predictors_names']][which( model[['predictors_names']] %notin% co )]
      model[['biodiversity']][[id]][['predictors_types']] <- model[['predictors_types']][model[['predictors_types']]$predictors %notin% co,]
  }

    # If the method of integration is weights and there are more than 2 datasets, combine
    if(method_integration == "weight" && length(model$biodiversity)>=2){
      if(getOption('ibis.setupmessages')) myLog('[Setup]','yellow','Experimental: Integration by weights assumes identical data parameters!')
      # Check that all types and families can be combined
      types <- as.character( sapply( model$biodiversity, function(x) x$type ) )
      fams <- as.character( sapply( model$biodiversity, function(z) z$family ) )
      assertthat::assert_that(length(unique(types))==1, length(unique(fams))==1,
                              msg = "Integration by weights requires identical biodiversity datasets!")
      obs <- lapply( model$biodiversity, function(x) {
        guess_sf( x$observations )
      } )
      obs <- do.call("rbind", obs)
      w <- lapply( model$biodiversity, function(x) x$expect ) |> unlist() |> unname()
      assertthat::assert_that(nrow(obs) == length(w))
      preds <- lapply( model$biodiversity, function(x) x$predictors )
      preds <- do.call("rbind", preds) |> unique()
      predn <- lapply( model$biodiversity, function(x) x$predictors_names ) |> unlist() |> unname() |> unique()
      predt <- lapply( model$biodiversity, function(x) x$predictors_types )
      predt <- do.call("rbind", predt) |> unique()
      # Now combine the biodiversity objects and create a new id
      new <- list(
        name = "Combined_data_weight",
        observations = obs |> sf::st_drop_geometry(),
        type = unique(types)[1], family = unique(fams)[1],
        equation = "<Default>", # Use default equation #FIXME This could be more cleverer
        use_intercept = TRUE, # Assume default
        expect = w,
        predictors = preds, predictors_names = predn, predictors_types = predt
      )
      model[['biodiversity']] <- list()
      model[['biodiversity']][[as.character(new_id())]] <- new
      rm(new, obs, w, preds, predn, predt)
    }

    # Warning if Np is larger than Nb
    if(settings$get("filter_predictors") == "none"){
      if( sum(x$biodiversity$get_observations() )-1 <= length(model$predictors_names)){
        if(getOption('ibis.setupmessages')) myLog('[Setup]','red', 'More predictors than observations! Consider settings optim_hyperparam or filter_predictors!')
      }
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
        "<GLMNET>" = x$priors$classes() == "GLMNETPrior",
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
      coords <- subset(coords, observed > 0) # Remove absences
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

      # Only if some points actually fall in the zones
      if(nrow(zones)>0){
        # Now clip all predictors and background to this
        model$background <- suppressMessages(
          suppressWarnings( sf::st_union( sf::st_intersection(zones, model$background), by_feature = TRUE) |>
                              sf::st_buffer(dist = 0)  |> # 0 distance buffer trick
                              sf::st_cast("MULTIPOLYGON")
          )
        )

        # Extract predictors and offsets again if set
        if(!is.Waiver(model$predictors_object)){
          # Using the raster operations is generally faster than point in polygon tests
          pred_ov <- model$predictors_object$get_data(df = FALSE)
          # Make a rasterized mask of the background
          pred_ov <- terra::mask( pred_ov, model$background )
          # Convert Predictors to data.frame, including error catching for raster errors
          # FIXME: This could be outsourced
          o <- try({ terra::as.data.frame(pred_ov, xy = TRUE, na.rm = FALSE) },silent = TRUE)
          if(inherits(o, "try-error")){
            o <- as.data.frame( cbind( terra::crds(pred_ov),
                                       as.matrix( pred_ov )) )
            if(any(is.factor(pred_ov))){
              o[names(pred_ov)[which(is.factor(pred_ov))]] <- factor(o[names(pred_ov)[which(is.factor(pred_ov))]] )
            }
          }
          model[['predictors']] <- o
          model[['predictors_object']]$data <- fill_rasters(o[,c(1,2)*-1], # Remove x and y coordinates for overwriting raster data
                                                            model$predictors_object$data)
          rm(pred_ov, o)
        } else {
          model$predictors[which( is.na(
            point_in_polygon(poly = model$background, points = model$predictors[,c('x','y')] )[['limit']]
          )),model$predictors_names] <- NA # Fill with NA
        }
        # The same with offset if specified, Note this operation below is computationally quite costly
        # MJ: 18/10/22 Removed below as (re)-extraction further in the pipeline makes this step irrelevant
        # if(!is.Waiver(x$offset)){
        #   model$offset[which( is.na(
        #     point_in_polygon(poly = zones, points = model$offset[,c('x','y')] )[['limit']]
        #   )), "spatial_offset" ] <- NA # Fill with NA
        # }
      }
    }
    # Messenger
    if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Adding engine-specific parameters.')

    # Basic consistency checks
    assertthat::assert_that(nrow(model$biodiversity[[1]]$observations)>0,
                            length(model[['biodiversity']][[1]][['expect']])>1,
                            all(c("predictors","background","biodiversity") %in% names(model) ),
                            length(model$biodiversity[[1]]$expect) == nrow(model$biodiversity[[1]]$predictors)
    )
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

      # Create the mesh if not already present
      x$engine$create_mesh(model = model)
      assertthat::assert_that(inherits(x$engine$get_data("mesh"), "inla.mesh"),
                              msg = "Something went wrong during mesh creation...")

      # If set specify a SPDE effect
      if((!is.Waiver(x$latentfactors))){
        if(attr(x$get_latent(),'method') == "spde"){
          x$engine$calc_latent_spatial(type = attr(x$get_latent(),'method'), priors = model[['priors']])
        }
      }

      # Process per supplied dataset
      for(id in ids) {

        # Update model formula in the model container
        model$biodiversity[[id]]$equation <- built_formula_inla(model = model,
                                                                id = id,
                                                                x = x,
                                                                settings = settings)

        # For each type include expected data
        # expectation vector (area for integration points/nodes and 0 for presences)
        if(model$biodiversity[[id]]$family == 'poisson') model$biodiversity[[id]][['expect']] <- rep(0, nrow(model$biodiversity[[id]]$predictors) )
        if(model$biodiversity[[id]]$family == 'binomial') model$biodiversity[[id]][['expect']] <- rep(1, nrow(model$biodiversity[[id]]$predictors) ) * model$biodiversity[[id]]$expect
      }
      # Run the engine setup script
      model <- x$engine$setup(model, settings)

      # Now train the model and create a predicted distribution model
      out <- x$engine$train(model, settings)

      # ----------------------------------------------------------- #
      #### INLABRU Engine ####
    } else if( inherits(x$engine,'INLABRU-Engine') ){

      # Create the mesh if not already present
      x$engine$create_mesh(model = model)
      assertthat::assert_that(inherits(x$engine$get_data("mesh"), "inla.mesh"),
                              msg = "Something went wrong during mesh creation...")

      # If set specify a SPDE effect
      if((!is.Waiver(x$latentfactors))){
        if(attr(x$get_latent(),'method') == "spde"){
          x$engine$calc_latent_spatial(type = attr(x$get_latent(),'method'), priors = model[['priors']])
        }
      }

      # Process per supplied dataset
      for(id in ids) {

        # Update model formula in the model container
        model$biodiversity[[id]]$equation <- built_formula_inla(model = model,
                                                                id = id,
                                                                x = x,
                                                                settings = settings)
        # For each type include expected data
        # expectation vector (area for integration points/nodes and 0 for presences)
        if(model$biodiversity[[id]]$family == 'poisson') model$biodiversity[[id]][['expect']] <- rep(0, nrow(model$biodiversity[[id]]$predictors) )
        if(model$biodiversity[[id]]$family == 'binomial') model$biodiversity[[id]][['expect']] <- rep(1, nrow(model$biodiversity[[id]]$predictors) ) * model$biodiversity[[id]]$expect
      }

      # Run the engine setup script
      x$engine$setup(model, settings)

      # Now train the model and create a predicted distribution model
      out <- x$engine$train(model, settings)

      # ----------------------------------------------------------- #
      #### GDB Engine ####
    } else if( inherits(x$engine,"GDB-Engine") ){

      # For each formula, process in sequence
      for(id in ids){

        model$biodiversity[[id]]$equation <- built_formula_gdb( model = model,
                                                                id = id,
                                                                x = x,
                                                                settings = settings)

        # Remove those not part of the modelling
        model2 <- model
        model2$biodiversity <- NULL; model2$biodiversity[[id]] <- model$biodiversity[[id]]

        # Run the engine setup script
        model2 <- x$engine$setup(model2, settings)

        # Now train the model and create a predicted distribution model
        settings2 <- settings
        if(id != ids[length(ids)] && method_integration == "prior") {
          # No need to make predictions if we use priors only
          settings2$set('inference_only', TRUE)
        } else if(id != ids[length(ids)]){
          # For predictors and offsets
          settings2$set('inference_only', FALSE)
        } else {
          settings2$set('inference_only', inference_only)
        }
        out <- x$engine$train(model2, settings2)
        rm(model2)
        # Add Prediction of model to next object if multiple are supplied
        if(length(ids)>1 && id != ids[length(ids)]){
          if(method_integration == "predictor"){
            # Add to predictors frame
            new <- out$get_data("prediction")
            pred_name <- paste0(model$biodiversity[[id]]$type, "_", make.names(model$biodiversity[[id]]$name),"_mean")
            names(new) <- pred_name

            # Add the object to the overall prediction object
            model$predictors_object$data <- c(model$predictors_object$get_data(), new)

            # Now for each biodiversity dataset and the overall predictors
            # extract and add as variable
            for(k in names(model$biodiversity)){
              env <- get_rastervalue(coords = guess_sf(model$biodiversity[[k]]$observations[,c('x','y')]),
                                     env = new)
              # Rename to current id dataset
              env <- env[names(new)]
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
            model$predictors <- cbind(model$predictors, as.data.frame(new, na.rm = FALSE) )
            model$predictors_names <- c(model$predictors_names, names(new))
            model$predictors_types <- rbind(model$predictors_types,
                                            data.frame(predictors = names(new), type = c('numeric')))
          } else if(method_integration == "offset"){
            # Adding the prediction as offset
            new <- out$get_data("prediction")
            # Back transforming offset to linear scale
            new[] <- switch (model$biodiversity[[id]]$family,
                             "binomial" = ilink(new[], link = "logit"),
                             "poisson" = ilink(new[], link = "log")
            )
            if(is.Waiver(model$offset)){
              ofs <- terra::as.data.frame(new, xy = TRUE, na.rm = FALSE)
              names(ofs)[which(names(ofs)==names(new))] <- "spatial_offset"
              model[['offset']] <- ofs
              # Also add offset object for faster extraction
              model[['offset_object']] <- new
            } else {
              # New offset
              news <- sum( model[['offset_object']], new, na.rm = TRUE)
              news <- terra::mask(news, x$background)
              model[['offset_object']] <- news
              ofs <- terra::as.data.frame(news, xy = TRUE, na.rm = FALSE)
              names(ofs)[which(names(ofs)=="layer")] <- "spatial_offset"
              model[['offset']] <- ofs
              rm(news)
            }
            rm(new)
          } else if(method_integration == "prior"){
            # Use the previous model to define and set priors
            po <- get_priors(out, x$engine$name)
            model$priors <- po
          }
        } # End of multiple ids
      }
      # ----------------------------------------------------------- #
      #### XGBoost Engine ####
    } else if( inherits(x$engine,"XGBOOST-Engine") ){
      # Create XGBboost regression and classification

      # TODO: Combine biodiversity datasets and add factor variable
      # Ideally figure out a convenient way to allow interactions. Maybe by just multiplying all predictors?
      if(method_integration == "interaction") stop("Not yet implemented")

      # Process per supplied dataset and in order of supplied data
      for(id in ids) {

        # Update model formula in the model container
        model$biodiversity[[id]]$equation <- built_formula_xgboost( model$biodiversity[[id]] )

        # Remove those not part of the modelling
        model2 <- model
        model2$biodiversity <- NULL; model2$biodiversity[[id]] <- model$biodiversity[[id]]

        # Run the engine setup script
        model2 <- x$engine$setup(model2, settings)

        # Now train the model and create a predicted distribution model
        settings2 <- settings
        if(id != ids[length(ids)] && method_integration == "prior") {
          # No need to make predictions if we use priors only
          settings2$set('inference_only', TRUE)
        } else if(id != ids[length(ids)]){
          # For predictors and offsets
          settings2$set('inference_only', FALSE)
        } else {
          settings2$set('inference_only', inference_only)
        }
        out <- x$engine$train(model2, settings2)
        rm(model2)

        # Add Prediction of model to next object if multiple are supplied
        if(length(ids)>1 && id != ids[length(ids)]){
          # Adding the prediction to the model object
          if(method_integration == "predictor"){
            # Add to predictors frame
            new <- out$get_data("prediction")
            pred_name <- paste0(model$biodiversity[[id]]$type, "_", make.names(model$biodiversity[[id]]$name),"_mean")
            names(new) <- pred_name
            # Add the object to the overall prediction object
            model$predictors_object$data <- c(model$predictors_object$get_data(), new)
            # Now for each biodiversity dataset and the overall predictors
            # extract and add as variable
            for(k in names(model$biodiversity)){
              env <- get_rastervalue(coords = guess_sf(model$biodiversity[[k]]$observations[,c('x','y')]),
                                     env = new)
              # Rename to current id dataset
              env <- env[names(new)]
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
            model$predictors <- cbind(model$predictors, as.data.frame(new, na.rm = FALSE) )
            model$predictors_names <- c(model$predictors_names, names(new))
            model$predictors_types <- rbind(model$predictors_types,
                                            data.frame(predictors = names(new), type = c('numeric')))

          } else if(method_integration == "offset"){
            # Adding the prediction as offset
            new <- out$get_data("prediction")
            # Back transforming offset to linear scale
            new[] <- switch (model$biodiversity[[id]]$family,
                             "binomial" = ilink(new[], link = "logit"),
                             "poisson" = ilink(new[], link = "log")
            )
            if(is.Waiver(model$offset)){
              ofs <- terra::as.data.frame(new, xy = TRUE, na.rm = FALSE)
              names(ofs)[which(names(ofs)==names(new))] <- "spatial_offset"
              model[['offset']] <- ofs
              # Also add offset object for faster extraction
              model[['offset_object']] <- new
            } else {
              # New offset
              news <- sum( model[['offset_object']], new, na.rm = TRUE)
              news <- terra::mask(news, x$background)
              model[['offset_object']] <- news
              ofs <- terra::as.data.frame(news, xy = TRUE, na.rm = FALSE)
              names(ofs)[which(names(ofs)=="layer")] <- "spatial_offset"
              model[['offset']] <- ofs
              rm(news)
            }
            rm(new)
          } else if(method_integration == "prior"){
            # Use the previous model to define and set priors
            po <- get_priors(out, x$engine$name)
            model$priors <- po
          }
        }
      }

      # ----------------------------------------------------------- #
      #### BART Engine ####
    } else if( inherits(x$engine,"BART-Engine") ){

      # TODO: Combine biodiversity datasets and add factor variable
      # Ideally figure out a convenient way to allow interactions. Maybe by just multiplying all predictors?
      if(method_integration == "interaction") stop("Not yet implemented")

      # Process each id
      for(id in ids){

        # Built formula
        model$biodiversity[[id]]$equation <- built_formula_bart( model$biodiversity[[id]] )

        # Model and estimate
        model2 <- model
        model2$biodiversity <- NULL; model2$biodiversity[[id]] <- model$biodiversity[[id]]

        # Run the engine setup script
        model2 <- x$engine$setup(model2, settings)

        # Now train the model and create a predicted distribution model
        settings2 <- settings
        if(id != ids[length(ids)] && method_integration == "prior") {
          # No need to make predictions if we use priors only
          settings2$set('inference_only', TRUE)
        } else if(id != ids[length(ids)]){
          # For predictors and offsets
          settings2$set('inference_only', FALSE)
        } else {
          settings2$set('inference_only', inference_only)
        }
        out <- x$engine$train(model2, settings2)
        rm(model2)

        # Add Prediction of model to next object if multiple are supplied
        if(length(ids)>1 && id != ids[length(ids)]){
          if(method_integration == "predictor"){
            # Add to predictors frame
            new <- out$get_data("prediction")[[1]] # Only take the mean!
            pred_name <- paste0(model$biodiversity[[id]]$type, "_", make.names(model$biodiversity[[id]]$name),"_mean")
            names(new) <- pred_name
            # Add the object to the overall prediction object
            model$predictors_object$data <- c(model$predictors_object$get_data(), new)

            # Now for each biodiversity dataset and the overall predictors
            # extract and add as variable
            for(k in names(model$biodiversity)){
              env <- get_rastervalue(coords = guess_sf(model$biodiversity[[k]]$observations[,c('x','y')]),
                                     env = new)
              # Rename to current id dataset
              env <- env[names(new)]

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
            model$predictors <- cbind(model$predictors, as.data.frame(new, na.rm = FALSE) )
            model$predictors_names <- c(model$predictors_names, names(new))
            model$predictors_types <- rbind(model$predictors_types,
                                            data.frame(predictors = names(new), type = c('numeric')))
          } else if(method_integration == "offset"){
            # Adding the prediction as offset
            new <- out$get_data("prediction")[["mean"]]
            # Back transforming offset to linear scale
            new[] <- switch (model$biodiversity[[id]]$family,
                             "binomial" = ilink(new[], link = "logit"),
                             "poisson" = ilink(new[], link = "log")
            )
            names(new) <- "spatial_offset"

            if(is.Waiver(model$offset)){
              ofs <- terra::as.data.frame(new, xy = TRUE, na.rm = FALSE)
              names(ofs)[which(names(ofs)==names(new))] <- "spatial_offset"
              model[['offset']] <- ofs
              # Also add offset object for faster extraction
              model[['offset_object']] <- new
            } else {
              # New offset
              news <- sum( model[['offset_object']], new, na.rm = TRUE)
              news <- terra::mask(news, x$background)
              names(news) <- "spatial_offset"
              model[['offset_object']] <- news
              ofs <- terra::as.data.frame(news, xy = TRUE, na.rm = FALSE)
              names(ofs)[which(names(ofs)=="layer")] <- "spatial_offset"
              model[['offset']] <- ofs
              rm(news)
            }
            rm(new)
          } else if(method_integration == "prior"){
            # Use the previous model to define and set priors
            po <- get_priors(out, x$engine$name)
            model$priors <- po
          }
        } # End of multiple likelihood function

      } # End of id loop
    } else if( inherits(x$engine,"STAN-Engine") ){
      # ----------------------------------------------------------- #
      #### STAN Engine ####
      # Process per supplied dataset
      for(id in ids) {
        # TODO
        if(length(model$biodiversity)>1) stop("Not yet implemented")

        # Update model formula in the model container
        model$biodiversity[[id]]$equation <- built_formula_stan(model = model,
                                                                id = id,
                                                                x = x,
                                                                settings = settings)
      }

      # Run the engine setup script
      model <- x$engine$setup(model, settings)

      # Now train the model and create a predicted distribution model
      out <- x$engine$train(model, settings)


    } else if (inherits(x$engine,"BREG-Engine") ){
      # ----------------------------------------------------------- #
      #### BREG Engine ####
      assertthat::assert_that(
        !(method_integration == "offset" && any(types == "poipa")),
        msg = "Due to engine limitations BREG models do not support offsets for presence-absence models!"
      )
      # For each formula, process in sequence
      for(id in ids){

        model$biodiversity[[id]]$equation <- built_formula_breg( model$biodiversity[[id]] )

        # Remove those not part of the modelling
        model2 <- model
        model2$biodiversity <- NULL; model2$biodiversity[[id]] <- model$biodiversity[[id]]

        # Run the engine setup script
        model2 <- x$engine$setup(model2, settings)

        # Now train the model and create a predicted distribution model
        settings2 <- settings
        if(id != ids[length(ids)] && method_integration == "prior") {
          # No need to make predictions if we use priors only
          settings2$set('inference_only', TRUE)
        } else if(id != ids[length(ids)]){
          # For predictors and offsets
          settings2$set('inference_only', FALSE)
        } else {
          settings2$set('inference_only', inference_only)
        }
        out <- x$engine$train(model2, settings2)

        # Add Prediction of model to next object if multiple are supplied
        if(length(ids)>1 && id != ids[length(ids)]){
          if(method_integration == "predictor"){
            # Add to predictors frame
            new <- out$get_data("prediction")[["mean"]]
            pred_name <- paste0(model$biodiversity[[id]]$type, "_", make.names(model$biodiversity[[id]]$name),"_mean")
            names(new) <- pred_name
            # Add the object to the overall prediction object
            model$predictors_object$data <- c(model$predictors_object$get_data(), new)

            # Now for each biodiversity dataset and the overall predictors
            # extract and add as variable
            for(k in names(model$biodiversity)){
              env <- get_rastervalue(coords = guess_sf(model$biodiversity[[k]]$observations[,c('x','y')]),
                                     env = new)
              # Rename to current id dataset
              env <- env[names(new)]
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
            model$predictors <- cbind(model$predictors, as.data.frame(new, na.rm = FALSE))
            model$predictors_names <- c(model$predictors_names, names(new))
            model$predictors_types <- rbind(model$predictors_types,
                                            data.frame(predictors = names(new), type = c('numeric')))

          } else if(method_integration == "offset"){
            # Adding the prediction as offset
            new <- out$get_data("prediction")
            # Back transforming offset to linear scale
            new[] <- switch (model$biodiversity[[id]]$family,
                             "binomial" = ilink(new[], link = "logit"),
                             "poisson" = ilink(new[], link = "log")
            )
            if(is.Waiver(model$offset)){
              ofs <- terra::as.data.frame(new, xy = TRUE, na.rm = FALSE)
              names(ofs)[which(names(ofs)==names(new))] <- "spatial_offset"
              model[['offset']] <- ofs
              # Also add offset object for faster extraction
              model[['offset_object']] <- new
            } else {
              # New offset
              news <- sum( model[['offset_object']], new, na.rm = TRUE)
              news <- terra::mask(news, x$background)
              model[['offset_object']] <- news
              ofs <- terra::as.data.frame(news, xy = TRUE, na.rm = FALSE)
              names(ofs)[which(names(ofs)=="layer")] <- "spatial_offset"
              model[['offset']] <- ofs
              rm(news)
            }
            rm(new)
          } else if(method_integration == "prior"){
            # Use the previous model to define and set priors
            po <- get_priors(out, x$engine$name)
            model$priors <- po
          }

        } # End of multiple ides
      }
    } else if (inherits(x$engine,"GLMNET-Engine") ){
      # ----------------------------------------------------------- #
      #### GLMNET Engine ####
      # For each formula, process in sequence
      for(id in ids){

        model$biodiversity[[id]]$equation <- built_formula_glmnet( model$biodiversity[[id]] )

        # Remove those not part of the modelling
        model2 <- model
        model2$biodiversity <- NULL; model2$biodiversity[[id]] <- model$biodiversity[[id]]

        # Run the engine setup script
        model2 <- x$engine$setup(model2, settings)

        # Now train the model and create a predicted distribution model
        settings2 <- settings
        if(id != ids[length(ids)] && method_integration == "prior") {
          # No need to make predictions if we use priors only
          settings2$set('inference_only', TRUE)
        } else if(id != ids[length(ids)]){
          # For predictors and offsets
          settings2$set('inference_only', FALSE)
        } else {
          settings2$set('inference_only', inference_only)
        }
        out <- x$engine$train(model2, settings2)

        # Add Prediction of model to next object if multiple are supplied
        if(length(ids)>1 && id != ids[length(ids)]){
          if(method_integration == "predictor"){
            # Add to predictors frame
            new <- out$get_data("prediction")[["mean"]]
            pred_name <- paste0(model$biodiversity[[id]]$type, "_", make.names(model$biodiversity[[id]]$name),"_mean")
            names(new) <- pred_name
            # Add the object to the overall prediction object
            model$predictors_object$data <- c(model$predictors_object$get_data(), new)

            # Now for each biodiversity dataset and the overall predictors
            # extract and add as variable
            for(k in names(model$biodiversity)){
              env <- get_rastervalue(coords =  guess_sf( model$biodiversity[[k]]$observations[,c('x','y')]),
                                     env = new)
              env <- env[names(new)]
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
            model$predictors <- cbind(model$predictors, as.data.frame(new, na.rm = FALSE))
            model$predictors_names <- c(model$predictors_names, names(new))
            model$predictors_types <- rbind(model$predictors_types,
                                            data.frame(predictors = names(new), type = c('numeric')))

          } else if(method_integration == "offset"){
            # Adding the prediction as offset
            new <- out$get_data("prediction")
            # Back transforming offset to linear scale
            new[] <- switch (model$biodiversity[[id]]$family,
                             "binomial" = ilink(new[], link = "logit"),
                             "poisson" = ilink(new[], link = "log")
            )
            if(is.Waiver(model$offset)){
              ofs <- terra::as.data.frame(new, xy = TRUE, na.rm = FALSE)
              names(ofs)[which(names(ofs)==names(new))] <- "spatial_offset"
              model[['offset']] <- ofs
              # Also add offset object for faster extraction
              model[['offset_object']] <- new
            } else {
              # New offset
              news <- sum( model[['offset_object']], new, na.rm = TRUE)
              news <- terra::mask(news, x$background)
              model[['offset_object']] <- news
              ofs <- terra::as.data.frame(news, xy = TRUE, na.rm = FALSE)
              names(ofs)[which(names(ofs)=="layer")] <- "spatial_offset"
              model[['offset']] <- ofs
              rm(news)
            }
            rm(new)
          } else if(method_integration == "prior"){
            # Use the previous model to define and set priors
            po <- get_priors(out, x$engine$name)
            model$priors <- po
          }
        } # End of multiple ides
      }
      # End of GLMNET engine
    } else { stop('Specified Engine not implemented yet.') }

    if(is.null(out)) return(NULL)

    if(getOption('ibis.setupmessages')) myLog('[Done]','green',paste0('Completed after ', round( as.numeric(out$settings$duration()), 2),' ',attr(out$settings$duration(),'units') ))

    # Clip to limits again to be sure
    if(!is.Waiver(x$limits)) {
      if(settings$get('inference_only')==FALSE){
        out <- out$set_data("prediction", terra::mask(out$get_data("prediction"), model$background))
      }
      out$settings$set("has_limits", TRUE)
    } else {
      out$settings$set("has_limits", FALSE)
    }

    # Stop logging if specified
    if(!is.Waiver(x$log)) x$log$close()

    # Return created object
    return(out)
  }
)
