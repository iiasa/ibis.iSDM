#' @include utils.R bdproto.R bdproto-biodiversitydistribution.R bdproto-predictors.R bdproto-biodiversityscenario.R
NULL

#' Add predictors to a Biodiversity distribution object
#'
#' @description
#' This function allows to add predictors to [distribution] or [BiodiversityScenario]
#' objects. Predictors are covariates that in spatial projection have to match
#' the geographic projection of the background layer in the [distribution] object.
#' This function furthermore allows to transform or create derivates of provided
#' predictors.
#'
#' A transformation takes the provided rasters and for instance rescales them or transforms
#' them through a principal component analysis ([prcomp]). In contrast, derivates leave
#' the original provided predictors alone, but instead create new ones, for instance by transforming
#' their values through a quadratic or hinge transformation. Note that this effectively
#' increases the number of predictors in the object, generally requiring stronger regularization by
#' the used [`engine`].
#' Both transformations and derivates can also be combined.
#' Available options for transformation are:
#' * \code{'none'} - Leaves the provided predictors in the original scale.
#' * \code{'pca'} - Converts the predictors to principal components. Note that this
#' results in a renaming of the variables to principal component axes!
#' * \code{'scale'} - Transforms all predictors by applying [scale] on them.
#' * \code{'norm'} - Normalizes all predictors by transforming them to a scale from 0 to 1.
#' * \code{'windsor'} - Applies a windsorization to the target predictors. By default
#' this effectively cuts the predictors to the 0.05 and 0.95, thus helping to remove
#' extreme outliers.
#'
#' Available options for creating derivates are:
#' * \code{'none'} - No additional predictor derivates are created.
#' * \code{'quad'} - Adds quadratic transformed predictors.
#' * \code{'interaction'} - Add interacting predictors. Interactions need to be specified (\code{"int_variables"})!
#' * \code{'thresh'} - Add threshold transformed predictors.
#' * \code{'hinge'} - Add hinge transformed predictors.
#' * \code{'bin'} - Add predictors binned by their percentiles.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param env A [`RasterStack-class`], [`RasterLayer-class`] or [`stars`] object.
#' @param names A [`vector`] of character names describing the environmental stack in case they should be renamed.
#' @param transform A [`vector`] stating whether predictors should be preprocessed in any way (Options: \code{'none'},\code{'pca'}, \code{'scale'}, \code{'norm'})
#' @param derivates A Boolean check whether derivate features should be considered (Options: \code{'none'}, \code{'thresh'}, \code{'hinge'}, \code{'quad'}) )
#' @param bgmask Check whether the environmental data should be masked with the background layer (Default: \code{TRUE})
#' @param harmonize_na A [`logical`] value indicating of whether NA values should be harmonized among predictors (Default: \code{FALSE})
#' @param explode_factors [`logical`] of whether any factor variables should be split up into binary variables (one per class). (Default: \code{FALSE}).
#' @param int_variables A [`vector`] with length greater or equal than \code{2} specifying the covariates  (Default: \code{NULL}).
#' @param priors A [`PriorList-class`] object. Default is set to \code{NULL} which uses default prior assumptions.
#' @param ... Other parameters passed down
#' @note
#' **Important:**
#' Not every [`engine`] supported by the \pkg{ibis.iSDM} R-package allows missing data points
#' among extracted covariates. Thus any observation with missing data is generally removed prior
#' from model fitting. Thus ensure that covariates have appropriate no-data settings (for instance setting \code{NA}
#' values to \code{0} or another out of range constant).
#'
#' Not every engine does actually need covariates. For instance it is perfectly legit
#' to fit a model with only occurrence data and a spatial latent effect ([add_latent]).
#' This correspondents to a spatial kernel density estimate.
#'
#' Certain names such \code{"offset"} are forbidden as predictor variable names. The function
#' will return an error message if these are used.
#' @aliases add_predictors
#' @examples
#' \dontrun{
#'  obj <- distribution(background) %>%
#'         add_predictors(covariates, transform = 'scale')
#'  obj
#' }
#' @name add_predictors
NULL

#' @name add_predictors
#' @rdname add_predictors
#' @exportMethod add_predictors
#' @export
methods::setGeneric(
  "add_predictors",
  signature = methods::signature("x", "env"),
  function(x, env, names = NULL, transform = 'scale', derivates = 'none', bgmask = TRUE,
           harmonize_na = FALSE, explode_factors = FALSE, int_variables = NULL, priors = NULL, ...) standardGeneric("add_predictors"))

#' @name add_predictors
#' @rdname add_predictors
#' @usage \S4method{add_predictors}{BiodiversityDistribution,RasterBrick}(x, env)
methods::setMethod(
  "add_predictors",
  methods::signature(x = "BiodiversityDistribution", env = "RasterBrick"),
  function(x, env, names = NULL, transform = 'scale', derivates = 'none', bgmask = TRUE, harmonize_na = FALSE,
           explode_factors = FALSE, int_variables = NULL, priors = NULL, ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            !missing(env))
    # Convert env to stack if it is a single layer only
    env = raster::stack(env)
    add_predictors(x, env, names, transform, derivates, bgmask, harmonize_na, explode_factors, int_variables, priors, ...)
  }
)

#' @name add_predictors
#' @rdname add_predictors
#' @usage \S4method{add_predictors}{BiodiversityDistribution,RasterLayer}(x, env)
methods::setMethod(
  "add_predictors",
  methods::signature(x = "BiodiversityDistribution", env = "RasterLayer"),
  function(x, env, names = NULL, transform = 'scale', derivates = 'none', bgmask = TRUE,
           harmonize_na = FALSE, explode_factors = FALSE, int_variables = NULL, priors = NULL, ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            !missing(env))
    # Convert env to stack if it is a single layer only
    env = raster::stack(env)
    add_predictors(x, env, names, transform, derivates, bgmask, harmonize_na, explode_factors, int_variables, priors, ...)
  }
)

# TODO: Support other objects other than Raster stacks such as data.frames and stars objects
#' @name add_predictors
#' @rdname add_predictors
#' @usage \S4method{add_predictors}{BiodiversityDistribution,RasterStack}(x, env)
methods::setMethod(
  "add_predictors",
  methods::signature(x = "BiodiversityDistribution", env = "RasterStack"),
  function(x, env, names = NULL, transform = 'scale', derivates = 'none', bgmask = TRUE,
           harmonize_na = FALSE, explode_factors = FALSE, int_variables = NULL, priors = NULL, ... ) {
    # Try and match transform and derivatives arguments
    transform <- match.arg(transform, c('none','pca', 'scale', 'norm', 'windsor') , several.ok = TRUE)
    derivates <- match.arg(derivates, c('none','thresh', 'hinge', 'quadratic', 'bin', 'interaction') , several.ok = TRUE)

    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.Raster(env),
                            all(transform == 'none') || all( transform %in% c('pca', 'scale', 'norm', 'windsor') ),
                            all(derivates == 'none') || all( derivates %in% c('thresh', 'hinge', 'quadratic', 'bin', 'interaction') ),
                            is.null(names) || assertthat::is.scalar(names) || is.vector(names),
                            is.logical(explode_factors),
                            is.null(priors) || inherits(priors,'PriorList')
    )
    assertthat::assert_that(sf::st_crs(x$background) == sf::st_crs(env@crs),
                            msg = 'Supplied environmental data not aligned with background.')
    # Messenger
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding predictors...')

    if(!is.null(names)) {
      assertthat::assert_that(nlayers(env)==length(names),
                              all(is.character(names)),
                                                msg = 'Provided names not of same length as environmental data.')
      # Set names of env
      names(env) <- names
    }

    # Check that all names allowed
    problematic_names <- grep("offset|w|weight|spatial_offset|Intercept|spatial.field", names(env),fixed = TRUE)
    if( length(problematic_names)>0 ){
      stop(paste0("Some predictor names are not allowed as they might interfere with model fitting:", paste0(names(env)[problematic_names],collapse = " | ")))
    }

    # If priors have been set, save them in the distribution object
    if(!is.null(priors)) {
      assertthat::assert_that( all( priors$varnames() %in% names(env) ) )
      x <- x$set_priors(priors)
    }
    # Harmonize NA values
    if(harmonize_na){
      if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Harmonizing missing values...')
      env <- predictor_homogenize_na(env, fill = FALSE)
    }

    # Don't transform or create derivatives of factor variables
    if(any(is.factor(env))){
      # Make subsets to join back later
      env_f <- raster::subset(env, which(is.factor(env)))
      env <- raster::subset(env, which(!is.factor(env)))
      if(explode_factors){
        # Refactor categorical variables
        if(inherits(env_f,'RasterLayer')){
          env_f <- explode_factorized_raster(env_f)
          env <- addLayer(env, env_f)
        } else {
          o <- raster::stack()
          for(layer in names(env_f)){
            o <- raster::addLayer(o, explode_factorized_raster(env_f[[layer]]))
          }
          env_f <- o;rm(o)
          # Joining back to full raster stack
          env <- raster::stack(env, env_f);rm(env_f)
        }
        has_factors <- FALSE # Set to false since factors have been exploded.
      } else { has_factors <- TRUE }
    } else { has_factors <- FALSE }

    # Standardization and scaling
    if('none' %notin% transform){
      if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Transforming predictors...')
      for(tt in transform) env <- predictor_transform(env, option = tt)
    }

    # Calculate derivates if set
    if('none' %notin% derivates){
      if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Creating predictor derivates...')
      # Specific condition for interaction
      if(any(derivates == "interaction")){
        assertthat::assert_that(is.vector(int_variables), length(int_variables)>=2)
        attr(env, "int_variables") <- int_variables
      }
      new_env <- raster::stack()
      for(dd in derivates) new_env <- raster::addLayer(new_env, predictor_derivate(env, option = dd) )

      # Add to env
      env <- raster::addLayer(env, new_env)
    }

    # Add factors back in if there are any.
    # This is to avoid that they are transformed or similar
    if(has_factors){
      env <- raster::addLayer(env, env_f)
    }
    attr(env, 'has_factors') <- has_factors

    # Assign an attribute to this object to keep track of it
    attr(env,'transform') <- transform

    # Mask predictors with existing background layer
    if(bgmask){
      env <- raster::mask(env, mask = x$background)
      # Reratify, work somehow only on stacks
      if(has_factors && any(is.factor(env)) ){
        new_env <- raster::stack(env)
        new_env[[which(is.factor(env))]] <- raster::ratify(env[[which(is.factor(env))]])
        env <- new_env;rm(new_env)
      } else env <- raster::stack(env)
    }

    # Check whether predictors already exist, if so overwrite
    if(!is.Waiver(x$predictors)) myLog('[Setup]','yellow','Overwriting existing predictors.')

    # Finally set the data to the BiodiversityDistribution object
    x$set_predictors(
        bdproto(NULL, PredictorDataset,
              id = new_id(),
              data = env,
              ...
        )
      )
  }
)

#' @name add_predictors
#' @rdname add_predictors
#' @usage \S4method{add_predictors}{BiodiversityDistribution, stars}(x, env)
methods::setMethod(
  "add_predictors",
  methods::signature(x = "BiodiversityDistribution", env = "stars"),
  function(x, env, names = NULL, transform = 'scale', derivates = 'none', bgmask = TRUE, harmonize_na = FALSE,
           explode_factors = FALSE, int_variables = NULL, priors = NULL, ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            !missing(env))
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Taking first time entry from object.')

    # Convert to raster
    env <- stars_to_raster(env, which = 1)
    if(is.list(env)) env <- env[[1]]
    add_predictors(x, env, names, transform, derivates, bgmask, harmonize_na, explode_factors, int_variables, priors, ...)
  }
)

# Add elevational delineation as predictor ----

#' Create lower and upper limits for an elevational range and add them as separate predictors
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param layer A [`character`] stating the elevational layer in the Distribution object or [`RasterLayer`] object.
#' @param lower [`numeric`] value for a lower elevational preference of a species.
#' @param upper [`numeric`] value for a upper elevational preference of a species.
#' @param transform [`character`] Any optional transformation to be applied. Usually not needed (Default: \code{"none"}).
#' @name add_predictor_elevationpref
NULL

#' @name add_predictor_elevationpref
#' @rdname add_predictor_elevationpref
#' @exportMethod add_predictor_elevationpref
#' @export
methods::setGeneric(
  "add_predictor_elevationpref",
  signature = methods::signature("x", "layer", "lower", "upper", "transform"),
  function(x, layer, lower, upper, transform = "none") standardGeneric("add_predictor_elevationpref"))

#' @name add_predictor_elevationpref
#' @rdname add_predictor_elevationpref
#' @usage \S4method{add_predictor_elevationpref}{BiodiversityDistribution, ANY, numeric, numeric, character}(x, layer, lower, upper, transform)
methods::setMethod(
  "add_predictor_elevationpref",
  methods::signature(x = "BiodiversityDistribution", layer = "ANY", lower = "numeric", upper = "numeric"),
  function(x, layer, lower, upper, transform = "none") {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.Raster(layer) || is.character(layer),
                            is.numeric(lower) || is.na(lower),
                            is.numeric(upper) || is.na(upper),
                            is.character(transform)
    )
    # Messager
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Formatting elevational preference predictors...')

    # If layer is a character, check that it is in the provided object
    if(is.character(layer)){
      assertthat::assert_that(layer %in% x$get_predictor_names())
      layer <- x$predictors$get_data()[[layer]]
    } else {
      # If it is a raster
      # Check that background and range align, otherwise raise error
      if(compareRaster(layer, x$background,stopiffalse = FALSE)){
        warning('Supplied range does not align with background! Aligning them now...')
        layer <- alignRasters(layer, x$background, method = 'bilinear', func = mean, cl = FALSE)
      }
    }

    # Format lower and upper preferences
    if(is.na(lower)) lower <- raster::cellStats(layer, "min")
    if(is.na(upper)) upper <- raster::cellStats(layer, "max")

    # Now create thresholded derivatives of lower and upper elevation
    ras1 <- layer
    # ras2[ras2 < lower] <- 0; ras2[ras2 > upper] <- 0; ras2[ras2 > 0] <- 1 # Both ways
    ras1[layer < lower] <- 0; ras1[ras1 > lower] <- 1
    ras2 <- layer
    ras2[ras2 < upper] <- 0; ras2[ras2 > 0] <- 1
    # If the minimum of those layers have equal min and max
    if(raster::cellStats(ras1, "min") == raster::cellStats(ras1, "max")){
      o <- ras2
      # Ensure that all layers have a minimum and a maximum
      o[is.na(o)] <- 0; o <- raster::mask(o, x$background)
      names(o) <- c('elev_high')
    } else {
      o <- raster::stack(ras1, ras2)
      # Ensure that all layers have a minimum and a maximum
      o[is.na(o)] <- 0; o <- raster::mask(o, x$background)
      names(o) <- c('elev_low', 'elev_high')
    }
    rm(ras1,ras2)

    # Add as predictor
    if(is.Waiver(x$predictors)){
      x <- add_predictors(x, env = o, transform = transform, derivates = 'none')
    } else {
      for(n in names(o)){
        r <- o[[n]]
        # If predictor transformation is specified, apply
        if(transform != "none") r <- predictor_transform(r, option = transform)
        x$predictors$set_data(n, r)
      }
    }
    return(x)
  }
)

# Add species ranges as predictor ----
#' Add a range of a species as predictor to a distribution object
#'
#' @description
#' This function allows to add a species range which is usually drawn by experts in a separate process
#' as spatial explicit prior. Both [`sf`] and [`Raster`]-objects are supported as input.
#'
#' Users are advised to look at the [`bossMaps`] R-package presented as part of Merow et al. (2017),
#' which allows flexible calculation of non-linear distance transforms from the boundary of the range.
#' Outputs of this package could be added directly to this function.
#' **Note that this function adds the range as predictor and not as offset. For this purpose a separate function [`add_offset_range()`] exists.**
#'
#' Additional options allow to include the range either as \code{"binary"} or as \code{"distance"} transformed
#' predictor. The difference being that the range is either directly included as presence-only predictor or
#' alternatively with a linear distance transform from the range boundary. The parameter
#' \code{"distance_max"} can be specified to constrain this distance transform.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param layer A [`sf`] or [`Raster`] object with the range for the target feature.
#' @param method [`character`] describing how the range should be included (\code{"binary"} | \code{"distance"}).
#' @param distance_max Numeric threshold on the maximum distance (Default: \code{NULL}).
#' @param priors A [`PriorList-class`] object. Default is set to NULL which uses default prior assumptions
#' @references
#' * Merow, C., Wilson, A. M., & Jetz, W. (2017). Integrating occurrence data and expert maps for improved species range predictions. Global Ecology and Biogeography, 26(2), 243–258. https://doi.org/10.1111/geb.12539
#' @name add_predictor_range
NULL

#' @name add_predictor_range
#' @rdname add_predictor_range
#' @exportMethod add_predictor_range
#' @export
methods::setGeneric(
  "add_predictor_range",
  signature = methods::signature("x", "layer", "method"),
  function(x, layer, method = 'distance', distance_max = NULL, priors = NULL) standardGeneric("add_predictor_range"))

#' Function for when distance raster is directly supplied (precomputed)
#' @name add_predictor_range
#' @rdname add_predictor_range
#' @usage \S4method{add_predictor_range}{BiodiversityDistribution, raster}(x, layer)
methods::setMethod(
  "add_predictor_range",
  methods::signature(x = "BiodiversityDistribution", layer = "RasterLayer"),
  function(x, layer, method = 'precomputed_range', priors = NULL) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.Raster(layer),
                            is.character(method)
    )
    # Messager
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding range predictors...')

    # Check that background and range align, otherwise raise error
    if(compareRaster(layer, x$background,stopiffalse = FALSE)){
      warning('Supplied range does not align with background! Aligning them now...')
      layer <- alignRasters(layer, x$background, method = 'bilinear', func = mean, cl = FALSE)
    }
    names(layer) <- method

    # Add as predictor
    if(is.Waiver(x$predictors)){
      x <- add_predictors(x, env = layer, transform = 'none',derivates = 'none', priors)
    } else {
      x$predictors$set_data('range_distance', layer)
      if(!is.null(priors)) {
        # FIXME: Ideally attempt to match varnames against supplied predictors vis match.arg or similar
        assertthat::assert_that( all( priors$varnames() %in% names(layer) ) )
        x <- x$set_priors(priors)
      }
    }
    return(x)
  }
)

#' @name add_predictor_range
#' @rdname add_predictor_range
#' @usage \S4method{add_predictor_range}{BiodiversityDistribution, sf, vector}(x, layer, method)
methods::setMethod(
  "add_predictor_range",
  methods::signature(x = "BiodiversityDistribution", layer = "sf", method = "character"),
  function(x, layer, method = 'distance', distance_max = Inf, priors = NULL ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.character(method),
                            inherits(layer, 'sf'),
                            method %in% c('binary','distance'),
                            is.null(distance_max) || is.numeric(distance_max) || is.infinite(distance_max),
                            is.null(priors) || inherits(priors,'PriorList')
    )
    # Messager
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding range predictors...')

    # Reproject if necessary
    if(sf::st_crs(layer) != sf::st_crs(x$background)) layer <- sf::st_transform(layer, sf::st_crs(x$background))

    # Template raster for background
    if(!is.Waiver(x$predictors)){
      temp <- emptyraster(x$predictors$get_data())
    } else {
      # TODO: Eventually make this work better
      myLog('[Setup]','red','CAREFUL - This might not work without predictors already in the model.')
      temp <- raster::raster(extent(x$background),resolution = 1)
    }

    # Rasterize the range
    if( 'fasterize' %in% installed.packages()[,1] ){
      ras_range <- fasterize::fasterize(layer, temp, field = NULL)
      if(inherits(ras_range,"try-error")){
        myLog('[Setup]','yellow','Fasterize package needs to be re-installed!')
        ras_range <- raster::rasterize(layer, temp, field = 1, background = NA)
      }
    } else {
      ras_range <- raster::rasterize(layer, temp, field = 1, background = NA)
    }

    # -------------- #
    if(method == 'binary'){
      dis <- ras_range
      dis[is.na(dis)] <- 0
      # Mask with temp again
      dis <- raster::mask(dis, x$background)
      names(dis) <- 'binary_range'
    } else if(method == 'distance'){
      # Calculate the linear distance from the range
      dis <- raster::gridDistance(ras_range, origin = 1)
      dis <- raster::mask(dis, x$background)
      # If max distance is specified
      if(!is.null(distance_max) && !is.infinite(distance_max)){
        dis[dis > distance_max] <- NA # Set values above threshold to NA
        attr(dis, "distance_max") <- distance_max
      }
      # Convert to relative for better scaling in predictions
      dis <- 1 - (dis / cellStats(dis,'max'))
      # Set NA to 0 and mask again
      dis[is.na(dis)] <- 0
      dis <- raster::mask(dis, x$background)
      names(dis) <- 'distance_range'
    }

    # If priors have been set, save them in the distribution object
    if(!is.null(priors)) {
      # FIXME: Ideally attempt to match varnames against supplied predictors vis match.arg or similar
      assertthat::assert_that( all( priors$varnames() %in% names(dis) ) )
      x <- x$set_priors(priors)
    }

    # Add as predictor
    if(is.Waiver(x$predictors)){
      x <- add_predictors(x,env = dis,transform = 'none',derivates = 'none')
    } else {
      x$predictors$set_data('range_distance', dis)
    }
    return(x)
  }
)

#' Remove specific predictors from a [distribution] object
#'
#' @description
#' Remove a particular variable from an [distribution] object with a
#' [`PredictorDataset-class`].
#' See Examples.
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param names [`vector`] A Vector of character names describing the environmental stack.
#' @examples
#' \dontrun{
#' distribution(background) %>%
#'  add_predictors(my_covariates) %>%
#'  rm_predictors(names = "Urban")
#' }
#' @name rm_predictors
NULL

#' @name rm_predictors
#' @rdname rm_predictors
#' @exportMethod rm_predictors
#' @export
methods::setGeneric(
  "rm_predictors",
  signature = methods::signature("x", "names"),
  function(x, names) standardGeneric("rm_predictors"))

#' @name rm_predictors
#' @rdname rm_predictors
#' @usage \S4method{rm_predictors}{BiodiversityDistribution,vector}(x, names)
methods::setMethod(
  "rm_predictors",
  methods::signature(x = "BiodiversityDistribution", names = "character"),
  # rm_predictors ----
  function(x, names ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.character(names) || assertthat::is.scalar(names) || is.vector(names)
                            )
    # TODO: Maybe implement a flexible wildcard, base::startsWith()
    # Is there anything to remove
    assertthat::assert_that(!is.Waiver(x$predictors),
                            all( names %in% x$get_predictor_names() ),
                            msg = 'Suggested variables not in model!')

    # Finally set the data to the BiodiversityDistribution object
    x$rm_predictors(names)
  }
)

#' Select specific predictors from a [distribution] object
#'
#' @description
#' This function allows - out of a [`character`] vector with the names
#' of an already added [`PredictorDataset-class`] object - to select a particular set of predictors.
#' See Examples.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param names [`vector`] A Vector of character names describing the environmental stack.
#' @examples
#' \dontrun{
#' distribution(background) %>%
#'  add_predictors(my_covariates) %>%
#'  sel_predictors(names = c("Forest", "Elevation"))
#' }
#' @name sel_predictors
NULL

#' @name sel_predictors
#' @rdname sel_predictors
#' @exportMethod sel_predictors
#' @export
methods::setGeneric(
  "sel_predictors",
  signature = methods::signature("x", "names"),
  function(x, names) standardGeneric("sel_predictors"))

#' @name sel_predictors
#' @rdname sel_predictors
#' @usage \S4method{sel_predictors}{BiodiversityDistribution,vector}(x, names)
methods::setMethod(
  "sel_predictors",
  methods::signature(x = "BiodiversityDistribution", names = "character"),
  # sel_predictors ----
  function(x, names ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.character(names) || assertthat::is.scalar(names) || is.vector(names)
    )
    # TODO: Maybe implement a flexible wildcard, base::startsWith()
    # Is there anything to remove
    assertthat::assert_that(!is.Waiver(x$predictors),
                            any( names %in% x$get_predictor_names() ),
                            msg = 'Suggested variables not in model!')

    # Get current predictors
    varnames <- x$get_predictor_names()
    varnames <- varnames[which(varnames %notin% names)]

    # Remove all predictors listed
    if(length(varnames)>=1) x$rm_predictors(varnames)
  }
)

# ---------------- #
# Add predictor actions for scenario objects ----
#' @name add_predictors
#' @rdname add_predictors
#' @usage \S4method{add_predictors}{BiodiversityScenario, stars}(x, env)
methods::setMethod(
  "add_predictors",
  methods::signature(x = "BiodiversityScenario", env = "stars"),
  function(x, env, names = NULL, transform = 'none', derivates = 'none', harmonize_na = FALSE, ... ) {
    # Try and match transform and derivatives arguments
    transform <- match.arg(transform, c('none','pca', 'scale', 'norm', 'windsor') , several.ok = TRUE)
    derivates <- match.arg(derivates, c('none','thresh', 'hinge', 'quadratic', 'bin') , several.ok = TRUE)

    assertthat::validate_that(inherits(env,'stars'),msg = 'Projection rasters need to be stars stack!')
    assertthat::assert_that(inherits(x, "BiodiversityScenario"),
                            transform == 'none' || all( transform %in% c('pca', 'scale', 'norm', 'windsor') ),
                            derivates == 'none' || all( derivates %in% c('thresh', 'hinge', 'quadratic', 'bin') ),
                            is.null(names) || assertthat::is.scalar(names) || is.vector(names),
                            is.logical(harmonize_na)
    )
    # Some stars checks
    assertthat::validate_that(length(env) >= 1)

    # Get model object
    obj <- x$get_model()

    # Messenger
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding scenario predictors...')

    # Rename attributes if names is specified
    if(!is.null(names)){
      assertthat::assert_that(length(names) == length(env))
      names(env) <- names
    }

    # Harmonize NA values
    if(harmonize_na){
      stop('Missing data harmonization for stars not yet implemented!') #TODO
      if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Harmonizing missing values...')
      env <- predictor_homogenize_na(env, fill = FALSE)
    }

    # Standardization and scaling
    if('none' %notin% transform){
      if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Transforming predictors...')
      for(tt in transform) env <- predictor_transform(env, option = tt)
    }

    # # Calculate derivates if set
    if('none' %notin% derivates){
      # Get variable names
      varn <- obj$get_coefficients()[['Feature']]
      # Are there any derivates present in the coefficients?
      if(any( length( grep("hinge__|bin__|quad__|thresh__", varn ) ) > 0 )){
          if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Creating predictor derivates...')
          for(dd in derivates){
            if(any(grep(dd, varn))){
              env <- predictor_derivate(env, option = dd, deriv = varn)
            } else {
              if(getOption('ibis.setupmessages')) myLog('[Setup]','red', paste0(derivates,' derivates should be created, but not found among coefficients!'))
            }
          }
      } else {
        if(getOption('ibis.setupmessages')) myLog('[Setup]','red','No derivates found among coefficients. None created for projection!')
      }
    }

    # Get, guess and format Time period
    env_dim <- stars::st_dimensions(env)
    timeperiod <- stars::st_get_dimension_values(env,
                                                 grep("year|time|date", names(env_dim), ignore.case = TRUE, value = TRUE)
                                                 )

    # Check whether predictors already exist, if so overwrite
    # TODO: In the future one could think of supplying predictors of varying grain
    if(!is.Waiver(x$predictors)) myLog('[Setup]','yellow','Overwriting existing predictors.')

    # Finally set the data to the BiodiversityScenario object
    x$set_predictors(
      bdproto(NULL, PredictorDataset,
              id = new_id(),
              data = env,
              timeperiod = timeperiod,
              ...
      )
    )
  }
)

# --------------------- #
#### GLOBIOM specific code ----

#' Add GLOBIOM-DownScaleR derived predictors to a Biodiversity distribution object
#'
#' @description
#' This is a customized function to format and add downscaled land-use shares from
#' the [Global Biosphere Management Model (GLOBIOM)](https://iiasa.github.io/GLOBIOM/) to a
#' [distribution] or [BiodiversityScenario] in ibis.iSDM. GLOBIOM is a partial-equilibrium model
#' developed at IIASA and represents land-use sectors with a rich set of environmental and
#' socio-economic parameters, where for instance the agricultural and forestry sector are estimated through
#' dedicated process-based models. GLOBIOM outputs are spatial explicit and usually at a half-degree resolution globally.
#' For finer grain analyses GLOBIOM outputs can be produced in a downscaled format with a
#' customized statistical [downscaling module](https://github.com/iiasa/DownScale).
#'
#' The purpose of this script is to format the GLOBIOM outputs of *DownScale* for the use in the
#' ibis.iSDM package.
#' @details
#' See [`add_predictors()`] for additional parameters and customizations.
#' For more (manual) control the function for formatting the GLOBIOM data can also be
#' called directly via `formatGLOBIOM()`.
#'
#' @param x A [`BiodiversityDistribution-class`] or [`BiodiversityScenario-class`] object.
#' @param fname A [`character`] pointing to a netCDF with the GLOBIOM data.
#' @param names A [`vector`] of character names describing the environmental stack in case they should be renamed (Default: \code{NULL}).
#' @param transform A [`vector`] stating whether predictors should be preprocessed in any way (Options: \code{'none'},\code{'pca'}, \code{'scale'}, \code{'norm'})
#' @param derivates A Boolean check whether derivate features should be considered (Options: \code{'none'}, \code{'thresh'}, \code{'hinge'}, \code{'quad'}) )
#' @param bgmask Check whether the environmental data should be masked with the background layer (Default: \code{TRUE})
#' @param harmonize_na A [`logical`] value indicating of whether NA values should be harmonized among predictors (Default: \code{FALSE})
#' @param priors A [`PriorList-class`] object. Default is set to \code{NULL} which uses default prior assumptions.
#' @param ... Other parameters passed down
#' @seealso [add_predictors]
#' @examples
#' \dontrun{
#'  obj <- distribution(background) %>%
#'         add_predictors_globiom(fname = "", transform = 'none')
#'  obj
#' }
#' @name add_predictors_globiom
NULL

#' @name add_predictors_globiom
#' @rdname add_predictors_globiom
#' @exportMethod add_predictors_globiom
#' @export
methods::setGeneric(
  "add_predictors_globiom",
  signature = methods::signature("x", "fname"),
  function(x, fname, names = NULL, transform = 'none', derivates = 'none', bgmask = TRUE, harmonize_na = FALSE,
           priors = NULL, ...) standardGeneric("add_predictors_globiom"))

#' @name add_predictors_globiom
#' @rdname add_predictors_globiom
#' @usage \S4method{add_predictors_globiom}{BiodiversityDistribution, character}(x, fname)
methods::setMethod(
  "add_predictors_globiom",
  methods::signature(x = "BiodiversityDistribution", fname = "character"),
  function(x, fname, names = NULL, transform = 'none', derivates = 'none', bgmask = TRUE, harmonize_na = FALSE,
           priors = NULL, ... ) {
    # Try and match transform and derivatives arguments
    transform <- match.arg(transform, c('none','pca', 'scale', 'norm', 'windsor'), several.ok = TRUE)
    derivates <- match.arg(derivates, c('none','thresh', 'hinge', 'quadratic', 'bin'), several.ok = TRUE)

    # Check that file exists and has the correct endings
    assertthat::assert_that(is.character(fname),
                            file.exists(fname),
                            assertthat::is.readable(fname),
                            assertthat::has_extension(fname, "nc"),
                            msg = "The provided path to GLOBIOM land-use shares could not be found or is not readable!"
                            )

    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.null(names) || assertthat::is.scalar(names) || is.vector(names),
                            is.null(priors) || inherits(priors,'PriorList')
    )

    # Messenger
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Formatting GLOBIOM inputs for species distribution modelling.')

    # Get and format the GLOBIOM data
    env <- formatGLOBIOM(fname = fname,
                         oftype = "raster",
                         period = "reference",
                         template = x$background
                         )

    if(is.list(env)) env <- env[[1]] # Take the first reference entry
    assertthat::assert_that(is.Raster(env),
                            raster::nlayers(env)>0)

    if(!is.null(names)) {
      assertthat::assert_that(nlayers(env)==length(names),
                              all(is.character(names)),
                              msg = 'Provided names not of same length as environmental data.')
      # Set names of env
      names(env) <- names
    }

    # Check that all names allowed
    problematic_names <- grep("offset|w|weight|spatial_offset|Intercept|spatial.field", names(env),fixed = TRUE)
    if( length(problematic_names)>0 ){
      stop(paste0("Some predictor names are not allowed as they might interfere with model fitting:", paste0(names(env)[problematic_names],collapse = " | ")))
    }

    # If priors have been set, save them in the distribution object
    if(!is.null(priors)) {
      assertthat::assert_that( all( priors$varnames() %in% names(env) ) )
      x <- x$set_priors(priors)
    }
    # Harmonize NA values
    if(harmonize_na){
      if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Harmonizing missing values...')
      env <- predictor_homogenize_na(env, fill = FALSE)
    }

    # Standardization and scaling
    if('none' %notin% transform){
      if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Transforming predictors...')
      for(tt in transform) env <- predictor_transform(env, option = tt)
    }

    # Calculate derivates if set
    if('none' %notin% derivates){
      if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Creating predictor derivates...')
      new_env <- raster::stack()
      for(dd in derivates) new_env <- raster::addLayer(new_env, predictor_derivate(env, option = dd) )

      # Add to env
      env <- raster::addLayer(env, new_env)
    }

    # Generally not relevant for GLOBIOM unless created as derivate
    attr(env, 'has_factors') <- FALSE

    # Assign an attribute to this object to keep track of it
    attr(env,'transform') <- transform

    # Mask predictors with existing background layer
    if(bgmask){
      env <- raster::mask(env, mask = x$background)
      env <- raster::stack(env)
    }

    # Check whether predictors already exist, if so overwrite
    if(!is.Waiver(x$predictors)) myLog('[Setup]','yellow','Overwriting existing predictors.')

    # Finally set the data to the BiodiversityDistribution object
    x$set_predictors(
      bdproto(NULL, PredictorDataset,
              id = new_id(),
              data = env,
              ...
      )
    )
  }
)

#' @name add_predictors_globiom
#' @rdname add_predictors_globiom
#' @usage \S4method{add_predictors_globiom}{BiodiversityScenario, character}(x, fname)
methods::setMethod(
  "add_predictors_globiom",
  methods::signature(x = "BiodiversityScenario", fname = "character"),
  function(x, fname, names = NULL, transform = 'none', derivates = 'none', harmonize_na = FALSE, ... ) {
    # Try and match transform and derivatives arguments
    transform <- match.arg(transform, c('none','pca', 'scale', 'norm', 'windsor') , several.ok = TRUE)
    derivates <- match.arg(derivates, c('none','thresh', 'hinge', 'quadratic', 'bin') , several.ok = TRUE)

    # Check that file exists and has the correct endings
    assertthat::assert_that(is.character(fname),
                            file.exists(fname),
                            assertthat::is.readable(fname),
                            assertthat::has_extension(fname, "nc"),
                            msg = "The provided path to GLOBIOM land-use shares could not be found or is not readable!"
    )
    assertthat::assert_that(inherits(x, "BiodiversityScenario"),
                            is.null(names) || assertthat::is.scalar(names) || is.vector(names),
                            is.logical(harmonize_na)
    )

    # Get model object
    obj <- x$get_model()

    # Messenger
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding GLOBIOM predictors to scenario object...')

    # Get and format the GLOBIOM data
    env <- formatGLOBIOM(fname = fname,
                         oftype = "stars",
                         period = "projection",
                         template = obj$model$background
    )
    assertthat::assert_that( inherits(env, "stars") )

    # Rename attributes if names is specified
    if(!is.null(names)){
      assertthat::assert_that(length(names) == length(env))
      names(env) <- names
    }

    # Harmonize NA values
    if(harmonize_na){
      stop('Missing data harmonization for stars not yet implemented!') #TODO
      if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Harmonizing missing values...')
      env <- predictor_homogenize_na(env, fill = FALSE)
    }

    # Standardization and scaling
    if('none' %notin% transform){
      if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Transforming predictors...')
      for(tt in transform) env <- predictor_transform(env, option = tt)
    }

    # # Calculate derivates if set
    if('none' %notin% derivates){
      # Get variable names
      varn <- obj$get_coefficients()[['Feature']]
      # Are there any derivates present in the coefficients?
      if(any( length( grep("hinge__|bin__|quad__|thresh__", varn ) ) > 0 )){
        if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Creating predictor derivates...')
        for(dd in derivates){
          if(any(grep(dd, varn))){
            env <- predictor_derivate(env, option = dd, deriv = varn)
          } else {
            if(getOption('ibis.setupmessages')) myLog('[Setup]','red', paste0(derivates,' derivates should be created, but not found among coefficients!'))
          }
        }
      } else {
        if(getOption('ibis.setupmessages')) myLog('[Setup]','red','No derivates found among coefficients. None created for projection!')
      }
    }

    # Get and format Time period
    env_dim <- stars::st_dimensions(env)
    timeperiod <- stars::st_get_dimension_values(env, "time", center = TRUE)
    if(is.numeric(timeperiod)){
      # Format to Posix. Assuming years only
      timeperiod <- as.POSIXct(paste0(timeperiod,"-01-01"))
    }
    if(anyNA(timeperiod)) stop('Third dimension is not a time value!')

    # Check whether predictors already exist, if so overwrite
    if(!is.Waiver(x$predictors)) myLog('[Setup]','yellow','Overwriting existing predictors.')

    # Finally set the data to the BiodiversityScenario object
    x$set_predictors(
      bdproto(NULL, PredictorDataset,
              id = new_id(),
              data = env,
              timeperiod = timeperiod,
              ...
      )
    )
  }
)

#' Function to format a prepared GLOBIOM netCDF file for use in Ibis.iSDM
#'
#' @description
#' This function expects a downscaled GLOBIOM output as created in the BIOCLIMA project.
#' Likely of little use for anyone outside IIASA.
#'
#' @param fname A filename in [`character`] pointing to a GLOBIOM output in netCDF format.
#' @param oftype A [`character`] denoting the output type (Default: \code{'raster'}).
#' @param ignore A [`vector`] of variables to be ignored (Default: \code{NULL}).
#' @param period A [`character`] limiting the period to be returned from the formatted data.
#' Options include \code{"reference"} for the first entry, \code{"projection"} for all entries but the first,
#' and \code{"all"} for all entries (Default: \code{"reference"}).
#' @param template An optional [`RasterLayer`] object towards which projects should be transformed.
#' @param verbose [`logical`] on whether to be chatty.
#'
#' @examples \dontrun{
#' # Expects a filename pointing to a netCDF file.
#' covariates <- formatBIOCLIMA(fname)
#' }
#' @keywords internal, utils
formatGLOBIOM <- function(fname, oftype = "raster", ignore = NULL,
                          period = "all", template = NULL,
                          verbose = getOption("ibis.setupmessages")){
  assertthat::assert_that(
    file.exists(fname),
    assertthat::has_extension(fname, "nc"),
    is.character(oftype),
    is.null(ignore) || is.character(ignore),
    is.character(period),
    is.character(fname),
    is.logical(verbose)
  )
  period <- match.arg(period, c("reference", "projection", "all"), several.ok = FALSE)
  check_package("stars")
  check_package("dplyr")
  check_package("cubelyr")
  check_package("ncdf4")

  # Try and load in the GLOBIOM file to get the attributes
  fatt <- ncdf4::nc_open(fname)
  if(verbose) myLog('[Setup]','green',"Found ", fatt$ndims, " dimensions and ", fatt$nvars, " variables")

  # Get all dimension names and variable names
  dims <- names(fatt$dim)
  vars <- names(fatt$var)
  if(!is.null(ignore)) assertthat::assert_that( all( ignore %in% vars ) )

  attrs <- list() # For storing the attributes
  sc <- vector() # For storing the scenario files

  # Now open the netcdf file with stats
  if( length( grep("netcdf", stars:::detect.driver(fname), ignore.case = TRUE) )>0 ){
    if(verbose){
      myLog('[Predictor]','green',"Loading in predictor file...")
      pb <- progress::progress_bar$new(total = length(vars),
                                       format = "Loading :variable (:spin) [:bar] :percent")
    }

    for(v in vars) {
      if(verbose) pb$tick(tokens = list(variable = v))
      if(!is.null(ignore)) if(ignore == v) next()

      # Get and save the attributes of each variable
      attrs[[v]] <- ncdf4::ncatt_get(fatt, varid = v, verbose = FALSE)

      # Load in the variable
      suppressWarnings(
        suppressMessages(
          ff <- stars::read_ncdf(fname,
                                 var = v,
                                 proxy = FALSE,
                                 make_time = TRUE, # Make time on 'time' band
                                 make_units = FALSE # To avoid unnecessary errors due to unknown units
          )
        )
      )

      # Sometimes variables don't seem to have a time dimension
      if(!"time" %in% names(stars::st_dimensions(ff))) next()

      # Crop to background extent if set
      if(!is.null(template)){
        # FIXME: Currently this code, while working clips too much of Europe.
        # Likely need to
        # bbox <- sf::st_bbox(template) |> sf::st_as_sfc() |>
        #   sf::st_transform(crs = sf::st_crs(ff))
        # suppressMessages(
        # ff <- ff |> stars:::st_crop.stars(bbox)
        # )
      }

      # Record dimensions for later
      full_dis <- stars::st_dimensions(ff)

      # Get dimensions other that x,y and time and split
      # Commonly used column names
      check = c("x","X","lon","longitude", "y", "Y", "lat", "latitude", "time", "Time", "year", "Year")
      chk <- which(!names(stars::st_dimensions(ff)) %in% check)

      if(length(chk)>0){
        for(i in chk){
          col_class <- names(stars::st_dimensions(ff))[i]
          # FIXME: Dirty hack to remove forest zoning
          if(length( grep("zone",col_class,ignore.case = T) )>0) next()

          # And class units as description from over
          class_units <- fatt$dim[[col_class]]$units
          class_units <-  class_units |>
            strsplit(";") |>
            # Remove emptyspace and special symbols
            sapply(function(y)  gsub("[^0-9A-Za-z///' ]", "" , y, ignore.case = TRUE) ) |>
            sapply(function(y)  gsub(" ", "" , y, ignore.case = TRUE) )
          # Convert to vector and make names
          class_units <- paste0(
            v, "__",
            make.names(unlist(class_units)) |> as.vector()
          )

          ff <- ff %>% stars:::split.stars(col_class) %>% setNames(nm = class_units)

          # FIXME: Dirty hack to deal with the forest zone dimension
          # If there are more dimensions than 3, aggregate over them
          if( length(stars::st_dimensions(ff)) >3){
            # Aggregate spatial-temporally
            ff <- stars::st_apply(ff, c("longitude", "latitude", "time"), sum, na.rm = TRUE)
          }
        }
      }

      # Finally aggregate
      if(!is.null(template) && is.Raster(template)){
        # FIXME:
        # MJ 14/11/2022 - The code below is buggy, resulting in odd curvilinear extrapolations for Europe
        # Hacky approach now is to convert to raster, crop, project and then convert back.
        ff <- hack_project_stars(ff, template)
        # Make background
        # bg <- stars::st_as_stars(template)
        #
        # # Get resolution
        # res <- sapply(stars::st_dimensions(bg), "[[", "delta")
        # res[1:2] = abs(res[1:2]) # Assumes the first too entries are the coordinates
        # assertthat::assert_that(!anyNA(res))
        #
        # # And warp by projecting and resampling
        # ff <- ff |> st_transform(crs = sf::st_crs(template)) |>
        #   stars::st_warp(crs = sf::st_crs(bg),
        #                           cellsize = res,
        #                           method = "near") |>
        #   stars:::st_transform.stars(crs = sf::st_crs(template))
        # Overwrite full dimensions
        full_dis <- stars::st_dimensions(ff)
      }
      # Now append to vector
      sc <- c(sc, ff)
      rm(ff)
    }
    invisible(gc())
    assertthat::assert_that(length(names(full_dis))>=3)

    # Format sc object as stars and set dimensions again
    sc <- stars::st_as_stars(sc)
    assertthat::assert_that(length(sc)>0)
    full_dis <- full_dis[c(
      grep("x|longitude",names(full_dis), ignore.case = TRUE,value = TRUE),
      grep("y|latitude",names(full_dis), ignore.case = TRUE,value = TRUE),
      grep("year|time",names(full_dis), ignore.case = TRUE,value = TRUE)
      )] # Order assumed to be correct
    stars:::st_dimensions(sc) <- full_dis # Target dimensions

  } else { stop("Fileformat not recognized!")}

  # Get time dimension (without applying offset) so at the centre
  times <- stars::st_get_dimension_values(sc, "time", center = TRUE)

  # Make checks on length of times and if equal to one, drop. check.
  if(length(times)==1){
    if(period == "projection") stop("Found only a single time slot. Projections not possible.")
    if(verbose) myLog('[Setup]','yellow','Found only a single time point in file. Dropping time dimension.')
    # Drop the time dimension
    sc <- stars:::adrop.stars(sc, drop = which(names(stars::st_dimensions(sc)) == "time") )
  }

  # Formate times unit and convert to posix if not already set
  if(is.numeric(times) && length(times) > 1){
    # Assume year and paste0 as properly POSIX formatted
    times <- as.POSIXct( paste0(times, "-01-01") )
    sc <- stars::st_set_dimensions(sc, "time", times)
  }

  # Depending on the period, slice the input data
  if(period == "reference"){
    # Get the first entry and filter
    if(length(times)>1){
      # In case times got removed
      times_first <- stars::st_get_dimension_values(sc, "time")[1]
      sc <- sc %>% stars:::filter.stars(time == times_first)
      times <- times_first;rm(times_first)
    }
  } else if(period == "projection"){
    # Remove the first time entry instead, only using the last entries
    times_allbutfirst <- stars::st_get_dimension_values(sc, "time")[-1]
    sc <- sc %>% stars:::filter.stars(time %in% times_allbutfirst)
    times <- times_allbutfirst; rm(times_allbutfirst)
  }
  assertthat::assert_that(length(times)>0,
                          length(sc)>=1)

  # Create raster template if set
  if(!is.null(template)){
    # Check that template is a raster, otherwise rasterize for GLOBIOM use
    if(inherits(template, "sf")){
      o <- sc %>% stars:::slice.stars("time" , 1) %>% as("Raster")
      if("fasterize" %in% installed.packages()[,1]){
        template <- fasterize::fasterize(sf = template, raster = o, field = NULL)
      } else {
        template <- raster::rasterize(template, o, field = 1)
      }
      rm(o)
    }
  }

  # Now format outputs depending on type, either returning the raster or the stars object
  if(oftype == "raster"){
    # Output type raster, use function from utils_scenario
    out <- stars_to_raster(sc, which = NULL, template = template)
    return(out)
  } else { return( sc ) }
}
