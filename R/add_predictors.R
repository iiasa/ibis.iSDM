#' @include utils.R bdproto.R bdproto-biodiversitydistribution.R bdproto-predictors.R bdproto-biodiversityscenario.R
NULL

#' Add predictors to a Biodiversity distribution object
#'
#' @description
#' This function allows to add predictors to [distribution] or [BiodiversityScenario]
#' objects.
#' @details
#' Predictors are covariates that in spatial projection have to match
#' the geographic projection of the background layer in the [distribution] object.
#' This function furthermore allows to transform or create derivates of provided
#' predictors.
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
#' * \code{'thresh'} - Add threshold transformed predictors.
#' * \code{'hinge'} - Add hinge transformed predictors.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param env A [`RasterStack-class`], [`RasterLayer-class`] or [`stars`] object.
#' @param names A [`vector`] of character names describing the environmental stack in case they should be renamed.
#' @param transform A [`vector`] stating whether predictors should be preprocessed in any way (Options: \code{'none'},\code{'pca'}, \code{'scale'}, \code{'norm'})
#' @param derivates A Boolean check whether derivate features should be considered (Options: \code{'none'}, \code{'thresh'}, \code{'hinge'}, \code{'quad'}) )
#' @param bgmask Check whether the environmental data should be masked with the background layer (Default: \code{TRUE})
#' @param harmonize_na A [`logical`] value indicating of whether NA values should be harmonized among predictors (Default: \code{FALSE})
#' @param explode_factors [`logical`] of whether any factor variables should be split up into binary variables (one per class). (Default: \code{FALSE}).
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
  function(x, env, names = NULL, transform = 'scale', derivates = 'none', bgmask = TRUE, harmonize_na = FALSE, explode_factors = FALSE, priors = NULL, ...) standardGeneric("add_predictors"))

#' @name add_predictors
#' @rdname add_predictors
#' @usage \S4method{add_predictors}{BiodiversityDistribution,RasterBrick}(x, env)
methods::setMethod(
  "add_predictors",
  methods::signature(x = "BiodiversityDistribution", env = "RasterBrick"),
  function(x, env, names = NULL, transform = 'scale', derivates = 'none', bgmask = TRUE, harmonize_na = FALSE, explode_factors = FALSE, priors = NULL, ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            !missing(env))
    # Convert env to stack if it is a single layer only
    env = raster::stack(env)
    add_predictors(x, env, names, transform, derivates, bgmask, harmonize_na, explode_factors = FALSE, priors, ...)
  }
)

#' @name add_predictors
#' @rdname add_predictors
#' @usage \S4method{add_predictors}{BiodiversityDistribution,RasterLayer}(x, env)
methods::setMethod(
  "add_predictors",
  methods::signature(x = "BiodiversityDistribution", env = "RasterLayer"),
  function(x, env, names = NULL, transform = 'scale', derivates = 'none', bgmask = TRUE, harmonize_na = FALSE, explode_factors = FALSE, priors = NULL, ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            !missing(env))
    # Convert env to stack if it is a single layer only
    env = raster::stack(env)
    add_predictors(x, env, names, transform, derivates, bgmask, harmonize_na, explode_factors = FALSE, priors, ...)
  }
)

# TODO: Support other objects other than Raster stacks such as data.frames and stars objects
#' @name add_predictors
#' @rdname add_predictors
#' @usage \S4method{add_predictors}{BiodiversityDistribution,RasterStack}(x, env)
methods::setMethod(
  "add_predictors",
  methods::signature(x = "BiodiversityDistribution", env = "RasterStack"),
  function(x, env, names = NULL, transform = 'scale', derivates = 'none', bgmask = TRUE, harmonize_na = FALSE, explode_factors = FALSE, priors = NULL, ... ) {
    # Try and match transform and derivatives arguments
    transform <- match.arg(transform, c('none','pca', 'scale', 'norm', 'windsor') , several.ok = TRUE)
    derivates <- match.arg(derivates, c('none','thresh', 'hinge', 'quadratic') , several.ok = TRUE)

    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.Raster(env),
                            all(transform == 'none') || all( transform %in% c('pca', 'scale', 'norm') ),
                            all(derivates == 'none') || all( derivates %in% c('thresh', 'hinge', 'quadratic') ),
                            is.null(names) || assertthat::is.scalar(names) || is.vector(names),
                            is.logical(explode_factors),
                            is.null(priors) || inherits(priors,'PriorList')
    )
    assertthat::assert_that(sf::st_crs(x$background) == sf::st_crs(env@crs),
                            msg = 'Supplied environmental data not aligned with background.')
    # Messager
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
      env_f <- raster::subset(env,which(is.factor(env)))
      env <- raster::subset(env, which(!is.factor(env)))
      if(explode_factors){
        # Refactor categorical variables
        if(inherits(env_f,'RasterLayer')){
          env_f <- explode_factorized_raster(env_f)
        } else {
          o <- raster::stack()
          for(layer in env_f){
            o <- raster::addLayer(o, explode_factorized_raster(layer))
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
    o <- raster::stack(ras1, ras2)
    names(o) <- c('elev_low', 'elev_high')
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
    } else {
      ras_range <- raster::rasterize(layer, temp,field = NULL)
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

# Add predictor actions for scenario objects ----
#' @name add_predictors
#' @rdname add_predictors
#' @usage \S4method{add_predictors}{BiodiversityScenario, stars}(x, env)
methods::setMethod(
  "add_predictors",
  methods::signature(x = "BiodiversityScenario", env = "stars"),
  function(x, env, names = NULL, transform = 'none', derivates = 'none', harmonize_na = FALSE, ... ) {
    # Try and match transform and derivatives arguments
    transform <- match.arg(transform, c('none','pca', 'scale', 'norm') , several.ok = TRUE)
    derivates <- match.arg(derivates, c('none','thresh', 'hinge', 'quadratic') , several.ok = TRUE)

    assertthat::validate_that(inherits(env,'stars'),msg = 'Projection rasters need to be stars stack!')
    assertthat::assert_that(inherits(x, "BiodiversityScenario"),
                            transform == 'none' || all( transform %in% c('pca', 'scale', 'norm') ),
                            derivates == 'none' || all( derivates %in% c('thresh', 'hinge', 'quadratic') ),
                            is.null(names) || assertthat::is.scalar(names) || is.vector(names),
                            is.logical(harmonize_na)
    )
    # Some stars checks
    assertthat::validate_that(length(env) >= 1)

    # Messager
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding scenario predictors...')

    # Rename attributes if names is specified
    if(!is.null(names)){
      assertthat::assert_that(length(names) == length(env))
      names(env) <- names
    }

    # Harmonize NA values
    if(harmonize_na){
      stop('Harmonization for stars not yet implemented!') #TODO
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
      stop('Derivate creation for stars not yet implemented!') #TODO
      if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Creating predictor derivates...')
      new_env <- raster::stack()
      for(dd in derivates) new_env <- raster::addLayer(new_env, predictor_derivate(env, option = dd) )

      # Add to env
      env <- addLayer(env, new_env)
    }

    # Get and format Time period
    env_dim <- stars::st_dimensions(env)
    timeperiod <- as.POSIXct(env_dim[[3]]$values$start) # Assumes the third dimension is time
    if(anyNA(timeperiod)) stop('Third dimension is not a time value!')

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

