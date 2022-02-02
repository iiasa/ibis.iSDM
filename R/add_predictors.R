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
#' @param priors A [`PriorList-class`] object. Default is set to \code{NULL} which uses default prior assumptions.
#' @param ... Other parameters passed down
#' @note
#' Not every engine does actually need covariates. For instance it is perfectly legit
#' to fit a model with only occurrence data and a spatial latent effect ([add_latent]).
#' This correspondents to a spatial kernel density estimate.
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
  function(x, env, names = NULL, transform = 'scale', derivates = 'none', bgmask = TRUE, harmonize_na = FALSE, priors = NULL, ...) standardGeneric("add_predictors"))

#' @name add_predictors
#' @rdname add_predictors
#' @usage \S4method{add_predictors}{BiodiversityDistribution,RasterBrick}(x, env)
methods::setMethod(
  "add_predictors",
  methods::signature(x = "BiodiversityDistribution", env = "RasterBrick"),
  function(x, env, names = NULL, transform = 'scale', derivates = 'none', bgmask = TRUE, harmonize_na = FALSE, priors = NULL, ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            !missing(env))
    # Convert env to stack if it is a single layer only
    env = raster::stack(env)
    add_predictors(x, env, names, transform, derivates, bgmask, harmonize_na, priors, ...)
  }
)

#' @name add_predictors
#' @rdname add_predictors
#' @usage \S4method{add_predictors}{BiodiversityDistribution,RasterLayer}(x, env)
methods::setMethod(
  "add_predictors",
  methods::signature(x = "BiodiversityDistribution", env = "RasterLayer"),
  function(x, env, names = NULL, transform = 'scale', derivates = 'none', bgmask = TRUE, harmonize_na = FALSE, priors = NULL, ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            !missing(env))
    # Convert env to stack if it is a single layer only
    env = raster::stack(env)
    add_predictors(x, env, names, transform, derivates, bgmask, harmonize_na, priors, ...)
  }
)

# TODO: Support other objects other than Raster stacks such as data.frames and stars objects
#' @name add_predictors
#' @rdname add_predictors
#' @usage \S4method{add_predictors}{BiodiversityDistribution,RasterStack}(x, env)
methods::setMethod(
  "add_predictors",
  methods::signature(x = "BiodiversityDistribution", env = "RasterStack"),
  function(x, env, names = NULL, transform = 'scale', derivates = 'none', bgmask = TRUE, harmonize_na = FALSE, priors = NULL, ... ) {
    # Try and match transform and derivatives arguments
    transform <- match.arg(transform, c('none','pca', 'scale', 'norm') , several.ok = TRUE)
    derivates <- match.arg(derivates, c('none','thresh', 'hinge', 'quadratic') , several.ok = TRUE)

    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.Raster(env),
                            transform == 'none' || all( transform %in% c('pca', 'scale', 'norm') ),
                            derivates == 'none' || all( derivates %in% c('thresh', 'hinge', 'quadratic') ),
                            is.null(names) || assertthat::is.scalar(names) || is.vector(names),
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

    # If priors have been set, save them in the distribution object
    if(!is.null(priors)) {
      # FIXME: Ideally attempt to match varnames against supplied predictors vis match.arg or similar
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

      # Refactor categorical variables
      if(inherits(env_f,'RasterLayer')){
        env_f <- explode_factorized_raster(env_f)
      } else {
        o <- raster::stack()
        for(layer in env_f){
          o <- raster::addLayer(o, explode_factorized_raster(layer))
        }
        env_f <- o;rm(o)
      }
      # Joing back to full raster stack
      env <- raster::stack(env, env_f);rm(env_f)
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
      env <- addLayer(env, new_env)
    }

    # Assign an attribute to this object to keep track of it
    attr(env,'transform') <- transform

    # Mask predictors with existing background layer
    if(bgmask){
      env <- raster::mask(env, mask = x$background)
      # Reratify, work somehow only on stacks
      if(any(is.factor(env))){
        new_env <- raster::stack(env)
        new_env[[which(is.factor(env))]] <- ratify(env[[which(is.factor(env))]])
        new_env <- env;rm(new_env)
      } else env <- raster::stack(env)
    }

    # Check whether predictors already exist, if so overwrite
    # TODO: In the future one could think of supplying predictors of varying grain
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

#' Add a range of a species as predictor to a distribution object
#'
#' As options allow specifying including the range either as binary or distance
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param range A [`sf`] object with the range for the target feature
#' @param method [`character`] describing how the range should be included (binary | distance)
#' @param distance_max Numeric threshold on the maximum distance (Default: NULL)
#' @param priors A [`PriorList-class`] object. Default is set to NULL which uses default prior assumptions
#' @name add_predictor_range
NULL

#' @name add_predictor_range
#' @rdname add_predictor_range
#' @exportMethod add_predictor_range
#' @export
methods::setGeneric(
  "add_predictor_range",
  signature = methods::signature("x", "range", "method"),
  function(x, range, method = 'distance', distance_max = NULL, priors = NULL) standardGeneric("add_predictor_range"))

#' Function for when distance raster is directly supplied (precomputed)
#' @name add_predictor_range
#' @rdname add_predictor_range
#' @usage \S4method{add_predictor_range}{BiodiversityDistribution, raster}(x, range)
methods::setMethod(
  "add_predictor_range",
  methods::signature(x = "BiodiversityDistribution", range = "RasterLayer"),
  function(x, range, method = 'precomputed_range', priors = NULL) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.Raster(range),
                            is.character(method)
    )
    # Messager
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding range predictors...')

    # Check that background and range align, otherwise raise error
    if(compareRaster(range, x$background,stopiffalse = FALSE)){
      warning('Supplied range does not align with background! Aligning them now...')
      range <- alignRasters(range, x$background, method = 'bilinear', func = mean, cl = FALSE)
    }
    names(range) <- method

    # Add as predictor
    if(is.Waiver(x$predictors)){
      x <- add_predictors(x, env = range,transform = 'none',derivates = 'none', priors)
    } else {
      x$predictors$set_data('range_distance', range)
      if(!is.null(priors)) {
        # FIXME: Ideally attempt to match varnames against supplied predictors vis match.arg or similar
        assertthat::assert_that( all( priors$varnames() %in% names(range) ) )
        x <- x$set_priors(priors)
      }
    }
    return(x)
  }
)

#' @name add_predictor_range
#' @rdname add_predictor_range
#' @usage \S4method{add_predictor_range}{BiodiversityDistribution,sf, vector}(x, range, method)
methods::setMethod(
  "add_predictor_range",
  methods::signature(x = "BiodiversityDistribution", range = "sf", method = "character"),
  function(x, range, method = 'distance', distance_max = NULL, priors = NULL ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.character(method),
                            inherits(range, 'sf'),
                            method %in% c('binary','distance'),
                            is.null(distance_max) || is.numeric(distance_max),
                            is.null(priors) || inherits(priors,'PriorList')
    )
    # Messager
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding range predictors...')

    # Reproject if necessary
    if(sf::st_crs(range) != sf::st_crs(x$background)) range <- sf::st_transform(range, sf::st_crs(x$background))

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
      ras_range <- fasterize::fasterize(range, temp, field = NULL)
    } else {
      ras_range <- raster::rasterize(range, temp,field = NULL)
    }

    # -------------- #
    if(method == 'binary'){
      dis <- ras_range
      # Probability for which the species is not in the expert range map
      dis <- dis + 0.001
      # Transform to log-scale
      dis <- log(dis)
      names(dis) <- 'binary_range'
    } else if(method == 'distance'){
      # TODO: The below can be much more sophisticated.
      # - For instance adding a exponential decay
      # Calculate the linear distance
      dis <- raster::distance(ras_range)
      dis <- raster::mask(dis, x$background)
      # Set areas not intersecting with range to 0
      dis <- raster::mask(dis,
                          x$background[unlist( st_intersects(st_buffer(range,0), x$background) ),]
      )
      # If max distance is specified
      if(!is.null(distance_max)) dis[dis > distance_max] <- NA # Set values above threshold to NA
      # Convert to relative for better scaling
      dis <- 1 - (dis / cellStats(dis,'max'))
      # Probability for which the species is not in the expert range map
      dis <- dis + 0.001
      # Transform to log-scale
      dis <- log(dis)
      names(dis) <- 'range_distance'
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
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param names [`vector`] A Vector of character names describing the environmental stack
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

#' Select specific predictors from a distribution object
#' For instance those previously selected
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param names [`vector`] A Vector of character names describing the environmental stack
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

