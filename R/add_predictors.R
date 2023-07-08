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
#' the used [`Engine`].
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
#' @param env A [`SpatRaster`] or [`stars`] object.
#' @param names A [`vector`] of character names describing the environmental stack in case they should be renamed.
#' @param transform A [`vector`] stating whether predictors should be preprocessed in any way (Options: \code{'none'},\code{'pca'}, \code{'scale'}, \code{'norm'})
#' @param derivates A Boolean check whether derivate features should be considered (Options: \code{'none'}, \code{'thresh'}, \code{'hinge'}, \code{'quad'}) )
#' @param derivate_knots A single [`numeric`] or [`vector`] giving the number of knots for derivate creation if relevant (Default: \code{4}).
#' @param int_variables A [`vector`] with length greater or equal than \code{2} specifying the covariates (Default: \code{NULL}).
#' @param bgmask Check whether the environmental data should be masked with the background layer (Default: \code{TRUE}).
#' @param harmonize_na A [`logical`] value indicating of whether NA values should be harmonized among predictors (Default: \code{FALSE}).
#' @param explode_factors [`logical`] of whether any factor variables should be split up into binary variables (one per class). (Default: \code{FALSE}).
#' @param priors A [`PriorList-class`] object. Default is set to \code{NULL} which uses default prior assumptions.
#' @param ... Other parameters passed down
#' @note
#' **Important:**
#' Not every [`Engine`] supported by the \pkg{ibis.iSDM} R-package allows missing data points
#' among extracted covariates. Thus any observation with missing data is generally removed prior
#' from model fitting. Thus ensure that covariates have appropriate no-data settings (for instance setting \code{NA}
#' values to \code{0} or another out of range constant).
#'
#' Not every engine does actually need covariates. For instance it is perfectly legit
#' to fit a model with only occurrence data and a spatial latent effect ([add_latent_spatial]).
#' This correspondents to a spatial kernel density estimate.
#'
#' Certain names such \code{"offset"} are forbidden as predictor variable names. The function
#' will return an error message if these are used.
#' @aliases add_predictors
#' @examples
#' \dontrun{
#'  obj <- distribution(background) |>
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
  function(x, env, names = NULL, transform = 'scale', derivates = 'none', derivate_knots = 4, int_variables = NULL, bgmask = TRUE,
           harmonize_na = FALSE, explode_factors = FALSE, priors = NULL, ...) standardGeneric("add_predictors"))

#' @name add_predictors
#' @rdname add_predictors
#' @usage \S4method{add_predictors}{BiodiversityDistribution,SpatRasterCollection, ANY, character, character, numeric, ANY, logical, logical, logical, ANY}(x, env, names, transform, derivates, derivate_knots, int_variables, bgmask, harmonize_na, explode_factors, priors, ...)
methods::setMethod(
  "add_predictors",
  methods::signature(x = "BiodiversityDistribution", env = "SpatRasterCollection"),
  function(x, env, names = NULL, transform = 'scale', derivates = 'none', derivate_knots = 4, int_variables = NULL,
           bgmask = TRUE, harmonize_na = FALSE, explode_factors = FALSE, priors = NULL, ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            !missing(env))
    # Convert env to stack if it is a single layer only
    if(!is.Raster(env)) env <- terra::rast(env)
    add_predictors(x, env, names, transform, derivates, derivate_knots, int_variables, bgmask, harmonize_na, explode_factors, priors, ...)
  }
)

#' @name add_predictors
#' @rdname add_predictors
#' @usage \S4method{add_predictors}{BiodiversityDistribution, SpatRaster}(x, env)
methods::setMethod(
  "add_predictors",
  methods::signature(x = "BiodiversityDistribution", env = "SpatRaster"),
  function(x, env, names = NULL, transform = 'scale', derivates = 'none', derivate_knots = 4, int_variables = NULL,
           bgmask = TRUE, harmonize_na = FALSE, explode_factors = FALSE, priors = NULL, ... ) {
    # Try and match transform and derivatives arguments
    transform <- match.arg(transform, c('none','pca', 'scale', 'norm', 'windsor') , several.ok = TRUE)
    derivates <- match.arg(derivates, c('none','thresh', 'hinge', 'quadratic', 'bin', 'interaction') , several.ok = TRUE)

    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.Raster(env),
                            all(transform == 'none') || all( transform %in% c('pca', 'scale', 'norm', 'windsor') ),
                            all(derivates == 'none') || all( derivates %in% c('thresh', 'hinge', 'quadratic', 'bin', 'interaction') ),
                            is.vector(derivate_knots) || is.numeric(derivate_knots),
                            is.null(names) || assertthat::is.scalar(names) || is.vector(names),
                            is.logical(explode_factors),
                            is.null(priors) || inherits(priors,'PriorList')
    )
    # Messenger
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding predictors...')

    if(!is.null(names)) {
      assertthat::assert_that(terra::nlyr(env) == length(names),
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
      env_f <- terra::subset(env, which(is.factor(env)))
      env <- terra::subset(env, which(!is.factor(env)))
      if(explode_factors){
        # Refactor categorical variables
        if(inherits(env_f,'SpatRaster')){
          env_f <- explode_factorized_raster(env_f)
          env <- c(env, env_f)
        } else {
          o <- terra::rast()
          for(layer in names(env_f)){
            suppressWarnings(
              o <-c(o, explode_factorized_raster(env_f[[layer]]))
            )
          }
          env_f <- o;rm(o)
          # Joining back to full raster stack
          env <- c(env, env_f);rm(env_f)
        }
        has_factors <- FALSE # Set to false since factors have been exploded.
      } else { has_factors <- TRUE }
    } else { has_factors <- FALSE }
    assertthat::assert_that(is.Raster(env),
                            is.logical(has_factors),
                            terra::nlyr(env)>=1)

    # Standardization and scaling
    if('none' %notin% transform){
      if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Transforming predictors...')
      for(tt in transform) env <- predictor_transform(env, option = tt)
    }
    assertthat::assert_that(is.Raster(env), terra::nlyr(env)>=1)

    # Calculate derivates if set
    if('none' %notin% derivates){
      if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Creating predictor derivates...')
      # Specific condition for interaction
      if(any(derivates == "interaction")){
        assertthat::assert_that(is.vector(int_variables), length(int_variables)>=2)
        attr(env, "int_variables") <- int_variables
      }
      new_env <- terra::rast()
      for(dd in derivates){
        suppressWarnings(
          new_env <- c(new_env, predictor_derivate(env, option = dd, nknots = derivate_knots, int_variables = int_variables) )
        )
      }

      # Add to env
      env <- c(env, new_env)
    }

    # This is to avoid that they are transformed or similar
    if(has_factors) env <- c(env, env_f)
    attr(env, 'has_factors') <- has_factors

    # Assign an attribute to this object to keep track of it
    attr(env,'transform') <- transform

    # Mask predictors with existing background layer
    if(bgmask){
      env <- terra::mask(env, mask = x$background)
      # Reratify, work somehow only on stacks
      if(has_factors && any(is.factor(env)) ){
        new_env <- env
        new_env[[which(is.factor(env))]] <- as.factor(env[[which(is.factor(env))]])
        env <- new_env;rm(new_env)
      }
    }

    # Check whether predictors already exist, if so overwrite
    if(!is.Waiver(x$predictors)) myLog('[Setup]','yellow','Overwriting existing predictors.')

    # Sanitize names if specified
    if(getOption('ibis.cleannames')) names(env) <- sanitize_names(names(env))

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
  function(x, env, names = NULL, transform = 'scale', derivates = 'none', derivate_knots = 4, int_variables = NULL,
           bgmask = TRUE, harmonize_na = FALSE, explode_factors = FALSE, priors = NULL, ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            !missing(env))
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Taking first time entry from object.')

    # Convert to raster
    env <- stars_to_raster(env, which = 1)
    if(is.list(env)) env <- env[[1]]
    x <- add_predictors(x, env, names, transform, derivates, derivate_knots, int_variables, bgmask, harmonize_na, explode_factors, priors, ...)
    return( x )
  }
)

# Add elevational delineation as predictor ----

#' Create lower and upper limits for an elevational range and add them as separate predictors
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param layer A [`character`] stating the elevational layer in the Distribution object or [`SpatRaster`] object.
#' @param lower [`numeric`] value for a lower elevational preference of a species.
#' @param upper [`numeric`] value for a upper elevational preference of a species.
#' @param transform [`character`] Any optional transformation to be applied. Usually not needed (Default: \code{"none"}).
#' @aliases add_predictor_elevationpref
#' @examples
#' \dontrun{
#' distribution(background) |>
#'   add_predictor_elevationpref(elevation, lower = 200, upper = 1000)
#' }
#'
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
      if(is.Raster(x$background)){
        # Check that background and range align, otherwise raise error
        if(is_comparable_raster(layer, x$background)){
          warning('Supplied range does not align with background! Aligning them now...')
          layer <- alignRasters(layer, x$background, method = 'bilinear', func = mean, cl = FALSE)
        }
      }
    }

    # Format lower and upper preferences
    if(is.na(lower)) lower <- terra::global(layer, "min", na.rm = TRUE)
    if(is.na(upper)) upper <- terra::global(layer, "max", na.rm = TRUE)

    # Now create thresholded derivatives of lower and upper elevation
    ras1 <- layer
    # ras2[ras2 < lower] <- 0; ras2[ras2 > upper] <- 0; ras2[ras2 > 0] <- 1 # Both ways
    ras1[layer < lower] <- 0; ras1[ras1 > lower] <- 1
    ras2 <- layer
    ras2[ras2 < upper] <- 0; ras2[ras2 > 0] <- 1
    # If the minimum of those layers have equal min and max
    if(terra::global(ras1, "min", na.rm = TRUE) == terra::global(ras1, "max", na.rm = TRUE)){
      o <- ras2
      # Ensure that all layers have a minimum and a maximum
      o[is.na(o)] <- 0; o <- terra::mask(o, x$background)
      names(o) <- c('elev_high')
    } else {
      o <- c(ras1, ras2)
      # Ensure that all layers have a minimum and a maximum
      o[is.na(o)] <- 0; o <- terra::mask(o, x$background)
      names(o) <- c('elev_low', 'elev_high')
    }
    rm(ras1,ras2)

    # Sanitize names if specified
    if(getOption('ibis.cleannames')) names(o) <- sanitize_names(names(o))

    # Add as predictor
    if(is.Waiver(x$predictors)){
      x <- add_predictors(x, env = o, transform = transform, derivates = 'none')
    } else {
      for(n in names(o)){
        r <- o[[n]]
        # If predictor transformation is specified, apply
        if(transform != "none") r <- predictor_transform(r, option = transform)
        x$predictors <- x$predictors$set_data(n, r)
        rm(r)
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
#' as spatial explicit prior. Both [`sf`] and [`SpatRaster`]-objects are supported as input.
#'
#' Users are advised to look at the \code{"bossMaps"} R-package presented as part of Merow et al. (2017),
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
#' @param layer A [`sf`] or [`SpatRaster`] object with the range for the target feature.
#' @param method [`character`] describing how the range should be included (\code{"binary"} | \code{"distance"}).
#' @param distance_max Numeric threshold on the maximum distance (Default: \code{NULL}).
#' @param fraction An optional [`SpatRaster`] object that is multiplied with digitized raster layer.
#' Can be used to for example to remove or reduce the expected value (Default: \code{NULL}).
#' @param priors A [`PriorList-class`] object. Default is set to NULL which uses default prior assumptions.
#' @aliases add_predictor_range
#' @examples
#' \dontrun{
#' distribution(background) |>
#'   add_predictor_range(range, method = "distance", distance_max = 2)
#' }
#'
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
  function(x, layer, method = 'distance', distance_max = NULL, fraction = NULL, priors = NULL) standardGeneric("add_predictor_range"))

#' Function for when distance raster is directly supplied (precomputed)
#' @name add_predictor_range
#' @rdname add_predictor_range
#' @usage \S4method{add_predictor_range}{BiodiversityDistribution, SpatRaster, character, ANY, ANY}(x, layer, method, fraction, priors)
methods::setMethod(
  "add_predictor_range",
  methods::signature(x = "BiodiversityDistribution", layer = "SpatRaster"),
  function(x, layer, method = 'precomputed_range', fraction = NULL, priors = NULL) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.Raster(layer),
                            is.Raster(fraction) || is.null(fraction),
                            is.character(method)
    )
    # Messenger
    if(getOption('ibis.setupmessages')) myLog('[Setup]', 'green', 'Adding range predictors...')

    # Check that background and range align, otherwise raise error
    if(is.Raster(layer)){
      if(is_comparable_raster(layer, x$background)){
        warning('Supplied range does not align with background! Aligning them now...')
        layer <- alignRasters(layer, x$background, method = 'bilinear', func = mean, cl = FALSE)
      }
    }
    names(layer) <- method

    # Multiply with fraction layer if set
    if(!is.null(fraction)){
      # Rescale if necessary and set 0 to a small constant 1e-6
      if(terra::global(fraction, "min") < 0) fraction <- predictor_transform(fraction, option = "norm")
      fraction[fraction==0] <- 1e-6
      layer <- layer * fraction
    }

    # Sanitize names if specified
    if(getOption('ibis.cleannames')) names(layer) <- sanitize_names(names(layer))

    # Add as predictor
    if(is.Waiver(x$predictors)){
      x <- add_predictors(x, env = layer, transform = 'none',derivates = 'none', priors)
    } else {
      x$predictors <- x$predictors$set_data('range_distance', layer)
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
#' @usage \S4method{add_predictor_range}{BiodiversityDistribution, sf}(x, layer)
methods::setMethod(
  "add_predictor_range",
  methods::signature(x = "BiodiversityDistribution", layer = "sf"),
  function(x, layer, method = 'distance', distance_max = Inf, fraction = NULL, priors = NULL ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.character(method),
                            inherits(layer, 'sf'),
                            method %in% c('binary','distance'),
                            is.null(fraction) || is.Raster(fraction),
                            is.null(distance_max) || is.numeric(distance_max) || is.infinite(distance_max),
                            is.null(priors) || inherits(priors,'PriorList')
    )
    # Messenger
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding range predictors...')

    # Reproject if necessary
    if(sf::st_crs(layer) != sf::st_crs(x$background)) layer <- sf::st_transform(layer, sf::st_crs(x$background))

    # Template raster for background
    if(!is.Waiver(x$predictors)){
      temp <- emptyraster(x$predictors$get_data())
    } else {
      # TODO: Eventually make this work better
      myLog('[Setup]','red','CAREFUL - This might not work without predictors already in the model.')
      temp <- terra::rast( terra::ext(x$background), resolution = 1)
    }

    # Rasterize the range
    # if( 'fasterize' %in% utils::installed.packages()[,1] ){
    #   ras_range <- try({ fasterize::fasterize(layer, temp, field = NULL) }, silent = TRUE)
    #   if(inherits(ras_range,"try-error")){
    #     myLog('[Setup]','yellow','Fasterize package needs to be re-installed!')
    #     ras_range <- raster::rasterize(layer, temp, field = 1, background = NA)
    #   }
    # } else {
    ras_range <- terra::rasterize(layer, temp, field = 1, background = 0)
    # }

    # -------------- #
    if(method == 'binary'){
      dis <- ras_range
      dis[is.na(dis)] <- 0
      # Mask with temp again
      dis <- terra::mask(dis, x$background)
      names(dis) <- 'binary_range'
    } else if(method == 'distance'){
      # Calculate the linear distance from the range
      dis <- terra::gridDist(ras_range, target = 1)
      dis <- terra::mask(dis, x$background)
      # If max distance is specified
      if(!is.null(distance_max) && !is.infinite(distance_max)){
        dis[dis > distance_max] <- NA # Set values above threshold to NA
        attr(dis, "distance_max") <- distance_max
      } else { distance_max <- terra::global(dis, "max", na.rm = TRUE)[,1] }
      # Grow baseline raster by using an exponentially weighted kernel
      alpha <- 1 / (distance_max / 4 ) # Divide by 4 for a quarter in each direction
      # Grow baseline raster by using an exponentially weighted kernel
      dis <- terra::app(dis, fun = function(x) exp(-alpha * x))
      # Convert to relative for better scaling in predictions
      dis <- (dis / terra::global(dis, 'max', na.rm = TRUE)[,1])

      # Set NA to 0 and mask again
      dis[is.na(dis)] <- 0
      dis <- terra::mask(dis, x$background)
      names(dis) <- 'distance_range'
    }

    # Multiply with fraction layer if set
    if(!is.null(fraction)){
      # Rescale if necessary and set 0 to a small constant 1e-6
      if( terra::global(fraction, "min", na.rm = TRUE)[,1] < 0) fraction <- predictor_transform(fraction, option = "norm")
      fraction[fraction==0] <- 1e-6
      layer <- layer * fraction
    }

    # If priors have been set, save them in the distribution object
    if(!is.null(priors)) {
      # FIXME: Ideally attempt to match varnames against supplied predictors vis match.arg or similar
      assertthat::assert_that( all( priors$varnames() %in% names(dis) ) )
      x <- x$set_priors(priors)
    }

    # Sanitize names if specified
    if(getOption('ibis.cleannames')) names(dis) <- sanitize_names(names(dis))

    # Add as predictor
    if(is.Waiver(x$predictors)){
      x <- add_predictors(x, env = dis, transform = 'none',derivates = 'none')
    } else {
      x$predictors <- x$predictors$set_data('range_distance', dis)
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
#' @aliases rm_predictors
#' @examples
#' \dontrun{
#' distribution(background) |>
#'  add_predictors(my_covariates) |>
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
#' @usage \S4method{rm_predictors}{BiodiversityDistribution, ANY}(x, names)
methods::setMethod(
  "rm_predictors",
  methods::signature(x = "BiodiversityDistribution", names = "character"),
  # rm_predictors ----
  function(x, names) {
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
#' @aliases sel_predictors
#' @examples
#' \dontrun{
#' distribution(background) |>
#'  add_predictors(my_covariates) |>
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
#' @usage \S4method{sel_predictors}{BiodiversityDistribution, ANY}(x, names)
methods::setMethod(
  "sel_predictors",
  methods::signature(x = "BiodiversityDistribution", names = "character"),
  # sel_predictors ----
  function(x, names) {
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
#' @usage \S4method{add_predictors}{BiodiversityScenario, SpatRaster, ANY, character, character, numeric, ANY, logical}(x, env, names, transform, derivates, derivate_knots, int_variables, harmonize_na, ...)
methods::setMethod(
  "add_predictors",
  methods::signature(x = "BiodiversityScenario", env = "SpatRaster"),
  function(x, env, names = NULL, transform = 'none', derivates = 'none',
           derivate_knots = 4, int_variables = NULL, harmonize_na = FALSE, ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityScenario"),
                            !missing(env))
    env <- raster_to_stars(env) # Convert to stars

    add_predictors(x, env, names = names, transform = transform, derivates = derivates,
                   derivate_knots = derivate_knots, int_variables = int_variables, harmonize_na = harmonize_na, ...)
  }
)

#' @name add_predictors
#' @rdname add_predictors
#' @usage \S4method{add_predictors}{BiodiversityScenario, stars, ANY, character, character, numeric, ANY, logical}(x, env,names, transform, derivates, derivate_knots, int_variables, harmonize_na, ...)
methods::setMethod(
  "add_predictors",
  methods::signature(x = "BiodiversityScenario", env = "stars"),
  function(x, env, names = NULL, transform = 'none', derivates = 'none',
           derivate_knots = 4, int_variables = NULL, harmonize_na = FALSE, ... ) {
    # Try and match transform and derivatives arguments
    transform <- match.arg(transform, c('none','pca', 'scale', 'norm', 'windsor', 'percentile'), several.ok = TRUE)
    derivates <- match.arg(derivates, c('none','thresh', 'hinge', 'quadratic', 'bin'), several.ok = TRUE)

    assertthat::validate_that(inherits(env,'stars'), msg = 'Projection rasters need to be stars stack!')
    assertthat::assert_that(inherits(x, "BiodiversityScenario"),
                            is.vector(derivate_knots) || is.numeric(derivate_knots),
                            is.null(int_variables) || is.character(int_variables),
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
              env <- predictor_derivate(env, option = dd, nknots = derivate_knots, deriv = varn, int_variables = int_variables)
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

    # Sanitize names if specified
    if(getOption('ibis.cleannames')) names(env) <- sanitize_names(names(env))

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
