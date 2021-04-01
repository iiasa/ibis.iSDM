#' @include utils.R bdproto.R bdproto-biodiversitydistribution.R bdproto-predictors.R
NULL

#' Add predictors to a Biodiversity distribution object
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param env [`Raster`] A [`RasterStack-class`] or [`RasterLayer-class`] object.
#' @param names [`vector`] A Vector of character names describing the environmental stack
#' @param transform [`vector`] A vector stating whether predictors should be preprocessed in any way (Options: 'none','pca', 'scale', 'norm')
#' @param derivates A Boolean check whether derivate features should be considered (Options: 'none', 'thresh', 'hinge', 'product') )
#' @param bgmask Check whether the environmental data should be masked with the background layer (Default: TRUE)
#' @param priors A [`PriorList-class`] object. Default is set to NULL which uses default priors
#' @param ... Other parameters passed down

#' @details TBD
#' @section Notes:
#' @aliases add_predictors
#' @references
#'
#' @examples
#' \dontrun{
#'  TBD
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
  function(x, env, names = NULL, transform = 'scale', derivates = 'none', bgmask = TRUE, priors = NULL, ...) standardGeneric("add_predictors"))

#' @name add_predictors
#' @rdname add_predictors
#' @usage \S4method{add_predictors}{BiodiversityDistribution,RasterBrick}(x, env)
methods::setMethod(
  "add_predictors",
  methods::signature(x = "BiodiversityDistribution", env = "RasterBrick"),
  function(x, env, ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"))
    # Convert env to stack if it is a single layer only
    env = raster::stack(env)
    add_predictors(x, env, ...)
  }
)

#' @name add_predictors
#' @rdname add_predictors
#' @usage \S4method{add_predictors}{BiodiversityDistribution,RasterLayer}(x, env)
methods::setMethod(
  "add_predictors",
  methods::signature(x = "BiodiversityDistribution", env = "RasterLayer"),
  function(x, env, ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"))
    # Convert env to stack if it is a single layer only
    env = raster::stack(env)
    add_predictors(x, env, ...)
  }
)

# TODO: Support other objects other than Raster stacks such as data.frames and stars objects
#' @name add_predictors
#' @rdname add_predictors
#' @usage \S4method{add_predictors}{BiodiversityDistribution,RasterStack}(x, env)
methods::setMethod(
  "add_predictors",
  methods::signature(x = "BiodiversityDistribution", env = "RasterStack"),
  function(x, env, names = NULL, transform = 'scale', derivates = 'none', bgmask = TRUE, priors = NULL, ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            inherits(env, 'Raster'),
                            transform == 'none' || all( transform %in% c('pca', 'scale', 'norm') ),
                            derivates == 'none' || all( derivates %in% c('thresh', 'hinge', 'quadratic') ),
                            is.null(names) || assertthat::is.scalar(names) || is.vector(names),
                            is.null(priors) || inherits(priors,'PriorList')
    )
    assertthat::assert_that(sf::st_crs(x$background) == sf::st_crs(env@crs),
                            msg = 'Supplied environmental data not aligned with background.')
    if(!is.null(names)) {
      assertthat::assert_that(nlayers(env)==length(names),
                                                msg = 'Provided names not of same length as environmental data.')
      # Set names of env
      names(env) <- names
    }

    # Standardization and scaling
    if('none' %notin% transform){
      for(tt in transform) env <- predictor_transform(env, option = tt)
    }

    # Calculate derivates if set
    if('none' %notin% derivates){
      new_env <- raster::stack()
      for(dd in derivates) new_env <- raster::addLayer(new_env, predictor_derivate(env, option = dd) )

      # Add to env
      env <- addLayer(env, new_env)
    }

    # Assign an attribute to this object to keep track of it
    attr(env,'transform') <- transform

    # Mask predictors with existing background layer
    if(bgmask){
      env <- raster::mask(env,mask = x$background)
    }

    # Check whether predictors already exist, if so overwrite
    # TODO: In the future one could think of supplying predictors of varying grain
    if(!is.Waiver(x$predictors)) message('Overwriting existing environmental predictors.')

    # If priors have been set, save them in distribution object
    if(!is.null(priors)) x$set_priors(priors)

    # Finally set the data to the BiodiversityDistribution object
    x$set_predictors(
        bdproto(NULL, PredictorDataset,
              id = new_id(),
              data = env,
              ...
        )
      )
    return(x)
  }
)

#' Remove specific predictors from a distribution object
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
    return(x)
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
    # Check wether any predictors are not in names and raise warning otherwise
    assertthat::validate_that(
      all( names %in% x$get_predictor_names()),msg = 'Not all predictors were present in distribution object.'
    )
    return(x)
  }
)
