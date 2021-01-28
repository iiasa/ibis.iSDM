#' @include utils.R bdproto.R bdproto-biodiversitydistribution.R bdproto-predictors.R
NULL

#' Add predictors to a Biodiversity distribution object
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param env [`Raster`] A [`RasterStack-class`] or [`RasterLayer-class`] object.
#' @param names [`Vector`] A Vector of character names describing the environmental stack
#' @param prep [`Vector`] A vector stating whether predictors should be preprocessed in any way (Options: 'none','pca', 'scale', 'norm')
#' @param derivates A Boolean check whether derivate features should be considered (Options: 'none', 'elu', 'hinge', 'product') )
#' @param bgmask Check whether the environmental data should be masked with the background layer (Default: TRUE)
#' @param ... Other parameters passed down
#'
#' @details TBD
#' @section Notes:
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
  function(x, env, names = NULL, prep = 'scale', derivates = 'none', bgmask = TRUE, ...) standardGeneric("add_predictors"))

# TODO: Support other objects other than Raster stacks
#' @name add_predictors
#' @rdname add_predictors
#' @usage \S4method{add_predictors}{BiodiversityDistribution,RasterStack}(x, env)
methods::setMethod(
  "add_predictors",
  methods::signature(x = "BiodiversityDistribution", env = "RasterStack"),
  function(x, env, names = NULL, prep = 'scale', derivates = 'none', bgmask = TRUE, ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            inherits(env, 'Raster'),
                            prep %in% c('none','pca', 'scale', 'norm'),
                            derivates %in% c('none', 'elu', 'hinge', 'product'),
                            is.null(names) || assertthat::is.scalar(names) || is.vector(names)
    )
    assertthat::assert_that(is_comparable_raster(x$background,env),
                            msg = 'Supplied environmental data not aligned with background.')
    if(!is.null(names)) {
      assertthat::assert_that(nlayers(env)==length(names),
                                                msg = 'Provided names not of same length as environmental data.')
      # Set names of env
      names(env) <- names
    }

    # Standardization and scaling
    # TODO: Think whether this makes sense for categorical predictors? Shouldn't affect the data I think
    if(prep != 'none'){
      # TODO: Possibly allow this for multiple correction at once?
      env <- adjustPredictors(env,option = prep)
    }
    # Assign an attribute to this object to keep track of it
    attr(env,'prep') <- prep

    if(derivates !='none'){
      # TODO: See maxnet etc code...
      message('Derivative features not yet implemented.')
    }
    # Mask predictors with existing background layer
    if(bgmask){
      env <- raster::mask(env,mask = x$background)
    }

    # Check whether predictors already exist, if so overwrite
    # TODO: In the future one could think of supplying predictors of varying grain
    if(!is.Waiver(x$predictors)) message('Overwriting existing environmental predictors.')

    # Finally set the data to the BiodiversityDistribution object
    x$set_predictors(
        bdproto(NULL, PredictorDataset,
              id = new_id(),
              data = env
        )
      )
    return(x)
  }
)
