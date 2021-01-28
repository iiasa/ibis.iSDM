#' @include bdproto-engine.R utils-inla.R
NULL

#' Use INLA as engine
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param optional_mesh A directly supplied [`INLA`] mesh (Default: NULL)
#' @param max_distance Range of maximum distance between two nodes is between 50 and 5000 meter
#' @param offset Offset for INLA mesh
#' @param ... Other variables
#' @name engine_inla
NULL
#' @import INLA
#' @rdname engine_inla
#' @export
engine_inla <- function(x, optional_mesh = NULL, max_distance = c(10,1000), offset = c(1,1), verbose = TRUE,...) {
  # TODO:
  # Find a better way to pass on parameters...
  # assert that arguments are valid
  assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                          inherits(x$background,'Raster'),
                          inherits(optional_mesh,'inla.mesh') || is.null(optional_mesh),
                          is.vector(max_distance),
                          assertthat::is.flag(verbose),
                          requireNamespace("INLA", quietly = TRUE))

  # Set the projection mesh
  if(inherits(optional_mesh,'inla.mesh')) {
    mesh <- optional_mesh
  } else {
    # Background points
    dat <- raster::rasterToPoints(x$background)[,c('x','y')]

    # Make a boundary from the background
    bounds <- raster::boundaries(background, type = 'outer', asNA = TRUE)

    bdry <- INLA::as.inla.mesh.segment(rasterToPolygons(bounds),join = TRUE)
    bdry$loc <- inla.mesh.map(bdry$loc)

    # Prepare the mesh
    suppressWarnings(
      mesh <- INLA::inla.mesh.2d(
        loc = dat[,c("x", "y")], # Define initial triangulation points
        max.edge = max_distance, # maximum distance between two nodes
        # FIXME: This still does not work correctly
#        boundary = bdry, # boundary
        offset = offset, # Offset of outer boundaries
        crs = inla.CRS(projargs = sp::proj4string(x$background))
      )
    )
    rm(dat)
  }

  # Print a message in case there is already an engine object
  if(!is.Waiver(x$engine)) message('Replacing currently selected engine.')

  # Set engine in distribution object
  x$set_engine(
    bdproto(
      "INLA-Engine",
      Engine,
      name = "<INLA>",
      data = list(
        'mesh' = mesh
      ),
      # parameters = parameters(
      #   numeric_parameter("gap", gap, lower_limit = 0),
      #   integer_parameter("time_limit", time_limit, lower_limit = -1L,
      #                     upper_limit = as.integer(.Machine$integer.max)),
      #   integer_parameter("presolve", presolve, lower_limit = -1L,
      #                     upper_limit = 2L),
      #   integer_parameter("threads", threads, lower_limit = 1L,
      #                     upper_limit = parallel::detectCores(TRUE)),
      #   binary_parameter("first_feasible", first_feasible),
      #   binary_parameter("numeric_focus", numeric_focus),
      #   binary_parameter("verbose", verbose)
      #   ),
      # Spatial latent function
      calc_latent_spatial = function(self, alpha = 2,...){
        # Define Matern SPDE model and save
        self$data$latentspatial <- inla.spde2.matern(self$data$mesh,alpha = alpha,...)
        # Make index for spatial field
        self$data$s.index <- inla.spde.make.index(name = "spatial.field",
                                                  n.spde = self$data$latentspatial$n.spde)
        assertthat::assert_that(
          inherits(self$data$latentspatial,'inla.spde'),
          length(self$data$s.index$spatial.field) == self$data$mesh$n
        )
        invisible()
      },
      # Get latent spatial equation bit
      get_equation_latent_spatial = function(self,spatial_object = 'spde'){
        assertthat::assert_that('latentspatial' %in% names(x$engine$data),
                                msg = 'latentspatial object not computed.')
        return(
          paste0('f(spatial.field, model = ',spatial_object,')')
        )
      },
      # Setup computation function
      setup = function(self, model, ...){
        # TODO: Some assert calls

        # Create Projection
        mat_proj <- inla.mesh.projector(
          mesh = self$get_data('mesh'),
          loc = as.matrix(model$data$poipo_values[,c('x','y')])
        )


      },
      # Main training function
      train = function(self, x) {
        stop('To be done.')
        # access input data and parameters

        model <- self$get_data("model")
        p <- self$get_data("parameters")
      }
      ))
}

# TODO:
# Other engines, e.g. STAN or simple estimators such as Maxent?
