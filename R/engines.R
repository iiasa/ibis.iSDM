#' @include bdproto-engine.R
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
        crs = inla.CRS(projargs = proj4string(x$background))
      )
    )
    rm(dat)
  }

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
      calculate = function(self, x, ...) {
        stop('To be done.')
        # return success
        invisible(TRUE)
      },
      set_variable_ub = function(self, index, value) {
        self$data$model$ub[index] <- value
        invisible(TRUE)
      },
      set_variable_lb = function(self, index, value) {
        self$data$model$lb[index] <- value
        invisible(TRUE)
      },
      run = function(self, x) {
        stop('To be done.')
        # access input data and parameters
        model <- self$get_data("model")
        p <- self$get_data("parameters")

        # result <- inla(formulaN,family="poisson",
        #                data=inla.stack.data(join.stack),
        #                control.predictor=list(A=inla.stack.A(join.stack), compute=TRUE),
        #                control.family = list(link = "log"),
        #                E = inla.stack.data(join.stack)$e,
        #                control.compute = list(cpo=TRUE, waic = TRUE, dic = TRUE)
        # )
        # out
      }))
}

# TODO:
# Other engines, e.g. STAN or simple estimators such as Maxent?
