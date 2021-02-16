#' Calculate area of each voronoi polygon in a INLA mesh
#'
#' @param mesh [`inla.mesh`] mesh object
#' @param region.poly A supplied [`region.poly`] object
#' @param variant A character to which type of area calculation (Default: 'gpc')
#' @returns A [`vector`] with the area of each polygon
#' @import deldir
#' @noRd

mesh_area = function(mesh, region.poly = NULL, variant = 'gpc'){
  assertthat::assert_that(inherits(mesh,'inla.mesh'),
                          is.null(region.poly) || inherits(region.poly,'Spatial'),
                          is.character(variant)
                          )
  # Precalculate the area of each
  # Get areas for Voronoi tiles around each integration point
  dd <- deldir::deldir(mesh$loc[,1], mesh$loc[,2])
  tiles <- deldir::tile.list(dd)

  if(variant == 'gpc'){
    # Try to convert to spatial already
    if(!inherits(region.poly, 'Spatial')) region.poly <- as(region.poly,'Spatial')

    poly.gpc <- as(region.poly@polygons[[1]]@Polygons[[1]]@coords,'gpc.poly')
    w <- sapply(tiles, function(p) rgeos::area.poly(rgeos::intersect(as(cbind(p$x, p$y), 'gpc.poly'), poly.gpc)))
  } else {
    # Convert to Spatial Polygons
    polys <- sp::SpatialPolygons(lapply(1:length(tiles), function(i) {
       p <- cbind(tiles[[i]]$x, tiles[[i]]$y)
       n <- nrow(p)
       sp::Polygons(list(sp::Polygon(p[c(1:n, 1), ])), i)
    }),proj4string = mesh$crs)

    # Calculate area of each polygon in km2
    w <- st_area(
       st_as_sf(polys)
    ) %>% units::set_units("km²") %>% as.numeric()
    w <- w / sum(w,na.rm = TRUE) # Relative area
  }

  assertthat::assert_that(is.vector(w))
  return(w)
}

#' Mesh to polygon script
#' @param mesh [`inla.mesh`] mesh object
#' @returns A [`sf`] object
#' @noRd

mesh_as_sf <- function(mesh) {
  assertthat::assert_that(inherits(mesh,'inla.mesh'),
                          mesh$manifold == 'R2' # Two-dimensional mesh
  )
  # Get the delaunay triangle points
  tv <- mesh$graph$tv

  # Now convert the delaunay triangles to spatial polygons
  dp <- lapply(X = 1:nrow(tv),
               FUN = function(index, points, pointindex) {
    # Retrieve the vertex coordinates of the current triangle
    cur <- pointindex[index, ]
    # Construct a Polygons object to contain the triangle
    Polygons(list(
             Polygon( points[c(cur, cur[1]), ], hole = FALSE)),
             ID = index
             )
  }, points = mesh$loc[, c(1, 2)], pointindex = tv) %>%
    # Convert the polygons to a SpatialPolygons object
    SpatialPolygons(., proj4string = mesh$crs) %>%
    # Convert to sf
    st_as_sf(.)
  # Calculate and add area to the polygon
  dp$areakm2 <- st_area(dp) %>% units::set_units("km²") %>% as.numeric()
  dp$relarea <- dp$areakm2 / sum(dp$areakm2,na.rm = TRUE)
  return(dp)
}

#' Extract boundary points from mesh
#' @param mesh A [`inla.mesh`] object
#' @noRd
mesh_boundary <- function(mesh){
  assertthat::assert_that(inherits(mesh,'inla.mesh'))
  # Mesh coordinates
  loc <- mesh$loc
  loc[mesh$segm$int$idx[,2],]
}

#' Make Integration stack
#' @param mesh The background projection mesh
#' @param mesh.area The area of the mesh, has to match the number of integration points
#' @param cov The covariate data stack
#' @param pred_names The names of the used predictors
#' @param bdry The boundary used to create the mesh
#' @param resp The column with the Observed value
#' @noRd
inla_make_integration_stack <- function(mesh, mesh.area, cov, pred_names, bdry, resp = 'Observed' ){
  assertthat::assert_that(
    inherits(mesh,'inla.mesh'),
    length(mesh.area) == mesh$n, is.vector(mesh.area),
    is.data.frame(cov), assertthat::has_name(cov, c('x','y')),
    is.data.frame(bdry), is.character(resp)
  )

  # Get nearest average environmental data
  all_env <- get_ngbvalue(
    coords = mesh$loc[,1:2],
    env = cov,
    field_space = c('x','y'),
    longlat = raster::isLonLat(bdry)
  )
  # Get only target variables
  all_env <-  subset(all_env, select = pred_names)
  all_env$intercept <- 1 #FIXME: Unsure yet whether this is even necessary. Might kick out

  # Add diagonal for integration points
  idiag <- Matrix::Diagonal(mesh$n, rep(1, mesh$n))

  # Make some assertions
  stopifnot(assertthat::assert_that(
    mesh$n == nrow(all_env),
    nrow(idiag) == nrow(all_env)
  ))

  # Single response
  ll_resp <- list()
  ll_resp[[ resp ]] <- cbind( rep(0, mesh$n) ) # Provide NA to make data.frame
  ll_resp[['e']] <- as.numeric( mesh.area )

  # Effects list
  ll_effects <- list()
  # Note, order adding this is important apparently...
  ll_effects[['predictors']] <- all_env
  ll_effects[['intercept']] <- list(intercept = seq(1,mesh$n) )

  # Build integration stack of nearest predictors
  stk_int <- INLA::inla.stack(
    data    = ll_resp,
    A       = list(1, idiag),
    tag     = 'stk_int',
    effects = ll_effects
  )
  return(stk_int)
}



#' Function for creating a joint fitted and prediction stack
#' @param stk_resp A inla stack object
#' @param cov The covariate data stack
#' @param pred.names The predictors to use
#' @param mesh The background projection mesh
#' @param mesh.area The area calculate for the mesh
#' @param type Name to use
#' @param spde An spde field if specified
#' @noRd

inla_make_prediction_stack <- function(stk_resp, cov, pred.names, mesh, mesh.area, type, spde = NULL){
  # Security checks
  assertthat::assert_that(
    inherits(stk_resp, 'inla.data.stack'),
    inherits(mesh,'inla.mesh'),
    is.character(type),
    is.data.frame(cov),
    is.null(spde)  || 'spatial.field' %in% names(spde)
  )
  # need to create a new A projection matrix for the prediction data
  mat_pred <- INLA::inla.spde.make.A(mesh,
                             loc = as.matrix(cov[,1:2]) )

  # Single response
  ll_pred <- list()
  ll_pred[[ names(stk_resp$data$names)[1] ]] <- cbind( rep(NA, nrow(cov)) ) # Set to NA to predict for fitted area
  ll_pred[['e']] <- rep(0, nrow(cov))

  # Effects
  ll_effects <- list()
  # Note, order adding this is important apparently...
  ll_effects[['predictors']] <- cov[,pred.names]
  ll_effects[['intercept']] <- list(intercept = seq(1,mesh$n) ) # FIXME: Potential source for bug. Think name of intersects need to differ if multiple likelihoods specified
  if(!is.null(spde)) ll_effects[['intercept']] <- c(ll_effects[['intercept']], spde)
  # Define A
  A <- list(1, mat_pred )

  # Create stack depending on the number of variables in response
  # if( length(stk_resp$data$names[[1]]) > 1) {
  #   stop('TBD')
  #   ys <- cbind(rep(NA, nrow(pred.grid)), rep(NA, nrow(pred.grid)))
  #   stack.pred.response <- inla.stack(data=list(y=ys),
  #                                     effects = list(list(data.frame(interceptA=rep(1,np))), env = pred.grid$cov, list(uns_field=1:spde$n.spde)),
  #                                     A = A,
  #                                     tag = paste0('pred_', type))
  # } else if("Ntrials" %in% stk_resp$data$names) {
  #   stop('TBD')
  #   stack.pred.response <- inla.stack(data=list(y=NA, Ntrials = rep(1,np)),
  #                                     effects = list(list(data.frame(interceptA=rep(1,np))), env = pred.grid$cov,
  #                                                    list(spatial.field=1:spde$n.spde)),
  #                                     A = A,
  #                                     tag = paste0('pred_', type) )
  # } else {
    stk_pred <-
      INLA::inla.stack(
        data =  ll_pred,              # Response
        A = A,                        # Predictor projection matrix
        effects = ll_effects,         # Effects matrix
        tag = paste0('pred_',type) # New description tag
      )
  # }

  # TODO: Insert some security checks that stk_resp and stK_pred are equal in variable names
  # Final security checks
  # assertthat::assert_that(
  #   names(stk_resp$effects$names) %in% names(stk_pred$effects$names),
  #   msg = 'Response stak and prediction stack mismatch.'
  # )
  return(stk_pred)
}

#' Tidy up summary information from a INLA model
#'
#' @param m A trained INLA model object
#' @param what A [`character`]
#' @param ... Other options to based on
#' @noRd

#TODO: Lot more to add here, including options on what to extract
tidy_inla_summary <- function(m, what = 'fixed',...){
  assertthat::assert_that(
    inherits(m,'inla'),
    is.character(what),length(what)==1,
    what %in% c('fixed','fitted','random','spde2')
  )

  w1 <- grep('summary',names(m),value = TRUE) # Grep summary objects
  w2 <- grep(what, w1,value = TRUE) # Grep specific summary
  assertthat::assert_that(length(w2)==1)

  # Format the output
  m[[w2]] %>%
    tibble::rownames_to_column('variable') %>%
    tibble::as_tibble()
}

#' Plot marginal distributions of effects or hyperparameters from INLA model
#' @param A INLA model
#' @param what Either 'fixed' or 'hyper'
#' @noRd
plot_inla_marginals = function(inla.model, what = 'fixed'){
  assertthat::assert_that(inherits(inla.model,'inla'),
                          is.character(what),
                          what %in% c('fixed','hyper'))
  par(mfrow = c(4,4))
  if(what == 'fixed'){
    varnames <- names(inla.model$marginals.fixed)
    for(i in 1: length(varnames)){
      var.mar <- data.frame(inla.model$marginals.fixed[i])
      plot(x = var.mar[,1], y=var.mar[, 2], type="l",
           xlab=paste(names(var.mar)[1]), ylab=paste(names(var.mar)[2]))
      abline(v=0, col="red")
    }
  } else {
    varnames <- names(inla.model$marginals.hyperpar)
    for(i in 1: length(varnames)){
      var.mar <- data.frame(inla.model$marginals.hyperpar[i])
      plot(x = var.mar[,1], y=var.mar[, 2], type="l",
           xlab=paste(names(var.mar)[1]), ylab=paste(names(var.mar)[2]))
    }
  }
}

