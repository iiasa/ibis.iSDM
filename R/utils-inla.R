#' Built formula for INLA model
#'
#' @description This function built a formula for a `engine_inla()` model.
#'
#' @param model A [`list()`] object containing the prepared model data.
#' @param id The id for the species formula.
#' @param x A [`BiodiversityDistribution`] object.
#' @param settings A [`Settings`] object.
#'
#' @note Function is not meant to be run outside the train() call.
#'
#' @author Martin Jung
#'
#' @noRd
#'
#' @keywords internal
built_formula_inla <- function(model, id, x, settings){
  assertthat::assert_that(
    is.list(model),
    length(model) > 0,
    assertthat::has_name(model, "biodiversity"),
    assertthat::has_name(model, "predictors_names"),
    inherits(x, "BiodiversityDistribution"),
    inherits(settings, 'Settings'),
    msg = "Error in model object. This function is not meant to be called outside ouf train()."
  )

  # check if all intercepts are identical
  assertthat::assert_that(length(unique(sapply(model$biodiversity, function(i) i$use_intercept))) == 1,
                          msg = "'separate_intercept' must be identical for all datasets.")

  # Get object from model
  obj <- model$biodiversity[[id]]

  # Extract basic stats from the model object
  types <- as.character( sapply( model$biodiversity, function(x) x$type ) )
  fams <- as.character( sapply( model$biodiversity, function(z) z$family ) )
  bionames = sapply(model$biodiversity, function(x) x$name)
  ids <- names(model$biodiversity)
  priors <- model$priors

  if(x$get_engine() == "<INLA>"){
    # Default equation found
    if(obj$equation =='<Default>' || is.Waiver(obj$equation)){
      # Check potential for rw1 fits
      if(settings$get('only_linear') == FALSE){
        # Get Numeric variables
        vf <- obj$predictors_types$predictors[obj$predictors_types$type=="numeric"]
        var_rw1 <- vf
        # Set remaining variables to linear, especially if they are factors
        var_lin <- c()
        if(any(obj$predictors_types$type=="factor")){
          vf <- obj$predictors_types$predictors[obj$predictors_types$type=="factor"]
          var_lin <- c(var_lin, names(explode_factor(obj$predictors[[vf]], vf)) )
        } else {
          var_lin <- obj[['predictors_names']][which( obj[['predictors_names']] %notin% var_rw1 )]
        }
      } else {
        var_rw1 <- c()
        # Set remaining variables to linear
        var_lin <- obj$predictors_types$predictors[obj$predictors_types$type=="numeric"]
        # If any factors are present, split them and add too
        if(any( obj$predictors_types$type=="factor")){
          vf <- obj$predictors_types$predictors[ obj$predictors_types$type=="factor"]
          var_lin <- c(var_lin, names(explode_factor( obj$predictors[[vf]], vf)) )
        }
      }

      # Check whether to use dataset specific intercepts
      if(length(types)>1 && obj$use_intercept){
        ii <- paste0('Intercept_', make.names(tolower(bionames)),'_', types,
                     collapse = " + ")
        ii <- paste0("0 + ", ii) # hack to remove global
      } else {ii <- ""}

      # Go through each variable and build formula for likelihood
      form <- paste0(paste0("observed ~ Intercept + ", ii))

      # remove global intercept term
      if(ii != "") form <- gsub(x = form, pattern = "Intercept \\+ ", replacement = "")

      # Check whether priors have been specified and if yes, use those
      if(!is.Waiver(priors)){
        # Loop through all provided INLA priors
        supplied_priors <- priors$ids()

        for(v in supplied_priors){
          # Prior variable name
          vn <- as.character( priors$varnames()[v] )
          if(vn == 'spde') next()
          # Prior variable type
          vt <- as.character( priors$types()[v] )
          if(vt == 'clinear'){
            # Constrained linear effect
            form <- paste(form, '+', paste0('f(', vn, ', model = \'clinear\', ',
                                            'range = c(', priors$get(vn)[1],',', priors$get(vn)[2],') )',
                                            collapse = ' + ' ) )
          } else if(vt %in% c('normal','gaussian')) {
            # Add linear effects
            form <- paste(form, '+', paste0('f(', vn, ', model = \'linear\', ',
                                            'mean.linear = ', priors$get(vn)[1],', ',
                                            'prec.linear = ', priors$get(vn)[2],')',
                                            collapse = ' + ' ) )
          } else if(vt == 'pc.prec' || vt == 'loggamma'){
            # Add RW effects with pc priors. PC priors is on the KL distance (difference between probability distributions), P(sigma >2)=0.05
            # Default is a loggamma prior with mu 1, 5e-05. Better would be 1, 0.5 following Caroll 2015
            form <- paste0(form, '+', paste0('f(INLA::inla.group(', vn, '), model = \'rw1\', ',
                                             # 'scale.model = TRUE,',
                                             'hyper = list(theta = list(prior = ',vt,', param = c(',priors$get(vn)[1],',',priors$get(vn)[2],')) )
                                                 )',collapse = ' + ')
            )
          }
        }
        # Add linear for those missed ones
        miss <- c(var_lin, var_rw1)[c(var_lin, var_rw1) %notin% priors$varnames()]
        if(length(miss)>0){
          if(any(miss %in% var_lin)){
            # Add linear predictors without priors
            form <- paste(form, ' + ',
                          paste('f(', miss[which(miss %in% var_lin)],', model = \'linear\')', collapse = ' + ')
            )
          }
          if(length(var_rw1)>0 & (any(miss %in% var_rw1))){
            # Random walk where feasible and not already included
            form <- paste(form, ifelse(length(var_lin) == 0,'+',''), paste('f(INLA::inla.group(', miss[which(miss %in% var_rw1)],'),',
                                                                           # 'scale.model = TRUE, ',
                                                                           # Add RW effects with pc priors. PC priors is on the KL distance (difference between probability distributions), P(sigma >2)=0.05
                                                                           # Default is a loggamma prior with mu 1, 5e-05. Better would be 1, 0.5 following Caroll 2015
                                                                           'hyper = list(theta = list(prior = \'loggamma\', param = c(1, 0.5))),',
                                                                           'model = \'rw1\')', collapse = ' + ' ) )
          }
        }
      } else {
        # No priors specified, simply add variables with default
        # Linear for those with few observations
        if(length(var_lin)>0){
          form <- paste0(form, paste0('f(', var_lin,', model = \'linear\')', collapse = ' + '))
        }
        # Random walk where feasible
        if(length(var_rw1)>0){
          form <- paste(form, ifelse(length(var_lin) == 0, '+', ''),
                        paste('f(INLA::inla.group(', var_rw1,'),',
                              # 'scale.model = TRUE,',
                              # Add RW effects with pc priors. PC priors is on the KL distance
                              # (difference between probability distributions), P(sigma >2)=0.05
                              # Default is a loggamma prior with mu 1, 5e-05. Better would be 1, 0.5 following Caroll 2015
                              'hyper = list(theta = list(prior = \'loggamma\', param = c(1, 0.5))),',
                              'model = \'rw1\')', collapse = ' + ' ) )
        }
      }
      form <- to_formula(form) # Convert to formula
      # Add offset if specified
      if(!is.Waiver(x$offset) ){ form <- stats::update.formula(form, paste0('~ . + offset(spatial_offset)') ) }
      if( length( grep('Spatial', x$get_latent() ) ) > 0 ){
        if(attr(x$get_latent(), "method") != "poly"){
          # Update with spatial term
          form <- stats::update.formula(form, paste0(" ~ . + ",
                                                     x$engine$get_equation_latent_spatial(
                                                       method = attr(x$get_latent(),'method'),
                                                       vars = which(ids == id),
                                                       separate_spde = attr(x$get_latent(),'separate_spde')
                                                     )
          )
          )
        }
      }
    } else{

      # MH: Add intercept stuff as well

      # If custom supplied formula, check that variable names match supplied predictors
      if(getOption('ibis.setupmessages', default = TRUE)) myLog('[Estimation]','yellow','Use custom model equation')
      form <- to_formula( obj$equation )
      # If response is missing, add manually
      if(attr(stats::terms(form), "response")==0){
        form <- stats::update.formula(form, "observed ~ .")
      }
      # Security checks
      assertthat::assert_that(
        is.formula(form),
        attr(stats::terms(form), "response")==1, # Has Response
        all( all.vars(form) %in% c('observed', obj[['predictors_names']]) )
      )
    }
    # -------------------- INLABRU formula----------------
  } else if(x$get_engine() == "<INLABRU>"){

    # Default equation found (e.g. no separate specification of effects)
    if(obj$equation=='<Default>'){
      # Check whether to use dataset specific intercepts
      if(length(types)>1 && obj$use_intercept){
        ii <- paste0('Intercept_', make.names(tolower(bionames)),'_', types,
                     collapse = " + ")
        ii <- paste0(" - 1 + ", ii) # hack to remove global
      } else {ii <- ""}

      # Go through each variable and build formula for likelihood
      form <- to_formula(paste0(paste0("observed ~ Intercept + ",
                                       paste0(obj$predictors_names, collapse = " + ")),
                                ii)) # adding separated intercepts or empty string

      # Add offset if specified TODO: Not quite sure if this formulation works
      # for inlabru predictor expressions
      if(!is.Waiver(x$offset) ){ form <- stats::update.formula(form, paste0('~ . + offset(spatial_offset)') ) }
      if( length( grep('Spatial', x$get_latent() ) ) > 0 ){
        if(attr(x$get_latent(), "method") != "poly"){
          # Update with spatial term
          form <- stats::update.formula(form, paste0(" ~ . + ",
                                                     # For SPDE components, simply add spatial.field
                                                     paste0("spatial.field", which(ids == id))
          )
          )
        }
      }
    } else {
      # If custom likelihood formula is provided, check that variable names
      # match supplied predictors
      form <- obj$equation
      assertthat::assert_that(
        all( all.vars(form) %in% c('observed', obj$predictors_names) )
      )

      # Convert to formula to be safe
      form <- to_formula( obj$equation )
      # Add generic Intercept if not set in formula
      if("Intercept" %notin% all.vars(form)) form <- stats::update.formula(form, ". ~ Intercept + .")
      # If length of ids is larger than 1, add dataset specific intercept too
      # Check whether to use dataset specific intercepts
      if(length(types)>1 && obj$use_intercept){
        ii <- paste0('Intercept_', make.names(tolower(bionames)),'_', types,
                     collapse = " + ")
        ii <- paste0(" - 1 + ", ii) # hack to remove global
      } else {ii <- ""}

      # get all predictors of formulat
      rhs <- all.vars(form)[all.vars(form) != "observed"]

      # update formula adding intercepts etc. (see above <Default>)
      form <- to_formula(paste0(paste0("observed ~ ", paste0(rhs, collapse = " + ")), ii))

      if( length( grep('Spatial',x$get_latent() ) ) > 0 ){
        if(attr(x$get_latent(), "method") != "poly"){
          # Update with spatial term
          form <- stats::update.formula(form, paste0(" ~ . + ",
                                                     # For SPDE components, simply add spatial.field
                                                     paste0("spatial.field",which(ids == id))
          )
          )
        }
      }

    }
  }
  return(form)
}

#' Calculate area of each voronoi polygon in a INLA mesh
#'
#' @param mesh [`inla.mesh`] mesh object.
#' @param region.poly A supplied [`region.poly`] object.
#' @param variant A character to which type of area calculation (Default: \code{'gpc'}).
#' @param relative Should the total amount of area converted to relatives (Default: \code{FALSE}).
#'
#' @returns A [`vector`] with the area of each polygon.
#'
#' @keywords utils
#'
#' @noRd
#'
#' @keywords internal
mesh_area = function(mesh, region.poly = NULL, variant = 'gpc', relative = FALSE){
  assertthat::assert_that(inherits(mesh,'inla.mesh'),
                          is.null(region.poly) || inherits(region.poly,'Spatial'),
                          is.character(variant)
  )
  check_package("deldir")
  # Function from SDraw
  voronoi.polygons <- function (x, bounding.polygon = NULL, range.expand = 0.1)
  {
    if (!inherits(x, "SpatialPoints")) {
      cli::cli_abort("Must pass a SpatialPoints* object to voronoi.polygons.")
    }
    crds = sp::coordinates(x)
    if (is.null(bounding.polygon)) {
      if (length(range.expand) == 1) {
        range.expand <- rep(range.expand, 2)
      }
      else if (length(range.expand) > 2) {
        cli::cli_alert_warning("Only first two elements of range.expand used in voronoi.polygons")
        range.expand <- range.expand[1:2]
      }
      dxy <- diff(c(t(sp::bbox(x))))[c(1, 3)]
      bb <- sp::bbox(x) + (matrix(dxy, nrow = 2, ncol = 1) %*%
                             matrix(c(-1, 1), nrow = 1, ncol = 2)) * abs(range.expand)
      bb <- c(t(bb))
    } else {
      bb = c(t(sp::bbox(bounding.polygon)))
    }
    z = deldir::deldir(crds[, 1], crds[, 2], rw = bb)
    w = deldir::tile.list(z)
    polys = vector(mode = "list", length = length(w))
    for (i in seq(along = polys)) {
      pcrds = cbind(w[[i]]$x, w[[i]]$y)
      # 01/11/2022 -> Added as bug fix
      if(nrow(pcrds)==0) next()
      pcrds = rbind(pcrds, pcrds[1, ])
      polys[[i]] = sp::Polygons(list(sp::Polygon(pcrds)), ID = as.character(i))
    }
    # 01/11/2022 -> Error removal
    if(length(which(sapply(polys, is.null)))>0) {
      crds <-  crds[-which(sapply(polys, is.null)),]
      polys[which(sapply(polys, is.null))] <- NULL
    }
    SP = sp::SpatialPolygons(polys, proj4string = sp::CRS(sp::proj4string(x)))
    voronoi = sp::SpatialPolygonsDataFrame(SP, data = data.frame(x = crds[,1],
                                                                 y = crds[, 2],
                                                                 area = sapply(methods::slot(SP, "polygons"), methods::slot, "area"),
                                                                 row.names = sapply(methods::slot(SP, "polygons"), methods::slot, "ID")))
    if (!is.null(bounding.polygon)) {
      bounding.polygon <- rgeos::gUnion(bounding.polygon, bounding.polygon)
      voronoi.clipped <- rgeos::gIntersection(voronoi, bounding.polygon,
                                              byid = TRUE, id = row.names(voronoi))
      df <- data.frame(voronoi)
      df$area <- sapply(methods::slot(voronoi.clipped, "polygons"),
                        methods::slot, "area")
      voronoi <- sp::SpatialPolygonsDataFrame(voronoi.clipped,
                                              df)
    }
    voronoi
  }
  # Precalculate the area of each
  # Get areas for Voronoi tiles around each integration point
  dd <- deldir::deldir(mesh$loc[,1], mesh$loc[,2])
  tiles <- deldir::tile.list(dd)

  if(variant == 'gpc'){
    # Try to convert to spatial already
    if(!inherits(region.poly, 'Spatial')) region.poly <- methods::as(region.poly,'Spatial')

    poly.gpc <- methods::as(region.poly@polygons[[1]]@Polygons[[1]]@coords,'gpc.poly')
    w <- sapply(tiles, function(p) rgeos::area.poly(rgeos::intersect(methods::as(cbind(p$x, p$y), 'gpc.poly'), poly.gpc)))
    if(relative) w <- w / sum(w)
  } else if (variant == 'gpc2'){
    # Try to convert to spatial already
    if(!inherits(region.poly, 'Spatial')) region.poly <- methods::as(region.poly,'Spatial')
    tiles <- voronoi.polygons(sp::SpatialPoints(mesh$loc[, 1:2]))
    w <- sapply(1:length(tiles), function(p) {
      aux <- tiles[p, ]

      if(rgeos::gIntersects(aux, region.poly) ) {
        return(rgeos::gArea(rgeos::gIntersection(aux, region.poly)))
      } else {
        return(0)
      }
    })
    assertthat::assert_that(assertthat::are_equal(round( sum(w) ), round( rgeos::gArea(region.poly) ) )) # Security check
    if(relative) w <- w / sum(w)
  } else {
    # Convert to Spatial Polygons
    polys <- sp::SpatialPolygons(lapply(1:length(tiles), function(i) {
      p <- cbind(tiles[[i]]$x, tiles[[i]]$y)
      n <- nrow(p)
      sp::Polygons(list(sp::Polygon(p[c(1:n, 1), ])), i)
    }),proj4string = mesh$crs)

    # Calculate area of each polygon in km2
    w <- sf::st_area(
      sf::st_as_sf(polys)
    ) |> units::set_units(km^2) |> as.numeric()
    # Relative area
    if(relative) w <- w / sum(w)
  }

  assertthat::assert_that(is.vector(w))
  return(w)
}

#' Mesh to polygon script
#'
#' @param mesh [`inla.mesh`] mesh object.
#'
#' @returns A [`sf`] object.
#'
#' @keywords utils
#'
#' @noRd
#'
#' @keywords internal
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
                 sp::Polygons(list(
                   sp::Polygon( points[c(cur, cur[1]), ], hole = FALSE)),
                   ID = index
                 )
               }, points = mesh$loc[, c(1, 2)], pointindex = tv) |>
    # Convert the polygons to a SpatialPolygons object
    sp::SpatialPolygons(., proj4string = mesh$crs) |>
    # Convert to sf
    sf::st_as_sf(.)
  # Calculate and add area to the polygon
  dp$areakm2 <- sf::st_area(dp) |> units::set_units(km^2) |> as.numeric()
  dp$relarea <- dp$areakm2 / sum(dp$areakm2,na.rm = TRUE)
  return(dp)
}

#' Extract boundary points from mesh
#'
#' @param mesh A [`inla.mesh`] object.
#'
#' @keywords utils
#'
#' @noRd
#'
#' @keywords internal
mesh_boundary <- function(mesh){
  assertthat::assert_that(inherits(mesh,'inla.mesh'))
  # Mesh coordinates
  loc <- mesh$loc
  loc[mesh$segm$int$idx[,2],]
}

#' Create a barrier representation of a mesh
#'
#' @description
#' **Work in progress* Creating a physical barrier model for INLA
#'
#' @param mesh A [`inla.mesh`] object.
#' @param region.poly A [`SpatialPolygons`] object.
#'
#' @source https://www.sciencedirect.com/science/article/pii/S221167531830099X
#'
#' @keywords utils
#'
#' @noRd
#'
#' @keywords internal
mesh_barrier <- function(mesh, region.poly){
  assertthat::assert_that(
    inherits(mesh,'inla.mesh'),
    inherits(region.poly,'SpatialPolygons')
  )
  # Get number of boundary trianbles on graph
  tl <- length(mesh$graph$tv[,1])
  # Define position matrix for triangles
  posTri <- matrix(0, tl, 2)

  # Loop through triangle
  for (t in 1:tl){
    temp <- mesh$loc[mesh$graph$tv[t, ], ]
    posTri[t,] <- colMeans(temp)[c(1,2)]
  }

  posTri <- sp::SpatialPoints(posTri)
  sp::proj4string(posTri) <- sp::proj4string(mesh$crs)

  # Overlay with background
  ovl <- sp::over(region.poly, posTri, returnList=T)
  ovl <- unlist(ovl)
  barrier.triangles <- setdiff(1:tl, ovl)

  # Define INLA barrier polygon
  suppressMessages(
    suppressWarnings(
      poly.barrier <- INLA::inla.barrier.polygon(mesh, barrier.triangles)
    )
  )

  return(poly.barrier)
}

#' Query if a point is inside the mesh boundary
#'
#' @param mesh A [`inla.mesh`] object.
#' @param coords Either a two-column [`data.frame`] or [`matrix`] of coordinates.
#' Alternatively a [`Spatial`] or [`sf`] object from which coordinates can be extracted.
#'
#' @return A [`vector`] of Boolean values indicating if a point is inside the mesh.
#'
#' @keywords utils
#'
#' @noRd
#'
#' @keywords internal
coords_in_mesh <- function(mesh, coords) {
  assertthat::assert_that(
    inherits(mesh,'inla.mesh'),
    inherits(coords,'sf') || inherits(region.poly,'Spatial') || inherits(coords, 'matrix') || inherits(coords, 'data.frame')
  )

  # Get coordinates depending on type
  if (inherits(coords, "Spatial")) {
    loc <- sp::coordinates(coords)
  } else if(inherits(coords, 'sf')){
    loc <- sf::st_coordinates(coords)
  } else if(inherits(coords, 'matrix') || inherits(coords, 'data.frame') ){
    assertthat::assert_that(ncol(coords)==2, msg = 'Supplied coordinate matrix is larger than 2 columns.')
    loc <- coords[,c(1,2)]
  }
  loc <- as.matrix(loc)

  loc_inside <- INLA::inla.fmesher.smorg(
    loc = mesh$loc,
    tv = mesh$graph$tv,
    points2mesh = loc
  )
  # Return vector with point not in
  return(loc_inside$p2m.t[,1]!=0)
}

#' Manual prediction by matrix multiplication
#'
#' @description
#' Spatial predictions with INLA can be quite computationally costly. Assuming
#' that model coefficients are fixed and linear, it is possible to obtain
#' comparable predictions simply by matrix multiplication.
#'
#' TODO: Switch to posterior sampling
#' https://groups.google.com/g/r-inla-discussion-group/c/y-rQlDVtzmM
#'
#' @param mesh x A [`distribution`] object used for fitting an INLA model.
#' @param mod A trained [`distribution`] model.
#' @param type The summary statistic to use.
#' @param backtransf Either NULL or a function.
#' @param coords A [matrix] with coordinates or \code{NULL}. If \code{NULL}
#' coordinates are recreated from predictors.
#'
#' @keywords utils
#'
#' @noRd
#'
#' @keywords internal
coef_prediction <- function(mesh, mod, type = 'mean',
                            backtransf = NULL,
                            coords = NULL){
  assertthat::assert_that(!missing(mesh), !missing(mod),
                          inherits(mesh,'inla.mesh'),
                          inherits(mod,'INLA-Model'),
                          (ncol(coords)==2) || is.null(coords),
                          msg = 'This function only works for INLA models')
  # Check whether data is all there
  assertthat::assert_that(
    utils::hasName(mod,'model'),
    # length( grep('stk',x$engine$list_data()) )>0,
    is.function(backtransf) || is.null(backtransf),
    msg = 'Not all data and parameters needed for prediction are present!'
  )

  # Get formula
  form <- mod$get_data('fit_best_equation')
  if(length(grep(pattern = '\\*',deparse(form)))) cli::cli_abort('Interactions are not supported!')
  # Check whether any rw1 effects are in the formula. If so return error
  te <- attr(stats::terms.formula(form),'term.label')
  if(length(grep(pattern = '\"rw',x = te))>0) cli::cli_abort('This function does not work with INLA rw effects!')

  # Manual prediction
  model <- mod$get_data('fit_best')
  mesh <- mod$get_data('mesh')
  # Covariates for prediction points
  preds <- mod$model$predictors_object
  ofs <- mod$model$offset

  # Some checks between models and data
  assertthat::assert_that(is.character(type),
                          type %in% names(model$summary.fixed),
                          all( rownames(model$summary.fixed) %in% names(preds) )
  )
  if(type !='mean') cli::cli_alert_warning('Predictions other than __mean__ unlikely to work well...!')

  # Output raster
  if(is.null(coords)) coords <- preds[,c('x','y')]
  if(nrow(coords)!= nrow(preds)){
    # Recalculate average predictors for new coordinates
    preds <- get_rastervalue(coords = coords,
                             env = preds,
                             rm.na = FALSE
    )
    # Set Intercept variables
    for(val in grep('Intercept',rownames(model$summary.fixed),value = TRUE)){
      preds[[val]] <- 1
    }
  }
  temp = terra::rast(coords, type = "xyz")

  # remake the A matrix for prediction
  Aprediction <- INLA::inla.spde.make.A(mesh = mesh,
                                        loc = as.matrix( coords ))

  # make empty matrix to fill predictions
  out <- matrix(0, nrow = dim(Aprediction)[1], ncol = 1)
  assertthat::assert_that(nrow(out)==nrow(preds))
  for( n in rownames(model$summary.fixed) ) {
    #TODO: If nrow(out) differ with preds, calculate nearest neighbours
    out[,1] <- out[,1] + (model$summary.fixed[n, type] %*% preds[,n])
  }

  # create the spatial structure if existing
  if( length(model$summary.random) >0){
    assertthat::assert_that(length(model$summary.random) == 1, # FIXME: If multiple spatial latent effects this needs adapting
                            'spatial.field1' %in% names(model$summary.random),
                            msg = 'Spatial random effect wrongly specified!')
    # Check that type is present, otherwise use 'mean'
    sfield_nodes <- model$summary.random[[1]][,type]
    field <- (Aprediction %*% as.data.frame(sfield_nodes)[, 1] )
    out <- out + field
  }

  # Add offset if specified
  if(length(grep('offset\\(', deparse1(form)))>0 && !is.Waiver(ofs)){
    if(nrow(coords)!= nrow(ofs)){
      # Recalculate average predictors for new coordinates
      ofs <- get_ngbvalue(coords = coords,
                          env = ofs,
                          longlat = terra::is.lonlat(mesh$crs),
                          field_space = c('x','y'))
    }
    # ofs[is.na(ofs[,3]),3] <- 0
    out <- out + ofs[,3]
  }

  # Fill output raster
  temp[] <- out[, 1]
  if(!is.null(backtransf)) temp <- terra::app(temp, backtransf)
  # plot(temp, col = cols)
  return( temp )
}

#' Direct prediction by posterior simulation
#'
#' @param mod A trained distribution model.
#' @param nsamples [`numeric`] on the number of samples to be taken from the posterior.
#' @param backtransf Either \code{NULL} or a function.
#' @param seed A random seed that can be specified.
#'
#' @keywords utils
#'
#' @noRd
#'
#' @keywords internal
post_prediction <- function(mod, nsamples = 100,
                            backtransf = NULL,
                            seed = 0){
  assertthat::assert_that(!missing(mod),
                          inherits(mod,'DistributionModel'),
                          inherits(mod,'INLA-Model'),
                          msg = 'This function only works for INLA models')
  # Check whether data is all there
  assertthat::assert_that(
    utils::hasName(mod,'model'),
    # length( grep('stk',x$engine$list_data()) )>0,
    is.function(backtransf) || is.null(backtransf),
    is.numeric(seed) || is.null(seed),
    msg = 'Not all data and parameters needed for prediction are present!'
  )
  check_package("inlabru")

  # Get data from object #
  model <- mod$get_data('fit_best')
  mesh <- mod$get_data('mesh')
  # Get formula
  form <- mod$get_data('fit_best_equation')
  # Some checks on the form
  if(length(grep(pattern = '\\*',deparse(form)))) cli::cli_abort('Interactions are not (yet) supported!')
  # Check whether any rw1 effects are in the formula. If so return error
  te <- attr(stats::terms.formula(form),'term.label')
  if(length(grep(pattern = '\"rw',x = te))>0) cli::cli_abort('This function does not work with INLA rw effects!')

  # Covariates for prediction points
  preds <- mod$model$predictors
  # Set any other existing intercept variables
  preds[,grep('Intercept',rownames(model$summary.fixed),value = TRUE)] <- 1
  preds <- sp::SpatialPixelsDataFrame(preds[,c('x','y')],data=preds)
  preds_names <- mod$model$predictors_names
  preds_types <- mod$model$predictors_types
  ofs <- mod$model$offset

  # See:
  # https://groups.google.com/g/r-inla-discussion-group/c/TjRwP6tB0nk/m/FSD39whVBgAJ
  # https://groups.google.com/g/r-inla-discussion-group/c/Lw2fI-u-EvU/m/rB_-gwWAAgAJ
  # --- #
  # Simulate from approximated posterior
  samples <- INLA::inla.posterior.sample(n = nsamples,
                                         result = model,
                                         # seed = seed,
                                         use.improved.mean = TRUE, # marginal means?
                                         num.threads = ifelse(getOption("ibis.nthread")>1, getOption("ibis.nthread"),NULL),
                                         parallel.configs = TRUE,
                                         verbose = TRUE
  )

  # inlabru extract entries function and names standardization functions
  standardise.names <- function (x) {
    new_names <- vapply(x, function(x) {
      gsub("[-() ]", "_", x = x, fixed = FALSE)
    }, "name")
    not_ok <- grepl("__", x = new_names)
    while (any(not_ok)) {
      new_names[not_ok] <- vapply(new_names[not_ok], function(x) {
        gsub("__", "_", x = x, fixed = FALSE)
      }, "name")
      not_ok <- grepl("__", x = new_names)
    }
    new_names
  }

  # Format the prior outputs to get the predictor values
  ssmpl <- list()
  for (i in seq_along(samples)) {
    smpl.latent <- samples[[i]]$latent
    smpl.hyperpar <- samples[[i]]$hyperpar
    vals <- list()

    # Extract simulated predictor and fixed effects
    for (name in unique(c("Predictor", preds_names))) {
      vals[[name]] <- inlabru:::extract.entries(name, smpl.latent)
    }
    # Remove any variables with no values (removed during model fit)
    vals <- vals[which(lapply(vals, length)>0)]

    # For fixed effects that were modelled via factors we attach an extra vector holding the samples
    fac.names <- subset(preds_types, type == 'factor')
    if(nrow(fac.names)>0){
      for (name in fac.names$predictors) {
        vals[[name]] <- smpl.latent[startsWith(rownames(smpl.latent), name), ]
        names(vals[[name]]) <- lapply(names(vals[[name]]), function(nm) {substring(nm, nchar(name) + 1)})
      }
    }

    # Extract simulated latent variables.
    # If the model is "clinear", however, we might extract the realisations
    # from the hyperpar field. TODO: check if all the special models now have
    # their results available as latent random effects, and avoid special code,
    # since the hyperpar name definition has changed
    if (length(model$summary.random) > 0) {
      for (k in seq_along(model$summary.random)) {
        name <- unlist(names(model$summary.random[k]))
        vals[[name]] <- inlabru:::extract.entries(name, smpl.latent)
      }
    }
    if (length(smpl.hyperpar) > 0) {
      ## Sanitize the variable names; replace problems with "_".
      ## Needs to handle whatever INLA uses to describe the hyperparameters.
      ## Known to include " " and "-" and potentially "(" and ")".
      names(smpl.hyperpar) <- standardise.names(names(smpl.hyperpar))
    }
    ssmpl[[i]] <- c(vals, smpl.hyperpar)
  }
  vals <- ssmpl;rm(ssmpl)
  # Equivalent of inlabru:::inla.posterior.sample.structured
  myLog('[Summary]','green',paste('Formatted', length(vals), 'posterior samples'))

  # evaluate_model Function
  A <- inlabru:::ibm_amatrix(model, data = preds)
  A <- x$engine$data$stk_pred$stk_proj$A

  effects <- inlabru::evaluate_effect_multi_state(
    model$effects[included],
    state = vals,
    data = preds,
    A = A
  )


  if (is.null(predictor)) {
    return(effects)
  }

  values <- inlabru:::evaluate_predictor(
    model,
    state = state,
    data = data,
    effects = effects,
    predictor = predictor,
    format = format
  )


  # --- #
  # Summarise Functions
  expand_to_dataframe <- function (x, data = NULL) {
    if (is.null(data)) {
      data <- data.frame(matrix(nrow = NROW(x), ncol = 0))
    }
    only_x <- setdiff(names(x), names(data))
    if (length(only_x) < length(names(x))) {
      x <- x[!(names(x) %in% names(data))]
    }
    if (inherits(x, "SpatialPixels") && !inherits(x, "SpatialPixelsDataFrame")) {
      result <- sp::SpatialPixelsDataFrame(x, data = data)
    }
    else if (inherits(x, "SpatialGrid") && !inherits(x, "SpatialGridDataFrame")) {
      result <- sp::SpatialGridDataFrame(x, data = data)
    }
    else if (inherits(x, "SpatialLines") && !inherits(x, "SpatialLinesDataFrame")) {
      result <- sp::SpatialLinesDataFrame(x, data = data)
    }
    else if (inherits(x, "SpatialPolygons") && !inherits(x, "SpatialPolygonsDataFrame")) {
      result <- sp::SpatialPolygonsDataFrame(x, data = data)
    }
    else if (inherits(x, "SpatialPoints") && !inherits(x, "SpatialPointsDataFrame")) {
      result <- sp::SpatialPointsDataFrame(x, data = data)
    }
    else if (inherits(x, "Spatial")) {
      result <- sp::cbind.Spatial(x, data)
    }
    else {
      result <- cbind(x, data)
    }
    result
  }
  post_summarize <- function(data, x = NULL, cbind.only = FALSE) {
    if (is.list(data)) {
      data <- do.call(cbind, data)
    }
    if (cbind.only) {
      smy <- data.frame(data)
      colnames(smy) <- paste0("sample.", 1:ncol(smy))
    }
    else {
      smy <- data.frame(apply(data, MARGIN = 1, mean, na.rm = TRUE),
                        apply(data, MARGIN = 1, stats::sd, na.rm = TRUE),
                        t(apply(data,MARGIN = 1, quantile, prob = c(0.025, 0.5, 0.975),na.rm = TRUE)),
                        apply(data, MARGIN = 1, min, na.rm = TRUE),
                        apply(data, MARGIN = 1, max, na.rm = TRUE))
      colnames(smy) <- c("mean", "sd", "q0.025",
                         "median", "q0.975", "smin", "smax")
      smy$cv <- smy$sd/smy$mean
      smy$var <- smy$sd^2
    }
    if (!is.null(x)) {
      smy <- expand_to_dataframe(x, smy)
    }
    return(smy)
  }

  drop <- FALSE # FIXME: Make parameter after finding out what this does

  if(is.data.frame(vals[[1]])){
    vals.names <- names(vals[[1]])
    covar <- intersect(vals.names, names(preds))
    estim <- setdiff(vals.names, covar)
    smy <- list()

    for (nm in estim) {
      smy[[nm]] <- post_summarize(
        lapply(
          vals,
          function(v) v[[nm]]
        ),
        x = vals[[1]][, covar, drop = FALSE]
      )
    }
    is.annot <- vapply(names(smy), function(v) all(smy[[v]]$sd == 0), TRUE)
    annot <- do.call(cbind, lapply(smy[is.annot], function(v) v[, 1]))
    smy <- smy[!is.annot]
    if (!is.null(annot)) {
      smy <- lapply(smy, function(v) cbind(data.frame(annot), v))
    }

    if (!drop) {
      smy <- lapply(
        smy,
        function(tmp) {
          if (NROW(preds) == NROW(tmp)) {
            expand_to_dataframe(preds, tmp)
          } else {
            tmp
          }
        }
      )
    }

    if (length(smy) == 1) smy <- smy[[1]]
  } else if(is.list(vals[[1]])) {
    vals.names <- names(vals[[1]])
    if (any(vals.names == "")) {
      cli::cli_alert_warning("Some generated list elements are unnamed")
    }
    smy <- list()
    for(nm in vals.names) {
      tmp <- post_summarize(
        lapply(
          vals,
          function(v) v[[nm]]
        )
      )
      if(!drop &&
         (NROW(preds) == NROW(tmp))) {
        smy[[nm]] <- expand_to_dataframe(preds, tmp)
      } else {
        smy[[nm]] <- tmp
      }
    }
  } else {
    tmp <- post_summarize(vals)
    if (!drop &&
        (NROW(preds) == NROW(tmp))) {
      smy <- expand_to_dataframe(preds, tmp)
    } else {
      smy <- tmp
    }
  }
  # Multiply ?

  # Output raster
  if(is.null(coords)) coords <- preds[,c('x','y')]
  if(nrow(coords)!= nrow(preds)){
    # Recalculate average predictors for new coordinates
    preds <- get_ngbvalue(coords = coords,
                          env = preds,
                          longlat = terra::is.lonlat(mesh$crs),
                          field_space = c('x','y'))
  }
  temp = terra::rast(coords, type = "xyz")

  # remake the A matrix for prediction
  Aprediction <- INLA::inla.spde.make.A(mesh = mesh,
                                        loc = as.matrix( coords ))

  # make empty matrix to fill predictions
  out <- matrix(0, nrow = dim(Aprediction)[1], ncol = 1)
  assertthat::assert_that(nrow(out)==nrow(preds))
  for( n in rownames(model$summary.fixed) ) {
    # If nrow(out) differ with preds, calculate nearest neighbours
    out[,1] <- out[,1] + (model$summary.fixed[n, type] %*% preds[,n])
  }

  # create the spatial structure if existing
  if( length(model$summary.random) >0){
    assertthat::assert_that(length(model$summary.random) == 1, # FIXME: If multiple spatial latent effects this needs adapting
                            'spatial.field1' %in% names(model$summary.random),
                            msg = 'Spatial random effect wrongly specified!')
    # Check that type is present, otherwise use 'mean'
    sfield_nodes <- model$summary.random[[1]][,type]
    field <- (Aprediction %*% as.data.frame(sfield_nodes)[, 1] )
    out <- out + field
  }

  # Add offset if specified
  if(length(grep('offset\\(', deparse1(form)))>0 && !is.Waiver(ofs)){
    if(nrow(coords)!= nrow(ofs)){
      # Recalculate average predictors for new coordinates
      ofs <- get_ngbvalue(coords = coords,
                          env = ofs,
                          longlat = terra::is.lonlat(mesh$crs),
                          field_space = c('x','y'))
    }
    out <- out + ofs[,"spatial_offset"]
  }

  # Fill output raster
  temp[] <- out[, 1]
  if(!is.null(backtransf)){
    temp <- terra::app(temp, backtransf)
  }
  # plot(temp, col = cols)
  return( temp )
}

#' Make Integration stack
#'
#' @param mesh The background projection mesh.
#' @param mesh.area The area of the mesh, has to match the number of integration
#' points.
#' @param model A prepared model object.
#' @param id A id supplied to name this object.
#' @param joint Whether a model with multiple likelihood functions is to be specified.
#'
#' @keywords utils
#'
#' @noRd
#'
#' @keywords internal
inla_make_integration_stack <- function(mesh, mesh.area, model, id, joint = FALSE){
  assertthat::assert_that(
    inherits(mesh,'inla.mesh'),
    length(mesh.area) == mesh$n, is.vector(mesh.area),
    is.list(model) && ("predictors_names" %in% names(model)),
    is.logical(joint)
  )
  # Get from model object everything we need
  cov <- model$predictors_object
  pred_names <- model$predictors_names
  bdry <- model$background

  # Get nearest average environmental data
  all_env <- get_rastervalue(
    coords = mesh$loc[,1:2],
    env = cov$get_data(df = FALSE),
    rm.na = FALSE
  )
  # Get only target variables
  all_env$Intercept <- 1

  # Add diagonal for integration points
  idiag <- Matrix::Diagonal(mesh$n, rep(1, mesh$n))

  # Make some assertions
  stopifnot(assertthat::assert_that(
    mesh$n == nrow(all_env),
    nrow(idiag) == nrow(all_env)
  ))

  # Single response
  ll_resp <- list()

  # Add the expected estimate and observed note
  # FIXME: Currently only two likelihoods are supported (binomial/poisson) with the NA order being the determining factor
  if(joint) ll_resp[[ 'observed' ]] <- cbind(rep(0, mesh$n), NA )
  if(!joint) ll_resp[[ 'observed' ]] <- cbind( rep(0, mesh$n) )
  ll_resp[[ 'e' ]] <- as.numeric( mesh.area )

  # Effects list
  ll_effects <- list()
  # Note, order adding this is important apparently...
  # ll_effects[['Intercept']] <- rep(1, nrow(all_env)) # Added 05/09
  ll_effects[['predictors']] <- all_env
  ll_effects[['spatial.field1']] <- list(spatial.field1 = seq(1, mesh$n) ) # Changed to spatial.field by default

  # Build integration stack of nearest predictors
  stk_int <- INLA::inla.stack(
    data    = ll_resp,
    A       = list(1, idiag),
    tag     = paste0('stk_int_',as.character(id)),
    effects = ll_effects
  )
  return(stk_int)
}

#' Create a projection stack
#'
#' @param stk_resp A inla stack object.
#' @param model A prepared model object.
#' @param mesh The background projection mesh.
#' @param mesh.area The area calculate for the mesh.
#' @param type Name to use.
#' @param background A [`sf`] formatted background layer.
#' @param spde An spde field if specified.
#' @param res Approximate resolution to the projection grid (default: \code{NULL}).
#' @param settings A settings object.
#' @param joint Whether more than 2 likelihoods are estimated.
#'
#' @keywords utils
#'
#' @noRd
#'
#' @keywords internal
inla_make_projection_stack <- function(stk_resp, model, mesh, mesh.area, type, background,
                                       res = NULL, spde = NULL, settings = NULL,joint = FALSE){
  # Security checks
  assertthat::assert_that(
    inherits(stk_resp, 'inla.data.stack'),
    inherits(mesh,'inla.mesh'),
    is.numeric(mesh.area),
    is.character(type),
    is.null(spde)  || 'spatial.field1' %in% names(spde),
    is.logical(joint)
  )
  cov        <- model$predictors
  cov_object <- model$predictors_object
  pred.names <- model$predictors_names
  background <- model$background
  offset     <- model$offset

  # New formulation of the projection matrix
  bdry <- mesh$loc[mesh$segm$int$idx[,2],]  # get the boundary of the region

  # Try to make up a projection grid size if not defined
  if(is.null(res)){
    res <- (diff(range(cov$x)) / 100)
  }

  # Approximate dimension of the projector matrix
  Nxy <- round( c(diff(range(bdry[,1])), diff(range(bdry[,2]))) / res )

  # Make a mesh projection
  suppressWarnings(
    projgrid <- INLA::inla.mesh.projector(mesh,
                                          xlim = range(bdry[,1]),
                                          ylim = range(bdry[,2]),
                                          dims = Nxy)
  )

  if(!inherits(background, "sf")) background <- sf::st_as_sf( background )

  # Buffer the region to be sure
  suppressMessages( suppressWarnings( background.g <- sf::st_buffer( background, dist = 0) ) )
  # # Get and append coordinates from each polygon
  # background.bdry <- unique(
  #   do.call('rbind', lapply(background.g@polygons[[1]]@Polygons, function(x) return(x@coords) ) )
  # )
  # background.bdry <- background.g@polygons[[1]]@Polygons[[1]]@coords # use the original polygon boundaries to avoid discrepancies in the 'true' points
  # cellsIn <- splancs::inout(projgrid$lattice$loc,
  #                           cbind(background.bdry[,1], background.bdry[,2]))

  # Get only those points from the projection grid that are on the background
  projpoints <- projgrid$lattice$loc  |> as.data.frame() |> sf::st_as_sf(coords = c(1,2),crs = sf::st_crs(background))
  assertthat::assert_that(inherits(projpoints, "sf"), nrow(projpoints)>0)

  suppressMessages(
    suppressWarnings(
      # cellsIn <- !is.na(sp::over(x = sp::SpatialPoints(projgrid$lattice$loc,
      #                                    proj4string = methods::as(background.g,'Spatial')@proj4string),
      #                            y = background.g))
      check <- sf::st_intersects(projpoints, background.g)
    )
  )
  assertthat::assert_that(is.list(check))
  # Get the intersecting cells/points
  cellsIn <- which( sapply(check, function(z) length(z)>0) |> unlist() )

  # Check for multipolygon and align grid if necessary
  # if(inherits(cellsIn,'matrix')){
  #   cellsIn <- which(apply(cellsIn,1,function(x) any(x == TRUE)))
  # } else { cellsIn <- which(cellsIn) }
  # assertthat::assert_that(length(cellsIn)>0)

  # Get prediction coordinates
  predcoords <- projgrid$lattice$loc[cellsIn,]
  colnames(predcoords) <- c('x','y')

  # Security check
  assertthat::assert_that(length(cellsIn) == nrow(predcoords))

  # get the points on the grid within the boundary
  # predcoords <- projgrid$lattice$loc[which(cellsIn),]
  Apred <- projgrid$proj$A[cellsIn, ]

  # For all coordinates get nearest value
  nearest_cov <- get_rastervalue(coords = predcoords,
                                 env = cov_object$get_data(df = FALSE),
                                 rm.na = FALSE)

  # Set target variables to bias_value for prediction if specified
  if(!is.Waiver(settings$get('bias_variable'))){
    for(i in 1:length(settings$get('bias_variable'))){
      if(settings$get('bias_variable')[i] %notin% names(nearest_cov)) next()
      nearest_cov[!is.na(nearest_cov[[settings$get('bias_variable')[i]]]), settings$get('bias_variable')[i]] <- settings$get('bias_value')[i]
    }
  }

  # Extract covariates for points
  if(!is.Waiver(offset)) {
    ofs <- get_ngbvalue(coords = predcoords,
                        env = offset,
                        field_space = c('x','y'))
    nearest_cov[["spatial_offset"]] <- ofs[["spatial_offset"]]
  }

  # Get from supplied stack this information
  nearest_cov[,'Intercept'] <- 1


  # Empty lists
  ll_pred <- list()
  ll_effects <- list()

  # Define INLA stack dependent on what should go in it
  if(joint){
    # Joint stack
    # Set to NA to predict for fitted area
    ll_pred[[ 'observed' ]] <- cbind(rep(NA, nrow(nearest_cov)), rep(NA, nrow(nearest_cov)) )
    ll_pred[['e']] <- rep(0, nrow(nearest_cov))
    ll_pred[['Ntrials']] <- rep(1, nrow(nearest_cov))

    # Note, order adding this is important apparently...
    # ll_effects[['Intercept']] <- rep(1, nrow(nearest_cov))
    ll_effects[['predictors']] <- nearest_cov
    ll_effects[['spatial.field1']] <- list(spatial.field1 = seq(1,mesh$n)) # Changed to spatial.field. intercept already included in neatest_cov
    # if(!is.null(spde)) ll_effects[['spatial.field']] <- c(ll_effects[['spatial.field']], spde)

    # Set A
    A = list(1, Apred)

  } else if("Ntrials" %in% stk_resp$data$names) {
    # Single stack and binomial

    # Set to NA to predict for fitted area
    ll_pred[[ 'observed' ]] <- cbind( rep(NA, nrow(nearest_cov)) )
    ll_pred[['Ntrials']] <- rep(1, nrow(nearest_cov))

    # Note, order adding this is important apparently...
    # ll_effects[['Intercept']] <- rep(1, nrow(nearest_cov))
    ll_effects[['predictors']] <- nearest_cov
    ll_effects[['spatial.field1']]  <- list(spatial.field1 = seq(1,mesh$n))
    # if(!is.null(spde)) ll_effects[['spatial.field']] <- c(ll_effects[['spatial.field']], spde)

    # Set A
    A = list(1, Apred)
  } else {
    # --- #
    ll_pred[[ 'observed' ]] <- cbind( rep(NA, nrow(nearest_cov)) ) # Set to NA to predict for fitted area
    ll_pred[['e']] <- rep(0, nrow(nearest_cov))

    # Note, order adding this is important apparently...
    # ll_effects[['Intercept']] <- rep(1, nrow(nearest_cov))
    ll_effects[['predictors']] <- nearest_cov
    ll_effects[['spatial.field1']] <- list(spatial.field1 = seq(1, mesh$n) )
    # if(!is.null(spde)) ll_effects[['spatial.field']] <- c(ll_effects[['spatial.field']], spde)

    # Set A
    A = list(1, Apred)
  }

  # Build stack
  stk_proj <-
    INLA::inla.stack(
      data =  ll_pred,              # Response
      A = A,                        # Predictor projection matrix
      effects = ll_effects,         # Effects matrix
      tag = paste0('stk_pred')      # New description tag
    )

  # --- #
  # Return a list with the projection grid
  return(list(
    stk_proj = stk_proj,
    predcoords = predcoords,
    cellsIn = cellsIn
  )
  )
}

#' Prediction coordinates for INLA
#'
#' @param mesh A \code{"INLA::inla.mesh"} object.
#' @param background A [sf] object containing the background region.
#' @param cov A [data.frame] or [matrix] with the covariates for the modelling.
#' @param proj_stepsize A numeric indication on the prediction stepsize to be used.
#' @param spatial A [logical] flag whether a spatialpoints [data.frame] should be returned.
#'
#' @keywords utils
#'
#' @noRd
#'
#' @keywords internal
inla_predpoints <- function( mesh, background, cov, proj_stepsize = NULL, spatial = TRUE){
  assertthat::assert_that(
    inherits(mesh,'inla.mesh'),
    inherits(background, 'sf'),
    is.data.frame(cov) || is.matrix(cov) || is.Raster(cov),
    is.null(proj_stepsize) || is.numeric(proj_stepsize),
    is.logical(spatial)
  )

  # Get the boundary region size
  bdry <- mesh$loc[mesh$segm$int$idx[,2],]

  # Approximate dimension of the projector matrix
  if(is.null(proj_stepsize)){
    proj_stepsize <- (diff(range(bdry[,1])) / 100)
  }
  # Calculate the approximate cell size
  Nxy <- round( c(diff(range(bdry[,1])),
                  diff(range(bdry[,2]))) / proj_stepsize )

  # Make a INLA projection grid
  projgrid <- INLA::inla.mesh.projector(mesh,
                                        xlim = range(bdry[,1]),
                                        ylim = range(bdry[,2]),
                                        dims = Nxy)
  # Convert background to buffered land
  suppressMessages(
    suppressWarnings(
      background.g <- sf::st_buffer(methods::as(background, 'Spatial') |> sf::st_as_sf(),
                                    dist = 0) |> methods::as("Spatial")
    )
  )
  suppressWarnings(
    cellsIn <- !is.na(sp::over(x = sp::SpatialPoints(projgrid$lattice$loc,
                                                     proj4string = methods::as(background.g,'Spatial')@proj4string),
                               y = background.g))
  )
  # Get the cells that are in
  if(inherits(cellsIn,'matrix')){
    cellsIn <- which(apply(cellsIn,1,function(x) any(x == TRUE)))
  } else { cellsIn <- which(cellsIn) }
  assertthat::assert_that(length(cellsIn)>0)

  # Get prediction coordinates
  predcoords <- projgrid$lattice$loc[cellsIn,]
  colnames(predcoords) <- c('x','y')

  # Get covariates
  if(is.Raster(cov)){
    preds <- get_rastervalue(coords = predcoords,
                             env = cov,
                             rm.na = FALSE)
  } else {
    preds <- get_ngbvalue(coords = predcoords,
                          env = cov,
                          longlat = terra::is.lonlat(background),
                          field_space = c('x','y'))
  }

  if(spatial){
    # Convert predictors to SpatialPixelsDataFrame as required for inlabru
    # preds <- sp::SpatialPointsDataFrame(coords = model$predictors[,c('x', 'y')],
    #                                     data = model$predictors[, which(names(model$predictors) %in% fit_bru$names.fixed)],
    #                                     proj4string = self$get_data('mesh')$crs
    # )
    preds <- sp::SpatialPointsDataFrame(coords = predcoords,
                                        data = preds,
                                        proj4string = sp::CRS(SRS_string = mesh$crs)
    )
    # Remove missing data
    preds <- subset(preds, stats::complete.cases(preds@data))
    preds <- methods::as(preds, 'SpatialPixelsDataFrame')
  } else {
    preds <- subset(preds, stats::complete.cases(preds))
  }

  return(preds)
}

#' Tidy up summary information from a INLA model
#' TODO: Lot more to add here, including options on what to extract
#' @param m A trained INLA model object.
#' @param what A [`character`].
#' @param ... Other options to based on.
#'
#' @keywords utils
#'
#' @noRd
#'
#' @keywords internal
tidy_inla_summary <- function(m, what = 'fixed',...){
  assertthat::assert_that(
    inherits(m,'inla'),
    is.character(what),
    length(what)==1,
    what %in% c('fixed','fitted','random','spde2')
  )

  w1 <- grep('summary',names(m),value = TRUE) # Grep summary objects
  w2 <- grep(what, w1,value = TRUE) # Grep specific summary
  assertthat::assert_that(length(w2)==1)

  # Format the output
  o <- m[[w2]]  |>
    tibble::rownames_to_column('variable') |>
    tibble::as_tibble()
  if(what == "fixed"){
    names(o) <- c("variable", "mean", "sd", "q05", "q50", "q95", "mode", "kld")
  }
  assertthat::assert_that(nrow(o)>0)
  return( o )
}

#' Plot marginal distributions of effects or hyperparameters from INLA model
#'
#' @param A INLA model.
#' @param what Either \code{'fixed'} or \code{'hyper'}.
#'
#' @keywords utils
#'
#' @noRd
#'
#' @keywords internal
plot_inla_marginals = function(inla.model, what = 'fixed'){
  assertthat::assert_that(inherits(inla.model,'inla'),
                          is.character(what),
                          what %in% c('fixed','hyper'))
  par.ori <- graphics::par(no.readonly = TRUE)
  graphics::par(mfrow = c(4,3), mar = c(3,3,1,0.3), mgp = c(2,1,0))
  if(what == 'fixed'){
    varnames <- names(inla.model$marginals.fixed)
    for(i in 1: length(varnames)){
      var.mar <- data.frame(inla.model$marginals.fixed[i])
      plot(x = var.mar[,1], y=var.mar[, 2], type="l",
           xlab=paste(names(var.mar)[1]), ylab=paste(names(var.mar)[2]))
      graphics::abline(v=0, col="red")
    }
  } else {
    varnames <- names(inla.model$marginals.hyperpar)
    for(i in 1: length(varnames)){
      var.mar <- data.frame(inla.model$marginals.hyperpar[i])
      plot(x = var.mar[,1], y=var.mar[, 2], type="l",
           xlab=paste(names(var.mar)[1]), ylab=paste(names(var.mar)[2]))
    }
  }
  graphics::par(par.ori)
}

#' Additional INLA priors not already available.
#'
#' @param prior Which prior to pick as [`character`].
#'
#' @source https://becarioprecario.bitbucket.io/inla-gitbook/ch-priors.html#sec:newpriors
#'
#' @keywords utils
#'
#' @noRd
#'
#' @keywords internal
manual_inla_priors <- function(prior){

  # Uniform prior on the standard deviation
  UN.prior = "expression:
  log_dens = 0 - log(2) - theta / 2;
  return(log_dens);"

  # Half-Normal
  # a normal with zero mean and precision truncated at 0
  HN.prior = "expression:
  tau0 = 0.001;
  sigma = exp(-theta/2);
  log_dens = log(2) - 0.5 * log(2 * pi) + 0.5 * log(tau0);
  log_dens = log_dens - 0.5 * tau0 * sigma^2;
  log_dens = log_dens - log(2) - theta / 2;
  return(log_dens);
  "

  # Half-Cauchy prior
  # Here, we have set the scale parameter lambda to 25 following A Gelman (2006).
  HC.prior = "expression:
              sigma = exp(-theta/2);
              gamma = 25;
              log_dens = log(2) - log(pi) - log(gamma);
              log_dens = log_dens - log(1 + (sigma / gamma)^2);
              log_dens = log_dens - log(2) - theta / 2;
              return(log_dens);"

  switch (prior,
          'halfcauchy' = return(HC.prior),
          'uniform'  = return(UN.prior),
          'halfnormal' = return(HN.prior)
  )
}

#' Backward variable selection using INLA
#'
#' @description Best model is assessed through their within-sample predictive
#' accuracy via conditional predictive ordinate (CPO) Ideally this procedure is
#' replaced by a proper regularizing prior at some point...
#'
#' @param form A supplied [`formula`] object.
#' @param stack_data_resp A list containing inla stack data.
#' @param stk_inference An inla.data.stack object.
#' @param fam A [`character`] indicating the distribution to be fitted.
#' @param cf List of link functions to be used.
#' @param li Internal indication for the link function (Default: \code{1}).
#' @param response The response variable. If not specified, extract from formula
#' (default: \code{NULL}).
#' @param keep A [`vector`] of variables that are to be removed from model
#' iterations (default: \code{NULL}).
#'
#' @keywords utils
#'
#' @noRd
#'
#' @keywords internal
inla.backstep <- function(master_form,
                          stack_data_resp, stk_inference,fam, cf, li = 1,
                          response = NULL, keep = NULL
){
  assertthat::assert_that(is.formula(master_form),
                          is.list(stack_data_resp),inherits(stk_inference,'inla.data.stack'),
                          is.character(fam), is.list(cf), !missing(li),
                          is.null(response) || response %in% all.vars(master_form),
                          is.null(keep) || is.character(keep) || is.vector(keep)
  )
  # Get response term
  if(is.null(response)) response <- all.vars(master_form)[1]
  # Formula terms
  te <- attr(stats::terms.formula(master_form),'term.label')
  te <- te[grep('Intercept',te,ignore.case = T,invert = T)] # remove intercept(s)
  # Remove variables that are never removed
  if(!is.null(keep)){
    te <- te[grep(pattern = paste0(keep,collapse = '|'),x = te, invert = TRUE, fixed = TRUE )]
    # Also remove keep from master_form as we won't use them below
    master_form <- stats::as.formula(paste0(response,' ~ ', paste0(te,collapse = " + ")," - 1"))
  }

  assertthat::assert_that(length(te)>0, !is.null(response), all(keep %notin% te ))
  # --- #
  # Iterate through unique combinations of variables backwards
  pb <- progress::progress_bar$new(total = length(te),format = "Backward eliminating variables... :spin [:elapsedfull]")
  test_form <- master_form
  not_found <- TRUE
  best_found <- NULL
  while(not_found) {
    pb$tick()
    # --- #
    # Base Model #
    fit <- try({INLA::inla(formula = test_form, # The specified formula
                           data = stack_data_resp,  # The data stack
                           E = INLA::inla.stack.data(stk_inference)$e, # Expectation (Eta) for Poisson model
                           Ntrials = INLA::inla.stack.data(stk_inference)$Ntrials,
                           family = fam,   # Family the data comes from
                           control.family = cf, # Control options
                           control.predictor=list(A = INLA::inla.stack.A(stk_inference),
                                                  link = li,
                                                  compute = FALSE),  # Compute for marginals of the predictors
                           control.compute = list(cpo = TRUE,dic = TRUE, waic = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                           INLA::control.inla(int.strategy = "eb"), # Empirical bayes for integration
                           num.threads = getOption('ibis.nthread'),
                           control.fixed = list(mean = 0),
                           verbose = FALSE # Verbose for variable selection
    )
    },silent = TRUE)
    if(inherits(fit, 'try-error')) {not_found <- FALSE;next()}

    o <- data.frame(form = deparse1(test_form),
                    converged = fit$ok,
                    waic = fit$waic$waic,
                    dic = fit$dic$dic,
                    # conditional predictive ordinate values
                    # The sum of the log CPO's and is an estimator for the log marginal likelihood.
                    # The factor -2 is included as that places the measure on the same scale as other commonly used information criteria such as the DIC or WAIC.
                    cpo = sum(log(fit$cpo$cpo)) * -2,
                    mean.deviance = fit$dic$mean.deviance )

    oo <- data.frame()

    te <- attr(stats::terms.formula(test_form),'term.label')
    te <- te[grep('Intercept',te,ignore.case = T,invert = T)] # remove intercept(s)

    # Now for each term in variable list
    for(vars in te){
      # New formula
      new_form <- stats::update.formula(test_form, paste0('. ~ . - ',vars ))

      fit <- try({INLA::inla(formula = new_form, # The specified formula
                             data = stack_data_resp,  # The data stack
                             E = INLA::inla.stack.data(stk_inference)$e, # Expectation (Eta) for Poisson model
                             Ntrials = INLA::inla.stack.data(stk_inference)$Ntrials,
                             family = fam,   # Family the data comes from
                             control.family = cf, # Control options
                             control.predictor=list(A = INLA::inla.stack.A(stk_inference),
                                                    link = li,
                                                    compute = FALSE),  # Compute for marginals of the predictors
                             control.compute = list(cpo = TRUE,dic = TRUE, waic = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                             INLA::control.inla(int.strategy = "eb"), # Empirical bayes for integration
                             num.threads = getOption('ibis.nthread'),
                             control.fixed = list(mean = 0),
                             verbose = FALSE # Verbose for variable selection
      )
      },silent = TRUE)
      if(inherits(fit,'try-error')) next()

      oo <- rbind(oo, data.frame(form = deparse1(new_form),
                                 converged = fit$ok,
                                 waic = fit$waic$waic,
                                 dic = fit$dic$dic,
                                 # conditional predictive ordinate values
                                 cpo = sum(log(fit$cpo$cpo)) * -2,
                                 mean.deviance = fit$dic$mean.deviance )
      )
      rm(fit)
    } # End of loop

    # Best model among competing models has the largest CPO. Ratio of CPO (LPML
    # really, e.g. the transformation above) is a surrogate for a Bayes Factor
    # In the code below we take the minimum since the CPO has been multiplied
    # with -2 to emulate comparison with AIC-like statistics

    if(!is.na(o$cpo) || nrow(oo) > 0) {
      # Now check whether any of the new models are 'better' than the full model
      # If yes, continue loop, if no stop
      if(o$cpo <= min(oo$cpo,na.rm = TRUE)){
        not_found <- FALSE
        best_found <- o
      } else {
        # Get best model
        test_form <- stats::as.formula(oo$form[which.min(oo$cpo)])
      }
      rm(o,oo)
    } else {
      # Check whether formula is empty, if yes, set to not_found to FALSE
      te <- attr(stats::terms.formula(test_form),'term.label')
      if(length(te)<=3){
        not_found <- FALSE
        best_found <- o
      }
    }

  } # End of While loop
  # Make sure to add kept variables back
  if(!is.null(keep)){
    best_found$form <- paste0(best_found$form," + ", paste0(keep,collapse = " + "))
  }
  return(best_found)
}
