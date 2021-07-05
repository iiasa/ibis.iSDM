#' @include bdproto-engine.R utils-inla.R bdproto-distributionmodel.R
NULL

#' Use INLA as engine
#'
#' References
#' https://becarioprecario.bitbucket.io/inla-gitbook/ch-spatial.html#ref-Simpsonetal:2016
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param optional_mesh A directly supplied [`INLA`] mesh (Default: NULL)
#' @param optional_projstk A directly supplied projection stack. Useful if projection stack is identical for multiple species (Default: NULL)
#' @param max.edge The largest allowed triangle edge length, must be in the same scale units as the coordinates
#' @param offset interpreted as a numeric factor relative to the approximate data diameter;
#' @param cutoff The minimum allowed distance between points on the mesh
#' @param proj_stepsize The stepsize in coordinate units between cells of the projection grid (Default: NULL)
#' @param barrier Should a barrier model be added to the model?
#' @param nonconvex.bdry Create a non-convex boundary hulls instead (Default: FALSE) TBD
#' @param nonconvex.convex Non-convex mimal extension radius for convex curvature TBD
#' @param nonconvex.concave Non-convex minimal extension radius for concave curvature TBD
#' @param nonconvex.res Computation resolution for nonconvex.hulls TBD
#' @param ... Other variables
#' @name engine_inla
NULL
#' @rdname engine_inla
#' @export
engine_inla <- function(x,
                        optional_mesh = NULL,
                        optional_projstk = NULL,
                        max.edge = c(1,5),
                        offset = c(1,1),
                        cutoff = 1,
                        proj_stepsize = NULL,
                        barrier = FALSE,
                        # nonconvex.bdry = FALSE,
                        # nonconvex.convex = -0.15,
                        # nonconvex.concave = -0.05,
                        # nonconvex.res = 40,
                        ...) {

  # Check whether INLA package is available
  check_package('INLA')
  if(!isNamespaceLoaded("INLA")) { attachNamespace("INLA");requireNamespace('INLA') }

  # TODO:
  # Find a better way to pass on parameters such as those related to the mesh size...
  # assert that arguments are valid
  assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                          inherits(x$background,'sf'),
                          inherits(optional_mesh,'inla.mesh') || is.null(optional_mesh),
                          is.list(optional_projstk) || is.null(optional_projstk),
                          is.vector(max.edge),
                          is.vector(offset) || is.numeric(offset),
                          is.numeric(cutoff),
                          is.null(proj_stepsize) || is.numeric(proj_stepsize)
                          )

  # Convert the study region
  region.poly <- as(sf::st_geometry(x$background), "Spatial")

  # Set the projection mesh
  if(inherits(optional_mesh,'inla.mesh')) {
    # Load a provided on
    mesh <- optional_mesh
    # Security check for projection and if not set, use the one from background
    if(is.null(mesh$crs))  mesh$crs <- sp::CRS( proj4string(region.poly) )
  } else {
    # Create a new mesh
    # Convert to boundary object for later
    suppressWarnings(
      bdry <- INLA::inla.sp2segment(
        sp = region.poly,
        join = TRUE,
        crs = INLA::inla.CRS(projargs = sp::proj4string(region.poly))
      )
    )
    bdry$loc <- INLA::inla.mesh.map(bdry$loc)

    # --- #
    # FIXME: Move this to the manual description above
    # Create the mesh
    # A good mesh needs to have triangles as regular as possible in size and shape: equilateral
    suppressWarnings(
      mesh <- INLA::inla.mesh.2d(
        # Boundary object
        boundary = bdry,
        # The largest allowed triangle edge length, must be in the same scale units as the coordinates
        # Lower bounds affect the density of triangles
        max.edge = max.edge,
        # The automatic extension distance.
        # If positive: same scale units.
        # If negative, interpreted as a factor relative to the approximate data diameter;
        #   i.e., a value of -0.10 will add a 10% of the data diameter as outer extension.
        offset = offset,
        # The minimum allowed distance between points,
        # it means that points at a closer distance than the supplied value are replaced by a single vertex.
        # it is critical when there are some points very close to each other,
        #   either for point locations or in the domain boundary.
        cutoff = cutoff,
        # Define the CRS
        crs = bdry$crs
      )
    )
  }

  # Get barrier from the region polygon
  # TODO: Add this in addition to spatial field below, possibly specify an option to calculate this
  if(barrier){
    mesh_bar <- mesh_barrier(mesh, region.poly)
  } else { mesh_bar <- new_waiver() }

  # Calculate area in km²
  ar <- suppressWarnings(
    mesh_area(mesh = mesh,region.poly = region.poly, variant = 'gpc2')
  )

  # Print a message in case there is already an engine object
  if(!is.Waiver(x$engine)) message('Replacing currently selected engine.')

  # Set engine in distribution object
  x$set_engine(
    bdproto(
      "INLA-Engine",
      Engine,
      name = "<INLA>",
      data = list(
        'mesh' = mesh,
        'mesh.area' = ar,
        'mesh.bar' = mesh_bar,
        'proj_stepsize' = proj_stepsize,
        'stk_pred' = optional_projstk
      ),
      # Generic plotting function for the mesh
      plot = function(self, assess = FALSE){
        if(assess){
          # For an INLA mesh assessment
          out <- INLA::inla.mesh.assessment(self$data$mesh,
                                      spatial.range = 3,
                                      alpha = 2,
                                      dims = c(300, 300))
          # Convert to raster stack
          out <- raster::stack(
              sp::SpatialPixelsDataFrame( sp::coordinates(out), data = as.data.frame(out),
                                      proj4string = self$data$mesh$crs )
            )

          plot(out[[c('sd','sd.dev','edge.len')]],
               col = c("#00204D","#00336F","#39486B","#575C6D","#707173","#8A8779","#A69D75","#C4B56C","#E4CF5B","#FFEA46")
               )
        } else {
          plot( self$data$mesh )
        }
      },
      # Spatial latent function
      # https://groups.google.com/g/r-inla-discussion-group/c/eqMhlbwChkQ/m/m0b0PuzL-PsJ
      # Default SPDE prior
      # It computes the approximate diameter of the mesh, multiplies by 0.2 to get a value for the prior median range, and then transforms it to log-kappa scale by the formula
      # log(sqrt(8*nu)/range) where nu is alpha-dim/2.
      calc_latent_spatial = function(self,type = 'spde', alpha = 2,
                                     priors = NULL,
                                     ...){
        # Catch prior objects
        if(is.null(priors) || is.Waiver(priors)) priors <- NULL

        # For calculating iCAR process
        if(type == 'iCAR'){
          stop('This is not implemented as of now.')
          # convert mesh to sf object
          ns <- mesh_as_sf(self$data$mesh)
          # Create adjacency matrix with queen's case
          nc.nb <- spdep::poly2nb(ns,queen = TRUE)
          #Convert the adjacency matrix into a file in the INLA format
          adjmat <- INLA::inla.graph2matrix(nc.nb)
          # Save the adjaceny matrix as output
          self$data$latentspatial <- adjmat
          self$data$s.index <- as.numeric(attr(nc.nb,'region.id'))
          # Security checks
          assertthat::assert_that(length( self$data$s.index ) == nrow(ns))
        }
        if(type=='spde'){
          # Check that everything is correctly specified
          if(!is.null(priors)) if('spde' %notin% priors$varnames() ) priors <- NULL

          # Use default spde
          if(is.null(priors) || is.Waiver(priors)){
            # Define PC Matern SPDE model and save
            self$data$latentspatial <- INLA::inla.spde2.matern(
              mesh = self$data$mesh,
              alpha = alpha
            )
          } else {
            # Get priors
            pr <- if(is.null(priors)) c(0.01, 0.05) else priors$get('spde','prior.range')
            ps <- if(is.null(priors)) c(10, 0.05) else priors$get('spde','prior.sigma')

            # Define PC Matern SPDE model and save
            self$data$latentspatial <- INLA::inla.spde2.pcmatern(
              mesh = self$data$mesh,
              alpha = alpha,
              # P(Range < 1°) = 0.001 and P(sigma > 0.5) = 0.05
              prior.range = pr, prior.sigma = ps
            )
          }
          # Make index for spatial field
          self$data$s.index <- INLA::inla.spde.make.index(name = "spatial.field",
                                                          n.spde = self$data$latentspatial$n.spde)
          # Security checks
          assertthat::assert_that(
            inherits(self$data$latentspatial,'inla.spde'),
            length(self$data$s.index$spatial.field) == self$data$mesh$n
          )
        }
        invisible()
      },
      # Get latent spatial equation bit
      get_equation_latent_spatial = function(self,spatial_model){
        if(spatial_model=='spde'){
          assertthat::assert_that(inherits(self$data$latentspatial, 'inla.spde'),
                                  msg = 'Latent spatial has not been calculated.')
          return(
            paste0('f(spatial.field, model = ',spatial_model,')')
          )
        } else if(spatial_model == 'iCAR'){
          assertthat::assert_that(inherits(self$data$latentspatial,'dgTMatrix'),
                                  msg = 'Neighborhood matrix has not been calculated.')
          return(
            paste0('f(','spatial.index',', model = "besag", graph = ','adjmat',')')
          )
        }
      },
      # Configure stack
      make_stack = function(self, model, id, intercept = TRUE, joint = FALSE) {
        assertthat::assert_that(
          is.list(model),
          is.character(id)
        )
        # Get Environment records
        env <- model$predictors

        # Include intercept in here
        # TODO: Note that this sets intercepts by type and not by dataset id
        if(intercept) {
          env$intercept <- 1 # Overall intercept
          env[[paste0('intercept',
                      ifelse(joint,paste0('_',
                                          make.names(tolower(model$name)),'_',
                                          model$type),''))]] <- 1 # Setting intercept to common type, thus sharing with similar types
        }
        # Set up projection matrix for the data
        suppressWarnings(
          mat_proj <- INLA::inla.spde.make.A(
            mesh = self$get_data('mesh'),
            loc = as.matrix(env[,c('x','y')])
          )
        )
        # Create INLA stack
        # The three main inla.stack() arguments are a vector list with the data (data),
        # a list of projector matrices (each related to one block effect,
        # A) and the list of effects (effects).

        # Response for inla stack
        ll_resp <- list()
        # Add the expected estimate and observed note
        # FIXME: Currently only two likelihoods are supported (binomial/poisson) with the NA order being the determining factor
        if(model$family == 'poisson') {
          if(joint) ll_resp[[ 'observed' ]] <- cbind(model$observations[['observed']], NA )
          if(!joint) ll_resp[[ 'observed' ]] <- cbind(model$observations[['observed']] )
          ll_resp[[ 'e' ]] <- model$expect
        }
        if(model$family == 'binomial') {
          if(joint) ll_resp[[ 'observed' ]] <- cbind(NA, model$observations[['observed']] )
          if(!joint) ll_resp[[ 'observed' ]] <- cbind( model$observations[['observed']] )
          ll_resp[[ 'Ntrials' ]] <- model$expect
        }

        # Effects matrix
        ll_effects <- list()
        # Note, order adding this is important apparently...
        ll_effects[['predictors']] <- env
        ll_effects[['spatial.field']] <- seq(1, self$get_data('mesh')$n)
        # ll_effects[['intercept']][[paste0('intercept',ifelse(joint,paste0('_',make.names(tolower(model$name)),'_',model$type),''))]]  <- seq(1, self$get_data('mesh')$n)

        # Add offset if specified
        if(!is.null(model$offset)){
         ll_effects[['predictors']] <- cbind( ll_effects[['predictors']],
                                              subset(model[['offset']],select = 3) # FIXME: If I want type-specific offsets, this would be the place
                                              )
        }

        # Check whether equation has spatial field and otherwise add
        # MJ 13/06: Spatial.field now set directly to effects
         # if( 'spde' %in% all.vars(model$equation) ){
         #   # Get Index Objects
         #   iset <- self$get_data('s.index')
         #   ll_effects[['spatial.field']] <- c(ll_effects[['spatial.field']], iset)
         # } else if ( 'adjmat' %in% all.vars(model$equation) ){
         #   iset <- self$get_data('s.index')
         #   ll_effects[['spatial.field']] <- c(ll_effects[['spatial.field']], data.frame(spatial.index = iset) )
         # }
        # Define A
        A <- list(1, mat_proj)

        # Define stack
        stk <- INLA::inla.stack(
          data     = ll_resp,
          A        = A,
          effects  = ll_effects,
          tag      = paste0('stk_',as.character(model$type),'_',id)
        )
        # Set the stack
        self$set_data(paste0('stk_',as.character(model$type),'_',id), stk)
        invisible()
      },
      # Setup computation function
      setup = function(self, model, ...){
        assertthat::assert_that(
          'background' %in% names(model),
          'biodiversity' %in% names(model),
          all( model$biodiversity[[1]]$predictors_names %in% model$predictors_names ),
          all( sapply(model$biodiversity, function(x) is.formula(x$equation)) ),
          length(model$biodiversity)>=1,
          msg = 'Some internal checks failed while setting up the model.'
        )

        # --- Prepare general inputs ---
        # Check whether spatial latent effects were added
        if( 'spde' %in% all.vars(model$biodiversity[[1]]$equation) ){
          # Get spatial index
          spde <- self$get_data('s.index')
        } else { spde <- NULL }

        # Check for existence of specified offset and use the full one in this case
        if(!is.Waiver(model$offset)) offset <- subset(model[['offset']],select = 3) else offset <- NULL

        # Projection stepsize
        if(is.null( self$get_data('proj_stepsize') )){
          # Set to stepsize equivalent of the resolution of the grid
          val <- max(diff(model[['predictors']]$x)) # TODO: Check that it works when dummy variable is used
          self$set_data('proj_stepsize', val )
          rm(val)
        }

        # Number of types to determine if a joint model is necessary
        nty <- length( unique( as.character(sapply(model$biodiversity, function(z) z$type)) ) )

        # Clean up previous data and integration stacks
        chk <- grep('stk_int|stk_poipo|stk_poipa|stk_polpo|stk_polpa|stk_pred', self$list_data())
        if(length(chk)>0) self$data[chk] <- NULL

        # Now for each dataset create a INLA stack
        for(id in 1:length(model$biodiversity) ){
          # Calculate observation stack INLA stack
          # Save stacks by id instead of type
          self$make_stack(model = model$biodiversity[[id]],
                          id = names(model$biodiversity)[id],
                          intercept = TRUE,
                          joint = ifelse(nty > 1, TRUE, FALSE)
                          )

          # Define mesh.area dependent on whether a single variable only is used or not
          if(model$biodiversity[[id]]$family == 'poisson'){
            # Only create on if not already existing
            chk <- grep('stk_int', self$list_data())
            if(length(chk)==0){
              # Make integration stack for given poisson model
              stk_int <- inla_make_integration_stack(
                mesh      = self$get_data('mesh'),
                mesh.area = self$get_data('mesh.area'),
                cov       = model$predictors,
                pred_names= model$predictors_names,
                bdry      = model$background,
                id        = names(model$biodiversity)[id],
                joint     = ifelse(nty > 1, TRUE, FALSE)
              )
              # Save integration stack
              self$set_data(paste0('stk_int_',names(model$biodiversity)[id]),stk_int)
            }
          }
        }

        # ------------------ #
        # Get all stacks defined so far and join them
        stk_inference <- lapply(
          self$list_data()[grep(paste(names(model$biodiversity),collapse = '|'), self$list_data())],
            function(x) self$get_data(x)
          )
        stk_inference <- do.call(INLA::inla.stack, stk_inference)

        # Make projection stack if not directly supplied
        if(is.null(self$data$stk_pred)){
          stk_pred <- inla_make_projection_stack(
            stk_resp   = stk_inference,
            cov        = model$predictors,
            pred.names = model$predictors_names,
            offset     = model$offset,
            mesh       = self$get_data('mesh'),
            mesh.area  = self$get_data('mesh.area'),
            background = model$background,
            res        = self$get_data('proj_stepsize'),
            type       = model$biodiversity[[id]]$type,
            spde       = spde,
            joint      = ifelse(nty > 1, TRUE, FALSE)
          )
          self$set_data('stk_pred', stk_pred)
        } else {
          # FIXME: Add some basic assertthat tests for when a prediction stack is directly supplied
          stk_pred <- self$get_data('stk_pred')
        }

        # Now join all stacks and save in full
        # Note: If integrated stack is included, E must be set to relative area (in mesh.area).
        self$set_data('stk_full',
                      INLA::inla.stack(stk_inference, stk_pred$stk_proj)
                      )
        invisible()
      },
      # Main INLA training function ----
      train = function(self, model, settings) {
        # Check that all inputs are there
        assertthat::assert_that(
          inherits(settings,'Settings'),
          is.list(model),length(model)>1,
          # Check that model id and setting id are identical
          settings$modelid == model$id,
          any(  (c('stk_full','stk_pred') %in% names(self$data)) ),
          inherits(self$get_data('stk_full'),'inla.data.stack')
        )

        # Get all datasets with id. This includes the data stacks and integration stacks
        stk_inference <- lapply(
          self$list_data()[grep(paste(names(model$biodiversity),collapse = '|'), self$list_data())],
                                function(x) self$get_data(x))
        stk_inference <- do.call(INLA::inla.stack, stk_inference)

        # Get full stack and projection grid
        stk_full <- self$get_data('stk_full')
        predcoords <- self$get_data('stk_pred')$predcoords

        # Get families
        fam <- unique( as.character( sapply(model$biodiversity, function(x) x$family) ) )
        # Define control family
        cf <- list()
        for(i in 1:length(fam)) cf[[i]] <- list(link = ifelse(fam[i] == 'poisson','log','cloglog' ))
        if(length(fam)==1 && fam == 'binomial') cf[[1]]$link <- 'logit'

        # Shared link?
        if(length(fam)==1) li <- 1 else li <- NULL

        if('spde' %in% all.vars(model$biodiversity[[1]]$equation) ){
            spde <- self$get_data('latentspatial')
            stack_data_resp <- INLA::inla.stack.data(stk_inference, spde = self$get_data('latentspatial'))
            stack_data_full <- INLA::inla.stack.data(stk_full, spde = self$get_data('latentspatial'))
        } else {
            # FIXME: Make sure this work for other types in the future
            spde <- self$get_data('latentspatial')
            stack_data_resp <- INLA::inla.stack.data(stk_inference)
            stack_data_full <- INLA::inla.stack.data(stk_full)
        }
        # ----------- #
        # Provided or default formula
        # TODO: Currently duplicate equations per dataset. Use only one here, but ideally support multiple?
        master_form <- model$biodiversity[[1]]$equation

        # Perform variable selection
        if( settings$get(what='varsel') ){
          message('Performing variable selection...')

          k <- NULL
          # Specify offsets and spde to be retained
          # FIXME: Also set priors and offsets here?
          if(is.Waiver(spde)) k <- NULL else k <- self$get_equation_latent_spatial('spde')

          # Use backward variable elimination
          vs <- inla.backstep(master_form = master_form,
                              stack_data_resp = stack_data_resp,
                              stk_inference = stk_inference,
                              fam = fam,
                              cf = cf,li = li,
                              response = 'observed',
                              keep = k
                              )
          master_form <- to_formula(vs$form)
        }

        # ------------------------------------------ #
        # Train the model on the response
        fit_resp <- INLA::inla(formula = master_form, # The specified formula
                          data  = stack_data_resp,  # The data stack
                          quantiles = c(0.05, 0.5, 0.95),
                          E = INLA::inla.stack.data(stk_inference)$e, # Expectation (Eta) for Poisson model
                          Ntrials = INLA::inla.stack.data(stk_inference)$Ntrials,
                          family = fam,   # Family the data comes from
                          control.family = cf, # Control options
                          control.predictor = list(A = INLA::inla.stack.A(stk_inference),
                                                   link = li, # Link to NULL for multiple likelihoods!
                                                   compute = TRUE),  # Compute for marginals of the predictors.
                          control.compute = list(cpo = FALSE, waic = FALSE, config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                          control.fixed = list(mean = 0),# prec = list( initial = log(0.000001), fixed = TRUE)), # Added to see whether this changes GMRFlib convergence issues
                          verbose = settings$get(what='verbose'), # To see the log of the model runs
                          INLA::control.inla(int.strategy = "eb"), # Empirical bayes for integration
                          #              strategy = 'simplified.laplace'
                          #              # https://groups.google.com/g/r-inla-discussion-group/c/hDboQsJ1Mls
                          # ), # To make it run faster...
                          num.threads = getOption('ibis.nthread')
        )

        # Predict spatially
        if(!settings$get(what='inference_only')){
          # Get thetas from initially fitted model as starting parameters
          thetas = fit_resp$mode$theta

          # Predict on full
          fit_pred <- try({INLA::inla(formula = master_form, # The specified formula
                                 data  = stack_data_full,  # The data stack
                                 quantiles = c(0.05, 0.5, 0.95),
                                 E = INLA::inla.stack.data(stk_full)$e,
                                 Ntrials = INLA::inla.stack.data(stk_full)$Ntrials,
                                 family= fam,   # Family the data comes from
                                 control.family = cf, # Control options
                                 control.predictor = list(A = INLA::inla.stack.A(stk_full),
                                                          link = li, # Link to NULL for multiple likelihoods!
                                                        compute = TRUE),  # Compute for marginals of the predictors.
                                 control.compute = list(cpo = TRUE, waic = TRUE, config = TRUE, openmp.strategy	= 'huge' ),
                                 # control.mode = list(theta = thetas, restart = FALSE), # To speed up use previous thetas
                                 verbose = settings$get(what='verbose'), # To see the log of the model runs
                                 control.results = list(return.marginals.random = FALSE,
                                                        return.marginals.predictor = FALSE), # Don't predict marginals to save speed
                                 # MJ: 15/6 -> Removed thetas as those cause SPDE convergence issues making the whole estimation slower
                                 # control.fixed = INLA::control.fixed(mean = 0),#, prec = list( initial = log(0.000001), fixed = TRUE)), # Added to see whether this changes GMRFlib convergence issues
                                 INLA::control.inla(int.strategy = "eb"), # Empirical bayes for integration
                                 #                    strategy = 'simplified.laplace'
                                 #                    # https://groups.google.com/g/r-inla-discussion-group/c/hDboQsJ1Mls
                                 # ),
                                 num.threads = getOption('ibis.nthread')
            )
          },silent = TRUE)
          if(class(fit_pred)=='try-error') stop('Model did not converge. Try to simplify structure and check priors!')
          # Create a spatial prediction
          index.pred <- INLA::inla.stack.index(stk_full, 'stk_pred')$data
          post <- fit_pred$summary.linear.predictor[index.pred, ]
          assertthat::assert_that(nrow(post)>0,
                                  nrow(post) == nrow(predcoords) ) # Check with cells in projection
          # TODO: Wrap Backtransform into a distribution dependent function
          if(length(fam)==1){
            if(fam == 'poisson') post[,c('mean','0.05quant','0.5quant','0.95quant','mode')] <- exp( post[,c('mean','0.05quant','0.5quant','0.95quant','mode')] )
            if(fam == 'binomial') post[,c('mean','0.05quant','0.5quant','0.95quant','mode')] <- logistic(post[,c('mean','0.05quant','0.5quant','0.95quant','mode')])
          } else{
            # Joint likelihood of Poisson log and binomial cloglog following Simpson et al.
            post[,c('mean','0.05quant','0.5quant','0.95quant','mode')] <- exp( post[,c('mean','0.05quant','0.5quant','0.95quant','mode')] )
          }
          post <- subset(post, select = c('mean','sd','0.05quant','0.5quant','0.95quant','mode','kld') )
          names(post) <- make.names(names(post))

          # Fill prediction
          prediction <- raster::stack(
            sp::SpatialPixelsDataFrame(
              points = predcoords,
              data = post,
              proj4string = CRS( self$get_data('mesh')$crs@projargs ) # x$engine$data$mesh$crs@projargs
            )
          )
          prediction <- raster::mask(prediction, model$background) # Mask with background
        } else {
          # No prediction to be conducted
          fit_pred <- NULL
          prediction <- NULL
        }

        # Compute end of computation time
        settings$set('end.time', Sys.time())

        # Definition of INLA Model object ----
        out <- bdproto(
              "INLA-Model",
              DistributionModel,
              id = new_id(),
              model = model,
              settings = settings,
              fits = list(
                "fit_best" = fit_resp,
                "fit_pred" = fit_pred,
                "fit_best_equation" = master_form,
                "mesh"     = self$get_data('mesh'),
                "spde"     = self$get_data('latentspatial'),
                "prediction" = prediction
                ),
              # Function to plot SPDE if existing
              plot_spatial = function(self, dim = c(300,300), kappa_cor = FALSE, dis = NULL,...){
                assertthat::assert_that(is.vector(dim), is.numeric(dis) || is.null(dis))
                if( length( self$fits$fit_best$size.spde2.blc ) == 1)
                {
                  # Get spatial projections from model
                  # FIXME: Potentially make the plotting of this more flexible
                  gproj <- INLA::inla.mesh.projector(self$get_data('mesh'),  dims = dim)
                  g.mean <- INLA::inla.mesh.project(gproj,
                                                    self$get_data('fit_pred')$summary.random$spatial.field$mean)
                  g.sd <- INLA::inla.mesh.project(gproj, self$get_data('fit_pred')$summary.random$spatial.field$sd)

                  # Convert to rasters
                  g.mean <- t(g.mean)
                  g.mean <- g.mean[rev(1:length(g.mean[,1])),]
                  r.m <- raster::raster(g.mean,
                                      xmn = range(gproj$x)[1], xmx = range(gproj$x)[2],
                                      ymn = range(gproj$y)[1], ymx = range(gproj$y)[2],
                                      crs = self$get_data('mesh')$crs
                  )
                  g.sd  <- t(g.sd)
                  g.sd <- g.sd[rev(1:length(g.sd[,1])),]
                  r.sd <- raster::raster(g.sd,
                                      xmn = range(gproj$x)[1], xmx = range(gproj$x)[2],
                                      ymn = range(gproj$y)[1], ymx = range(gproj$y)[2],
                                      crs = self$get_data('mesh')$crs
                  )

                  spatial_field <- raster::stack(r.m, r.sd);names(spatial_field) <- c('SPDE_mean','SPDE_sd')
                  # Mask with prediction
                  spatial_field <- raster::resample(spatial_field, self$get_data('prediction')[[1]])
                  spatial_field <- raster::mask(spatial_field, self$get_data('prediction')[[1]])

                  # -- #
                  if(kappa_cor){
                    # Also build correlation fun
                    # Get SPDE results
                    spde_results <- INLA::inla.spde2.result(
                      inla = self$get_data('fit_pred'),
                      name = 'spatial.field',
                      spde = self$get_data('spde'),
                      do.transfer = TRUE)

                    # Large kappa (inverse range) equals a quick parameter change in space.
                    # Small kappa parameter have much longer, slower gradients.
                    Kappa <- lapply(spde_results,  function(j) INLA::inla.emarginal(function(x) x, j$marginals.kappa[[1]] ))

                    d.vec <- seq(0, MaxRange, length = Resolution)
                    # Distance vector. Maximum distance being defined by half of extent
                    if(is.null(dis)) dis <- abs((xmin(self$get_data('prediction')) - xmax(self$get_data('prediction')) ) / 4)  # Take a quarter of the max distance

                    # Modified Bessel function to get coorelation strength
                    dis.cor <- data.frame(distance = seq(0, dis, length = dim[1]))
                    dis.cor$cor <- as.numeric((Kappa * dis.cor$distance) * besselK(Kappa * dis.cor$distance, 1))
                    dis.cor$cor[1] <- 1
                  # ---
                  # Build plot
                  cols <- c("#00204DFF","#00336FFF","#39486BFF","#575C6DFF","#707173FF","#8A8779FF","#A69D75FF","#C4B56CFF","#E4CF5BFF","#FFEA46FF")

                  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))

                  plot(dis.cor$cor ~ dis.cor$distance, type = 'l', lwd = 3,
                       xlab = 'Distance', ylab = 'Correlation', main = paste0('Kappa: ', round(Kappa,2) ) )
                  plot(spatial_field[[1]],col = cols, main = 'mean spatial effect')
                  plot(spatial_field[[2]], main = 'sd spatial effect')
                  } else {
                  # Just plot the SPDE
                    cols <- c("#00204DFF","#00336FFF","#39486BFF","#575C6DFF","#707173FF","#8A8779FF","#A69D75FF","#C4B56CFF","#E4CF5BFF","#FFEA46FF")
                    par(mfrow=c(1,2))
                    plot(spatial_field[[1]],col = cols, main = 'mean spatial effect')
                    plot(spatial_field[[2]], main = 'sd spatial effect')
                    # And return
                    return(spatial_field)
                }


                } else {
                  message(text_red('No spatial covariance in model specified.'))
                }
              }
        )
        return(out)
      }
      ))
}
