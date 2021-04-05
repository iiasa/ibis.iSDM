#' @include bdproto-engine.R utils-inla.R bdproto-distributionmodel.R
NULL

#' Use INLA as engine
#'
#' References
#' https://becarioprecario.bitbucket.io/inla-gitbook/ch-spatial.html#ref-Simpsonetal:2016
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param optional_mesh A directly supplied [`INLA`] mesh (Default: NULL)
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
    # FIXME: Move this to the man description above
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
        'proj_stepsize' = proj_stepsize
      ),
      # Generic plotting function for the mesh
      plot = function(self){
        plot( self$get_data('mesh') )
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
          # Get prior
          pr <- if(is.null(priors)) c(0.01, 0.05) else priors$get('spde','prior.range')
          ps <- if(is.null(priors)) c(10, 0.05) else priors$get('spde','prior.sigma')
          # Use default spde
          if(is.null(priors)){
            # Define PC Matern SPDE model and save
            self$data$latentspatial <- INLA::inla.spde2.matern(
              mesh = self$data$mesh,
              alpha = alpha
            )
          } else {
            # Define PC Matern SPDE model and save
            self$data$latentspatial <- INLA::inla.spde2.pcmatern(
              mesh = self$data$mesh,
              alpha = alpha,
              # P(Range < 1°) = 0.001 and P(sigma > 0.5) = 0.05
              prior.range = pr,prior.sigma = ps
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
      # Calculate stack for presence only records
      calc_stack_poipo = function(self, model, intercept = TRUE) {
        assertthat::assert_that(
          is.list(model),
          'poipo' %in% names(model$data),'poipo_response' %in% names(model$data),
          'poipo_values' %in% names(model$data), 'poipo_expect' %in% names(model$data)
        )
        # Get Environment records
        env <- model$data$poipo_values

        # Include intercept in here
        if(intercept) env$intercept <- 1 # FIXME: Rename intercepts when using varying likelihoods

        # Set up projection matrix for the data
        mat_proj <- INLA::inla.spde.make.A(
          mesh = self$get_data('mesh'),
          loc = as.matrix(env[,c('x','y')])
        )
        # Create INLA stack
        # The three main inla.stack() arguments are a vector list with the data (data),
        # a list of projector matrices (each related to one block effect,
        # A) and the list of effects (effects).

        # Response for inla stack
        resp <- model$data$poipo_response
        ll_resp <- list()
        ll_resp[[ resp ]] <- cbind( rep(1, nrow(env)) ) #FIXME Set to NA if multiple likelihoods
        # Add the expect
        ll_resp[[ 'e' ]] <- model$data$poipo_expect

        # Effects matrix
        ll_effects <- list()
        pred_names <- model$predictors_names
        # Note, order adding this is important apparently...
        ll_effects[['predictors']] <- env[,pred_names]
        ll_effects[['intercept']] <- list(intercept = seq(1, self$get_data('mesh')$n) )

        # Add offset if specified
        if('offset_poipo' %in% names(model)){
         ll_effects[['predictors']] <- cbind( ll_effects[['predictors']],
                                              subset(model[['offset_poipo']],select = 3)
                                              )
        }

        # Check whether equation has spatial field
         if( 'spde' %in% all.vars(model$equation$poipo) ){
           # Get Index Objects
           iset <- self$get_data('s.index')
           ll_effects[['intercept']] <- c(ll_effects[['intercept']], iset)
         } else if ( 'adjmat' %in% all.vars(model$equation$poipo) ){
           iset <- self$get_data('s.index')
           ll_effects[['intercept']] <- c(ll_effects[['intercept']], data.frame(spatial.index = iset) )
         }
        # Define A
        A <- list(1, mat_proj)

        # Define stack
        stk <- INLA::inla.stack(
          data     = ll_resp,
          A        = A,
          effects  = ll_effects,
          tag      = 'stk_poipo'
        )
        # Set the stack
        self$set_data('stk_poipo',stk)
        invisible()
      },
      # Setup computation function
      setup = function(self, model, ...){
        # TODO: Some assert calls
        # TODO: Potentially add a setup log later for this function

        # Calculate observation stack INLA stack
        self$calc_stack_poipo(model, intercept = TRUE)
        stk_poipo <- self$get_data('stk_poipo')

        # Make integration stack
        stk_int <- inla_make_integration_stack(
          mesh      = self$get_data('mesh'),
          mesh.area = self$get_data('mesh.area'),
          cov       = model$predictors,
          pred_names= model$predictors_names,
          bdry      = model$background,
          resp      = model$data$poipo_response
        )
        self$set_data('stk_int',stk_int)

        # ------------------ #
        if( 'spde' %in% all.vars(model$equation$poipo) ){
          # Get spatial index
          spde <- self$get_data('s.index')
        } else { spde <- NULL}

        # Check for existence of specified offset and use the full one in this case
        if('offset' %in% names(model)) offset <- subset(model[['offset']],select = 3) else offset <- NULL

        # Projection stepsize
        if(is.null( self$get_data('proj_stepsize') )){
          # Set to stepsize equivalent of the resolution of the grid
          val <- max(diff(model[['predictors']]$x)) # TODO: Check that it works when dummy variable is used
          self$set_data('proj_stepsize', val )
          rm(val)
        }

        # Make projection stack
        # stk_pred_poipo <- inla_make_prediction_stack(
        #   stk_resp = stk_poipo,
        #   cov       = model$predictors,
        #   pred.names= model$predictors_names,
        #   offset    = offset,
        #   mesh      = self$get_data('mesh'),
        #   mesh.area = self$get_data('mesh.area'),
        #   type = 'poipo',
        #   spde = spde
        # )
        stk_pred_poipo <- inla_make_projection_stack(
          stk_resp   = stk_poipo,
          cov        = model$predictors,
          pred.names = model$predictors_names,
          offset     = offset,
          mesh       = self$get_data('mesh'),
          mesh.area  = self$get_data('mesh.area'),
          background = model$background,
          res        = self$get_data('proj_stepsize'),
          type       = 'poipo',
          spde       = spde
        )
        self$set_data('stk_pred_poipo',stk_pred_poipo)

        # Now join all stacks and save in full
        # Note: If integrated stack is included, E must be set to relative area (in mesh.area).
        self$set_data('stk_full', INLA::inla.stack(stk_poipo, stk_int,
                                                   stk_pred_poipo$stk_proj
                                                   )
                      )
        invisible()
      },
      # Main INLA training function ----
      train = function(self, model, varsel = TRUE, inference_only = FALSE, verbose = FALSE,...) {
        # Check that all inputs are there
        assertthat::assert_that(
          is.list(model),length(model)>1,
          any(  (c('stk_poipo','stk_polpo','stk_polpa','stk_poipa') %in% names(self$data)) ),
          'stk_int' %in% names(self$data), inherits(self$get_data('stk_int'),'inla.data.stack'),
          'stk_full' %in% names(self$data), inherits(self$get_data('stk_full'),'inla.data.stack')
        )
        # Get the datasets
        stk_poipo <- self$get_data('stk_poipo')
        stk_int <- self$get_data('stk_int')
        # Get full stack and projection grid
        stk_full <- self$get_data('stk_full')
        predcoords <- self$get_data('stk_pred_poipo')$predcoords
        # Make joint stack of data and integrations
        # TODO: Needs to be done for all other types
        stk_var <- INLA::inla.stack(stk_poipo, stk_int)

        if('spde' %in% all.vars(model$equation$poipo) ){
            spde <- self$get_data('latentspatial')
            stack_data_resp <- INLA::inla.stack.data(stk_var, spde = self$get_data('latentspatial'))
            stack_data_full <- INLA::inla.stack.data(stk_full, spde = self$get_data('latentspatial'))
          } else {
            # FIXME: Make sure this work for other types in the future
            stack_data_resp <- INLA::inla.stack.data(stk_var)
            stack_data_full <- INLA::inla.stack.data(stk_full)
          }
        # ----------- #
        # Provided or default formula
        master_form <- model$equation$poipo

        # Perform variable selection
        if(varsel){
          message('Performing variable selection...')
          #
          if('spde' %in% all.vars(model$equation$poipo) ) speq <- self$get_equation_latent_spatial('spde') else speq <- NULL
          # Get formula list
          lf <- formula_combinations(varnames = model$predictors_names,
                               response = all.vars(master_form)[1],
                               InterceptOnly = FALSE,
                               spde_term = speq,
                               type = 'forward'
                               )

          # Set progress bar
          pb <- progress::progress_bar$new(total = length(lf))

          results <- data.frame(stringsAsFactors = FALSE)
          # Now loop through the formulas
          for(k in 1:length(lf)){
            pb$tick()

            # Train the model on the response
            fit <- INLA::inla(formula = as.formula(lf[k]), # The specified formula
                                   data = stack_data_resp,  # The data stack
                                   E = INLA::inla.stack.data(stk_var)$e, # Expectation (Eta) for Poisson model
                                   family= 'poisson',   # Family the data comes from
                                   control.family = list(link = "log"), # Control options
                                   control.predictor=list(A = INLA::inla.stack.A(stk_var),
                                                          link = NULL, compute = TRUE),  # Compute for marginals of the predictors
                                   control.compute = list(cpo = TRUE,dic = TRUE, waic = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                                   INLA::control.inla(#int.strategy = "eb", # Empirical bayes for integration
                                                strategy = 'simplified.laplace', huge = TRUE), # To make it run faster...
                                   num.threads = parallel::detectCores() - 1
            )
            # Add results
            results <- rbind(results,
                             data.frame(form = lf[k],
                                        converged = fit$ok,
                                        waic = fit$waic$waic,
                                        dic = fit$dic$dic,
                                        # conditional predictive ordinate values
                                        cpo = sum(log(fit$cpo$cpo)) * -2,
                                        mean.deviance = fit$dic$mean.deviance ) )
            rm(fit)
          }
          rm(pb)
          # Determine best model by DIC
          # Alternative: Negative sum of the log CPO
          #   - sum(log(m$cpo$cpo), na.rm = na.rm)
          results <- results[order(results$dic,decreasing = FALSE),]
          results$diff_dic <- c(0, diff(results$dic) )
          # Set new master form
          master_form <- as.formula(
            as.character(results$form[results$diff_dic==0])
          )
        }

        # Train the model on the response
        fit_resp <- INLA::inla(formula = master_form, # The specified formula
                          data  = stack_data_resp,  # The data stack
                          quantiles = c(0.05, 0.5, 0.95),
                          E = INLA::inla.stack.data(stk_var)$e, # Expectation (Eta) for Poisson model
                          family = 'poisson',   # Family the data comes from
                          control.family = list(link = "log"), # Control options
                          control.predictor=list(A = INLA::inla.stack.A(stk_var),link = 1, compute = TRUE),  # Compute for marginals of the predictors. Link to NULL for multiple likelihoods!
                          control.compute = list(cpo = TRUE, waic = TRUE, config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                          # control.fixed = INLA::control.fixed(prec = list( initial = log(0.000001), fixed = TRUE)), # Added to see whether this changes GMRFlib convergence issues
                          verbose = verbose, # To see the log of the model runs
                          INLA::control.inla(int.strategy = "eb", # Empirical bayes for integration
                                       strategy = 'simplified.laplace'
                                       # https://groups.google.com/g/r-inla-discussion-group/c/hDboQsJ1Mls
                          ), # To make it run faster...
                          num.threads = parallel::detectCores()-1
        )

        # Predict spatially
        if(!inference_only){
          # Get theta from initiall fitted model as starting parameters
          thetas = fit_resp$internal.summary.hyperpar$mean

          # Predict on full
          fit_pred <- INLA::inla(formula = master_form, # The specified formula
                                 data  = stack_data_full,  # The data stack
                                 quantiles = c(0.05, 0.5, 0.95),
                                 E = INLA::inla.stack.data(stk_full)$e,
                                 family= 'poisson',   # Family the data comes from
                                 control.family = list(link = "log"), # Control options
                                 control.predictor = list(A = INLA::inla.stack.A(stk_full), link = 1, compute = TRUE),  # Compute for marginals of the predictors.  Link to NULL for multiple likelihoods!
                                 control.compute = list(cpo = FALSE, waic = FALSE, config = TRUE, openmp.strategy	= 'huge' ), #model diagnostics and config = TRUE gives you the GMRF
                                 control.mode = list(theta = thetas, restart = TRUE), # To speed up use previous thetas
                                 # control.fixed = INLA::control.fixed(prec = list( initial = log(0.000001), fixed = TRUE)), # Added to see whether this changes GMRFlib convergence issues
                                 verbose = verbose, # To see the log of the model runs
                                 control.results = list(return.marginals.random = FALSE,
                                                        return.marginals.predictor = FALSE), # Don't predict marginals to save speed
                                 INLA::control.inla(int.strategy = "eb", # Empirical bayes for integration
                                                    strategy = 'simplified.laplace'
                                                    # https://groups.google.com/g/r-inla-discussion-group/c/hDboQsJ1Mls
                                 ),
                                 num.threads = parallel::detectCores() - 1
          )
          # Create a spatial prediction
          index.pred <- INLA::inla.stack.index(stk_full, 'pred_poipo')$data
          post <- fit_pred$summary.linear.predictor[index.pred, ]
          assertthat::assert_that(nrow(post)>0,
                                  nrow(post) == nrow(predcoords) ) # Check with cells in projection
          # TODO: Wrap Backtransform into a distribution dependent function
          post[,c('mean','0.05quant','0.5quant','0.95quant','mode')] <-
            exp(
              post[,c('mean','0.05quant','0.5quant','0.95quant','mode')]
            )
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

        # Definition of INLA Model object ----
        out <- bdproto(
              "INLA-Model",
              DistributionModel,
              model = model,
              fits = list(
                "fit_best" = fit_resp,
                "fit_pred" = fit_pred,
                "fit_best_equation" = master_form,
                "mesh"     = self$get_data('mesh'),
                "prediction" = prediction
                ),
              # Function to plot SPDE if existing
              plot_spde = function(self,...){
                if( length( self$fits$fit_best$size.spde2.blc ) == 1)
                {
                  # Get spatial projections from model
                  # FIXME: Potentially make the plotting of this more flexible
                  gproj <- INLA::inla.mesh.projector(self$get_data('mesh'),  dims = c(300, 300))
                  g.mean <- INLA::inla.mesh.project(gproj,
                                                    self$get_data('fit_pred')$summary.random$spatial.field$mean)
                  g.sd <- INLA::inla.mesh.project(gproj, self$get_data('fit_pred')$summary.random$spatial.field$sd)

                  # Out
                  r.m <- rasterFromXYZ(xyz = cbind(gproj$x,gproj$y ))
                  r.m[] <- as.vector(g.mean)
                  r.sd <- rasterFromXYZ(xyz = cbind(gproj$x,gproj$y ))
                  r.sd[] <- as.vector(g.sd)

                  # Plot
                  cols <- c("#00204DFF","#00336FFF","#39486BFF","#575C6DFF","#707173FF","#8A8779FF","#A69D75FF","#C4B56CFF","#E4CF5BFF","#FFEA46FF")
                  par(mfrow=c(1,2))
                  plot(raster::flip(r.m,direction = 'y'),col = cols, main = 'mean spatial effect')
                  plot(raster::flip(r.sd,direction = 'y'), main = 'sd spatial effect')
                  # raster::image(g.mean,col = cols, main = 'mean spatial effect')
                  # raster::image(r.sd, main = 'sd spatial effect')
                } else {
                  message('No spatial covariance in model specified.')
                }
              }
        )
        return(out)
      }
      ))
}
