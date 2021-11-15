#' @include bdproto-engine.R utils-inla.R bdproto-distributionmodel.R
NULL

#' Use inlabru as engine
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param optional_mesh A directly supplied [`INLA`] mesh (Default: NULL)
#' @param optional_projstk A directly supplied projection stack. Useful if projection stack is identical for multiple species (Default: NULL)
#' @param max.edge The largest allowed triangle edge length, must be in the same scale units as the coordinates
#' @param offset interpreted as a numeric factor relative to the approximate data diameter;
#' @param cutoff The minimum allowed distance between points on the mesh
#' @param proj_stepsize The stepsize in coordinate units between cells of the projection grid (Default: NULL)
#' @param timeout Specify a timeout for INLA models in sec. Afterwards it passed.
#' @param ... Other variables
#' @references https://becarioprecario.bitbucket.io/inla-gitbook/ch-spatial.html#ref-Simpsonetal:2016
#' @name engine_inlabru
NULL
#' @rdname engine_inlabru
#' @export
engine_inlabru <- function(x,
                        optional_mesh = NULL,
                        optional_projstk = NULL,
                        max.edge = c(1,5),
                        offset = c(1,1),
                        cutoff = 1,
                        proj_stepsize = NULL,
                        timeout = NULL,
                        ...) {

  # Check whether INLA package is available
  check_package('inlabru')
  if(!isNamespaceLoaded("inlabru")) { attachNamespace("inlabru");requireNamespace('inlabru') }
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
                          is.null(timeout) || is.numeric(timeout),
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

  # If time out is specified
  if(!is.null(timeout)) INLA::inla.setOption(fmesher.timeout = timeout)

  # Get barrier from the region polygon
  # TODO: Add this in addition to spatial field below, possibly specify an option to calculate this

  # Calculate area in km²
  ar <- suppressWarnings(
    mesh_area(mesh = mesh,region.poly = region.poly, variant = 'gpc2')
  )

  # Print a message in case there is already an engine object
  if(!is.Waiver(x$engine)) myLog('[Setup]','yellow','Replacing currently selected engine.')

  # Set engine in distribution object
  x$set_engine(
    bdproto(
      "INLABRU-Engine",
      Engine,
      name = "<INLABRU>",
      data = list(
        'mesh' = mesh,
        'mesh.area' = ar,
        'proj_stepsize' = proj_stepsize,
        'stk_pred' = optional_projstk
      ),
      # Generic plotting function for the mesh
      plot = function(self, assess = FALSE){
        if(assess){
          # For an INLA mesh assessment
          out <- INLA:::inla.mesh.assessment(
            mesh = self$get_data('mesh'),
            spatial.range = 3,
            alpha = 2,
            dims = c(300, 300)
          )
          # Convert to raster stack
          out <- raster::stack(
            sp::SpatialPixelsDataFrame( sp::coordinates(out), data = as.data.frame(out),
                                        proj4string = self$get_data('mesh')$crs )
          )

          raster::plot(out[[c('sd','sd.dev','edge.len')]],
                       col = c("#00204D","#00336F","#39486B","#575C6D","#707173","#8A8779","#A69D75","#C4B56C","#E4CF5B","#FFEA46")
          )
        } else {
          INLA:::plot.inla.mesh( self$get_data('mesh') )
        }
      },
      # Spatial latent function
      # https://groups.google.com/g/r-inla-discussion-group/c/eqMhlbwChkQ/m/m0b0PuzL-PsJ
      # Default SPDE prior
      # It computes the approximate diameter of the mesh, multiplies by 0.2 to get a value for the prior median range, and then transforms it to log-kappa scale by the formula
      # log(sqrt(8*nu)/range) where nu is alpha-dim/2.
      calc_latent_spatial = function(self,type = 'spde', alpha = 2,
                                     priors = NULL,
                                     polynames = NULL,
                                     varname = "spatial.field",
                                     ...){
        # Catch prior objects
        if(is.null(priors) || is.Waiver(priors)) priors <- NULL

        # For calculating iCAR process
        if(type == 'car'){
          # convert mesh to sf object
          ns <- mesh_as_sf(self$data$mesh)
          # Create adjacency matrix with queen's case
          nc.nb <- spdep::poly2nb(ns, queen = TRUE)
          #Convert the adjacency matrix into a file in the INLA format
          adjmat <- spdep::nb2mat(nc.nb,style = "B")
          adjmat <- as(adjmat, "dgTMatrix")
          # adjmat <- INLA::inla.graph2matrix(nc.nb)
          # Save the adjaceny matrix as output
          self$data$latentspatial <- adjmat
          self$data$s.index <- as.numeric(attr(nc.nb,varname))
        } else if(type=='spde'){
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
          self$data$s.index <- INLA::inla.spde.make.index(name = varname,
                                                          n.spde = self$data$latentspatial$n.spde,
                                                          n.group = 1,
                                                          n.repl = 1)
          # Security checks
          assertthat::assert_that(
            inherits(self$data$latentspatial,'inla.spde'),
            length(self$data$s.index$spatial.field) == self$data$mesh$n
          )
        } else if(type == 'poly'){
          # Save column names of polynomial transformed coordinates
          assertthat::assert_that(!is.null(polynames))
          self$data$latentspatial <- polynames
        }
        invisible()
      },
      # Get latent spatial equation bit
      # Set vars to 2 or larger to get copied spde's
      get_equation_latent_spatial = function(self, method, vars = 1){
        assertthat::assert_that(is.numeric(vars))
        if(method == 'spde'){
          assertthat::assert_that(inherits(self$data$latentspatial, 'inla.spde'),
                                  msg = 'Latent spatial has not been calculated.')
          # SPDE string
          if(vars == 1){
            ss <- paste0('f(spatial.field, model = ',method,')')
          } else {
            ss <- paste0("f(spatial.field",vars,", copy = \'spatial.field\', model = ",method,", fixed = TRUE)")
          }
          return(ss)

        } else if(method == 'car'){
          assertthat::assert_that(inherits(self$data$latentspatial,'dgTMatrix'),
                                  msg = 'Neighborhood matrix has not been calculated.')
          return(
            # BESAG model or BYM model to specify
            # BYM found to be largely similar to SPDE https://onlinelibrary.wiley.com/doi/pdf/10.1002/ece3.3081
            paste0('f(','spatial.field',', model = "bym", graph = ','adjmat',')')
          )
        }
      },
      # Main INLA training function ----
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
        # Messager
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Engine setup.')

        # Construct likelihoods
        lhl <- list()
        for(j in 1:length(model$biodiversity)){
          # Combine observed and predictors for the
          df <- cbind(
            data.frame(observed = model$biodiversity[[i]]$observations[['observed']]),
            model$biodiversity[[i]]$predictors
          )
          # Convert to Spatial points
          df <- sp::SpatialPointsDataFrame(coords = df[,c('x', 'y')],
                                              data = df[, names(df) %notin% c('x','y')],
                                              proj4string = self$data$mesh$crs
          )

          # Options for specifying link function of likelihood
          o <- inlabru::bru_options_get()
          o[['control.family']] <- list(link = ifelse(model$biodiversity[[i]]$family=='binomial', 'cloglog', 'default'))

          lh <- inlabru::like(formula = model$biodiversity[[i]]$equation,
                              family = model$biodiversity[[i]]$family,
                              data = df,
                              mesh = self$data$mesh,
                              # E = 1e6,
                              # ips = data@ips,
                              # Ntrials = attributes(data_points[[i]])$Ntrials,
                              # include = include[[i]]
                              options = o
                              )
          lhl[[j]] <- lh
        }

        # List of likelihoods
        self$set_data("likelihoods", inlabru::like_list(lhl) )

        # Set number of threads via set.Options
        inlabru::bru_safe_inla(quietly = TRUE)
        INLA::inla.setOption(num.threads = getOption('ibis.nthread'),
                             blas.num.threads = getOption('ibis.nthread'))

        # Set any bru options via verbosity of fitting
        inlabru::bru_options_set(bru_verbose = getOption('ibis.setupmessages'))
        inlabru::bru_options_set(control.inla = list(int.strategy = "eb"))
        invisible()
      },
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


        # Convert predictors to SpatialPixelsDataFrame as required for
        if(class(model$predictors) == 'data.frame') {
          preds <- sp::SpatialPointsDataFrame(coords = model$predictors[,c('x', 'y')],
                                                          data = model$predictors[, names(model$predictors) %notin% c('x','y')],
                                                          proj4string = self$data$mesh$crs
                                          )
          preds <- as(preds, 'SpatialPixelsDataFrame')
        } else stop('Predictors not in right format.')

        # Get likelihood
        likelihoods <- self$get_data("likelihoods")

        # Get options
        options <- inlabru::bru_options_get()
        # -------- #
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Starting fitting.')

        # Compute end of computation time
        settings$set('end.time', Sys.time())

        comp <- as.formula("~ Intercept(1) ")
        # Fit
        mod_bru <- inlabru::bru(components = comp,
                                likelihoods, options = options)


        # Definition of INLA Model object ----
        out <- bdproto(
          "INLABRU-Model",
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
          # Projection function
          project = function(self, newdata, mode = 'coef', backtransf = NULL){
            assertthat::assert_that('fit_best' %in% names(self$fits),
                                    is.data.frame(newdata) || is.matrix(newdata),
                                    mode %in% c('coef','sim','full'),
                                    assertthat::has_name(newdata,c('x','y'))
            )
            # Try and guess backtransformation
            if(is.null(backtransf)){
              fam <- self$get_data('fit_best')$.args$family
              backtransf <- ifelse(fam == 'poisson', exp, logistic)
            }

            if(mode == 'coef'){
              # We use the coefficient prediction
              out <- coef_prediction(mesh = self$get_data('mesh'),
                                     mod = self,
                                     type = 'mean',
                                     backtransf = backtransf
              )
            } else if(mode == 'sim'){
              # Simulate from posterior. Not yet coded
              stop('Simulation from posterior not yet there.')
            } else {
              stop('Full prediction not yet added.')
            }
            # Return result
            return(out)
          },
          # Partial response
          # FIXME: Create external function
          partial = function(self, x, variable, constant = NULL, variable_length = 100, plot = FALSE){
            # Goal is to create a sequence of value and constant and append to existing stack
            # Alternative is to create a model-matrix through INLA::inla.make.lincomb() and
            # model.matrix(~ vars, data = newDummydata) fed to make.lincomb
            # provide via lincomb = M to an INLA call.
            # Both should be identical

            # Check that provided model exists and variable exist in model
            mod <- self$get_data('fit_best')
            assertthat::assert_that(inherits(mod,'inla'),
                                    'model' %in% names(self),
                                    inherits(x,'BiodiversityDistribution'),
                                    length(variable)==1, is.character(variable),
                                    is.null(constant) || is.numeric(constant)
            )
            varn <- mod$names.fixed
            variable <- match.arg(variable, varn, several.ok = FALSE)
            assertthat::assert_that(variable %in% varn, length(variable)==1,!is.null(variable))

            # ------------------ #
            # Get all datasets with id in model. This includes the data stacks and integration stacks
            stk_inference <- lapply(
              x$engine$list_data()[grep(paste(names(model$biodiversity),collapse = '|'), x$engine$list_data())],
              function(z) x$engine$get_data(z))
            stk_inference <- do.call(INLA::inla.stack, stk_inference)
            # FIXME: Test that this works with SPDE present
            stack_data_resp <- INLA::inla.stack.data(stk_inference)
            # ------------------ #

            # If constant is null, calculate average across other values
            if(is.null(constant)){
              constant <- lapply(stack_data_resp, function(x) mean(x,na.rm = T))[varn[varn %notin% variable]]
            }
            # For target variable calculate range
            variable_range <- range(stack_data_resp[[variable]],na.rm = TRUE)

            stop('Not yet done!')
            # Create dummy data.frame
            dummy <- data.frame(observed = rep(NA, variable_length))
            dummy

            seq(variable_range[1],variable_range[2],length.out = variable_length)

            # # add sequence of data and na to data.frame. predict those
            # control.predictor = list(A = INLA::inla.stack.A(stk_full),
            #                          link = li, # Link to NULL for multiple likelihoods!
            #                          compute = TRUE),  # Compute for marginals of the predictors.

            print('Refitting model for partial effect')
            ufit <- INLA::inla(formula = as.formula(mod$.args$formula), # The specified formula
                               data  = stk_inference,  # The data stack
                               quantiles = c(0.05, 0.5, 0.95),
                               E = INLA::inla.stack.data(stk_inference)$e, # Expectation (Eta) for Poisson model
                               Ntrials = INLA::inla.stack.data(stk_inference)$Ntrials,
                               family = mod$.args$family,   # Family the data comes from
                               control.family = mod$.args$control.family, # Control options
                               control.predictor = mod$.args$control.predictor,  # Compute for marginals of the predictors.
                               control.compute = mod$.args$control.compute,
                               control.fixed = mod$.args$control.fixed,
                               verbose = FALSE, # To see the log of the model runs
                               control.inla = mod$.args$control.inla,
                               num.threads = getOption('ibis.nthread')
            )
            control.predictor = list(A = INLA::inla.stack.A(stk_inference))

            # Plot and return result
          },
          # Function to plot SPDE if existing
          plot_spatial = function(self, dim = c(300,300), kappa_cor = FALSE, ...){
            assertthat::assert_that(is.vector(dim))
            if( length( self$fits$fit_best$size.spde2.blc ) == 1)
            {
              # Get spatial projections from model
              # FIXME: Potentially make the plotting of this more flexible
              gproj <- INLA::inla.mesh.projector(self$get_data('mesh'),  dims = dim)
              g.mean <- INLA::inla.mesh.project(gproj,
                                                self$get_data('fit_best')$summary.random$spatial.field$mean)
              g.sd <- INLA::inla.mesh.project(gproj, self$get_data('fit_best')$summary.random$spatial.field$sd)

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
              # Mask with prediction if exists
              if(!is.null(self$get_data('prediction'))){
                spatial_field <- raster::resample(spatial_field, self$get_data('prediction')[[1]])
                spatial_field <- raster::mask(spatial_field, self$get_data('prediction')[[1]])
              }

              # -- #
              if(kappa_cor){
                # Also build correlation fun
                # Get SPDE results
                spde_results <- INLA::inla.spde2.result(
                  inla = self$get_data('fit_best'),
                  name = 'spatial.field',
                  spde = self$get_data('spde'),
                  do.transfer = TRUE)

                # Large kappa (inverse range) equals a quick parameter change in space.
                # Small kappa parameter have much longer, slower gradients.
                Kappa <- INLA::inla.emarginal(function(x) x, spde_results$marginals.kappa[[1]])
                sigmau <- INLA::inla.emarginal(function(x) sqrt(x), spde_results$marginals.variance.nominal[[1]])
                r <- INLA::inla.emarginal(function(x) x, spde_results$marginals.range.nominal[[1]])

                # Get Mesh and distance between points
                mesh <- self$get_data('mesh')
                D <- as.matrix( dist(mesh$loc[, 1:2]) )

                # Distance vector.
                dis.cor <- data.frame(distance = seq(0, max(D), length = 100))
                # Maximum distance by quarter of extent
                dis.max <- abs((xmin(self$get_data('prediction')) - xmax(self$get_data('prediction')) ) / 2)  # Take a quarter of the max distance

                # Modified Bessel function to get correlation strength
                dis.cor$cor <- as.numeric((Kappa * dis.cor$distance) * base::besselK(Kappa * dis.cor$distance, 1))
                dis.cor$cor[1] <- 1
                # ---
                # Build plot
                layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
                plot(dis.cor$cor ~ dis.cor$distance, type = 'l', lwd = 3,
                     xlab = 'Distance (proj. unit)', ylab = 'Correlation', main = paste0('Kappa: ', round(Kappa,2) ) )
                abline(v = dis.max,lty = 'dotted')
                plot(spatial_field[[1]],col = ibis_colours[['viridis_cividis']], main = 'mean spatial effect')
                plot(spatial_field[[2]], main = 'sd spatial effect')
              } else {
                # Just plot the SPDE
                par(mfrow=c(1,2))
                plot(spatial_field[[1]],col = ibis_colours[['viridis_cividis']], main = 'mean spatial effect')
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
