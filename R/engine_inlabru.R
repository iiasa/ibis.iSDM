#' @include bdproto-engine.R utils-inla.R bdproto-distributionmodel.R
NULL

#' Use inlabru as engine
#'
#' @description Model components are specified with general inputs and mapping methods to the
#' latent variables, and the predictors are specified via general R expressions,
#' with separate expressions for each observation likelihood model in multi-likelihood models.
#' The inlabru engine - similar as the [`engine_inla`] function acts a wrapper for INLA,
#' albeit \code{"inlabru"} has a number of convenience functions implemented that make in particular predictions
#' with new data much more straight forward (e.g. via posterior simulation instead of fitting).
#' Since more recent versions \code{"inlabru"} also supports the addition of multiple likelihoods, therefore
#' allowing full integrated inference.
#' @details
#' All \code{INLA} engines require the specification of a mesh that needs to be provided to the
#' \code{"optional_mesh"} parameter. Otherwise the mesh will be created based on best guesses of the
#' data spread. A good mesh needs to have triangles as regular as possible in size and shape: equilateral.
#'
#' [*] \code{"max.edge"}: The largest allowed triangle edge length, must be in the same scale units as the coordinates
#' Lower bounds affect the density of triangles
#' [*] \code{"offset"}: The automatic extension distance of the mesh
#' If positive: same scale units. If negative, interpreted as a factor relative to the approximate data diameter
#' i.e., a value of -0.10 will add a 10% of the data diameter as outer extension.
#' [*] \code{"cutoff"}: The minimum allowed distance between points,
#' it means that points at a closer distance than the supplied value are replaced by a single vertex.
#' it is critical when there are some points very close to each other, either for point locations or in the
#' domain boundary.
#' [*] \code{"proj_stepsize"}: The stepsize for spatial predictions, which affects the spatial grain of any outputs
#' created.
#'
#' Priors can be set via [INLAPrior].
#' @note
#' **How INLA Meshes are generated, substantially influences prediction outcomes. See Dambly et al. (2023).**
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param optional_mesh A directly supplied \code{"INLA"} mesh (Default: \code{NULL})
#' @param max.edge The largest allowed triangle edge length, must be in the same scale units as the coordinates.
#' Default is an educated guess (Default: \code{NULL}).
#' @param offset interpreted as a numeric factor relative to the approximate data diameter.
#' Default is an educated guess (Default: \code{NULL}).
#' @param cutoff The minimum allowed distance between points on the mesh.
#' Default is an educated guess (Default: \code{NULL}).
#' @param proj_stepsize The stepsize in coordinate units between cells of the projection grid (Default: \code{NULL})
#' @param strategy Which approximation to use for the joint posterior. Options are \code{"auto"} ("default"), \code{"adaptative"},
#'  \code{"gaussian"}, \code{"simplified.laplace"} & \code{"laplace"}.
#' @param int.strategy Integration strategy. Options are \code{"auto"},\code{"grid"}, \code{"eb"} ("default") & \code{"ccd"}.
#' @param area Accepts a [`character`] denoting the type of area calculation to be done on the mesh (Default: \code{'gpc2'}).
#' @param timeout Specify a timeout for INLA models in sec. Afterwards it passed.
#' @param type The mode used for creating posterior predictions. Either summarizing the linear \code{"predictor"} or \code{"response"} (Default: \code{"response"}).
#' @param ... Other variables
#' @references
#' * Bachl, F. E., Lindgren, F., Borchers, D. L., & Illian, J. B. (2019). inlabru: an R package for Bayesian spatial modelling from ecological survey data. Methods in Ecology and Evolution, 10(6), 760-766.
#' * Simpson, Daniel, Janine B. Illian, S. H. Sørbye, and Håvard Rue. 2016. “Going Off Grid: Computationally Efficient Inference for Log-Gaussian Cox Processes.” Biometrika 1 (103): 49–70.
#' * Dambly, L. I., Isaac, N. J., Jones, K. E., Boughey, K. L., & O'Hara, R. B. (2023). Integrated species distribution models fitted in INLA are sensitive to mesh parameterisation. Ecography, e06391.
#' @source [https://inlabru-org.github.io/inlabru/articles/](https://inlabru-org.github.io/inlabru/articles/)
#' @family engine
#' @returns An engine.
#' @examples
#' \dontrun{
#' # Add inlabru as an engine
#' x <- distribution(background) |> engine_inlabru()
#' }
#' @name engine_inlabru
NULL
#' @rdname engine_inlabru
#' @export
engine_inlabru <- function(x,
                        optional_mesh = NULL,
                        max.edge = NULL,
                        offset = NULL,
                        cutoff = NULL,
                        proj_stepsize = NULL,
                        strategy = "auto",
                        int.strategy = "eb",
                        area = "gpc2",
                        timeout = NULL,
                        type = "response",
                        ...) {

  # Check whether INLA package is available
  check_package('inlabru')
  if(!isNamespaceLoaded("inlabru")) { attachNamespace("inlabru");requireNamespace('inlabru') }
  check_package('INLA')
  if(!isNamespaceLoaded("INLA")) { attachNamespace("INLA");requireNamespace('INLA') }

  # assert that arguments are valid
  assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                          inherits(x$background,'sf'),
                          inherits(optional_mesh,'inla.mesh') || is.null(optional_mesh),
                          is.vector(max.edge) || is.null(max.edge),
                          (is.vector(offset) || is.numeric(offset)) || is.null(offset),
                          is.numeric(cutoff) || is.null(cutoff),
                          is.null(timeout) || is.numeric(timeout),
                          is.character(strategy),
                          is.character(int.strategy),
                          is.character(area),
                          is.character(type),
                          is.null(proj_stepsize) || is.numeric(proj_stepsize)
  )

  # Match strategy
  strategy <- match.arg(strategy, c('auto', 'gaussian', 'simplified.laplace', 'laplace', 'adaptive'), several.ok = FALSE)
  int.strategy <- match.arg(int.strategy, c('auto', 'ccd', 'grid', 'eb'), several.ok = FALSE)
  area <- match.arg(area, c("gpc", "gpc2", "km"), several.ok = FALSE)
  type <- match.arg(type, c("response", "predictor"), several.ok = FALSE)

  # Set inlabru options for strategies. These are set globally
  inlabru::bru_options_set(control.inla = list(strategy = strategy,
                                               int.strategy = int.strategy))

  # Set the projection mesh
  if(inherits(optional_mesh,'inla.mesh')) {
    # Load a provided on
    mesh <- optional_mesh
    # Security check for projection and if not set, use the one from background
    if(is.null(mesh$crs))  mesh$crs <- sp::CRS( sp::proj4string(region.poly) )

    # Convert the study region
    region.poly <- methods::as(sf::st_geometry(x$background), "Spatial")

    # Calculate area
    ar <- suppressWarnings(
      mesh_area(mesh = mesh, region.poly = region.poly, variant = area)
    )
  } else {
    mesh <- new_waiver()
    ar <- new_waiver()
  }

  # Collate other parameters in a specific object
  params <- list(
    max.edge = max.edge,
    offset = offset,
    cutoff = cutoff,
    proj_stepsize = proj_stepsize,
    type = type,
    area = area,
    strategy = strategy,
    int.strategy = int.strategy,
    ...
  )

  # If time out is specified
  if(!is.null(timeout)) INLA::inla.setOption(fmesher.timeout = timeout)

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
        'params' = params
      ),
      # Function to create a mesh
      create_mesh = function(self, model){
        assertthat::assert_that(is.list(model),
                                "background" %in% names(model))
        # Check if mesh is already present, if so use it
        if(!is.Waiver(self$get_data("mesh"))) return()
        # Create a new mesh based on the available data

        # Get parameters
        params <- self$get_data("params")

        # Convert the study region
        region.poly <- methods::as(sf::st_geometry(model$background), "Spatial")

        # Convert to boundary object for later
        suppressWarnings(
          bdry <- INLA::inla.sp2segment(
            sp = region.poly,
            join = TRUE,
            crs = INLA::inla.CRS(projargs = sp::proj4string(region.poly))
          )
        )
        bdry$loc <- INLA::inla.mesh.map(bdry$loc)

        # Get all coordinates of observations
        locs <- collect_occurrencepoints(model, include_absences = FALSE)

        # Try and infer mesh parameters if not set
        if(is.null(params$max.edge)){
          # A good guess here is usally a max.edge of between 1/3 to 1/5 of the spatial range.
          max.edge <- c(diff(range(locs[,1]))/(3*5) , diff(range(locs[,1]))/(3*5) * 2)
          params$max.edge <- max.edge
        }
        if(is.null(params$offset)){
          # Check whether the coordinate system is longlat
          if( sf::st_is_longlat(bdry$crs) ){
            # Specify offset as 1/100 of the boundary distance
            offset <- c( diff(range(bdry$loc[,1]))*0.01,
                         diff(range(bdry$loc[,1]))*0.01)
          } else {
            offset <- c( diff(range(bdry$loc[,1]))*0.01,
                         diff(range(bdry$loc[,1]))*0.01)
          }
          params$offset <- offset
        }

        if(is.null(params$cutoff)){
          # Specify as minimum distance between y coordinates
          # Thus capturing most points on this level
          # otherwise set to default
          val <- min(abs(diff(locs[,2])))
          cutoff <- ifelse(val == 0, 1e-12, val)
          params$cutoff <- cutoff
        }

        suppressWarnings(
          mesh <- INLA::inla.mesh.2d(
            # Point localities
            loc = locs,
            # Boundary object
            boundary = bdry,
            # Mesh Parameters
            max.edge = params$max.edge,
            offset = params$offset,
            cutoff = params$cutoff,
            # Define the CRS
            crs = bdry$crs
          )
        )
        # Calculate area
        # ar <- suppressMessages(
        #         suppressWarnings(
        #           mesh_area(mesh = mesh, region.poly = region.poly, variant = params$area)
        #         )
        #       )
        # 06/01/2023: This should work and is identical to inlabru::ipoints
        ar <- suppressWarnings(
          inlabru::ipoints(samplers = mesh)$weight |> as.vector()
        )
        assertthat::assert_that(length(ar) == mesh$n)

        # Now set the output
        self$set_data("mesh", mesh)
        self$set_data("mesh.area", ar)

        invisible()
      },
      # Generic plotting function for the mesh
      plot = function(self, assess = FALSE){
        if(is.Waiver(self$get_data('mesh'))) stop("No mesh found!")
        if(assess){
          # For an INLA mesh assessment
          out <- INLA::inla.mesh.assessment(
            mesh = self$get_data('mesh'),
            spatial.range = 3,
            alpha = 2,
            dims = c(300, 300)
          )
          # Convert to raster stack
          out <- terra::rast(
            sp::SpatialPixelsDataFrame( sp::coordinates(out), data = as.data.frame(out),
                                        proj4string = terra::crs(self$get_data('mesh')$crs) )
          )

          terra::plot(out[[c('sd','sd.dev','edge.len')]],
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
                                     varname = "spatial.field1",
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
          adjmat <- methods::as(adjmat, "dgTMatrix")
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
            length(self$data$s.index[[1]]) == self$data$mesh$n
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
      get_equation_latent_spatial = function(self, method, vars = 1, separate_spde = FALSE){
        assertthat::assert_that(is.numeric(vars))
        if(method == 'spde'){
          assertthat::assert_that(inherits(self$data$latentspatial, 'inla.spde'),
                                  msg = 'Latent spatial has not been calculated.')
          # SPDE string
          if(separate_spde){
            ss <- paste0("f(spatial.field",vars,", model = ",method,")")
          } else {
            if(vars >1){
              ss <- paste0("f(spatial.field",vars,", copy = \'spatial.field1\', model = ",method,", fixed = TRUE)")
            } else {
              ss <- paste0("f(spatial.field",vars,", model = ",method,")")
            }
          }
          return(ss)

        } else if(method == 'car'){
          assertthat::assert_that(inherits(self$data$latentspatial,'dgTMatrix'),
                                  msg = 'Neighborhood matrix has not been calculated.')
          return(
            # BESAG model or BYM model to specify
            # BYM found to be largely similar to SPDE https://onlinelibrary.wiley.com/doi/pdf/10.1002/ece3.3081
            paste0('f(','spatial.field1',', model = "bym", graph = ','adjmat',')')
          )
        }
      },
      calc_integration_points = function(self, model, mode = 'stack'){
        # Mode cox process integration
        if(mode == 'cp'){
          # Create integration points by using the mesh as sampler
          suppressWarnings(
            ips <- inlabru::ipoints(
              samplers = self$get_data('mesh')
            )
          )
          # Extract predictors add to integration point data
          d <- get_rastervalue(coords = ips@coords,
                               env = model$predictors_object$get_data(df = FALSE),
                               rm.na = FALSE)
          for (cov in model$predictors_names) ips@data[,cov] <- d[,cov]
          ips@data$Intercept <- 1
          ips <- subset(ips, stats::complete.cases(ips@data)) # Necessary as some integration points can fall outside land area
          # Return results
          return(ips)
        } else if(mode == 'stack'){
          # Use INLA make stack function instead. Useful for poisson created data so that
          # integration points are created as absence. Numerically inefficient though compared to cp
          # FIXME: Ideally sample from a provided pseudo-background
          istk <- inla_make_integration_stack(mesh = self$get_data('mesh'),
                                              mesh.area = self$get_data("mesh.area"),
                                              model = model,
                                              id = "istack",
                                              joint = FALSE)
          ips <- cbind(istk$data$data, istk$effects$data) # Combine observations and stack
          ips <- subset(ips, stats::complete.cases(ips[,c("x", "y")])) # Remove NA coordinates
          # Convert to sp
          ips <- sp::SpatialPointsDataFrame(coords = ips[,c('x', 'y')],
                                            data = ips[, names(ips) %notin% c('x','y')],
                                            proj4string = self$get_data('mesh')$crs
          )
          # Select only the predictor names
          ips <- subset(ips, select = c("observed", "Intercept", "e", model$predictors_names))
          ips <- subset(ips, stats::complete.cases(ips@data))
          abs_E <- ips$e; ips$e <- NULL
          # Return list of result
          return(list(ips = ips, E = abs_E))
        }
      },
      # Main inlabru setup ----
      # Setup computation function
      setup = function(self, model, settings, ...){
        assertthat::assert_that(
          'background' %in% names(model),
          'biodiversity' %in% names(model),
          all( sapply(model$biodiversity, function(x) is.formula(x$equation)) ),
          length(model$biodiversity)>=1,
          msg = 'Some internal checks failed while setting up the model.'
        )
        # Messenger
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Engine setup.')

        # Construct likelihoods for each entry in the dataset
        lhl <- list()
        for(j in 1:length(model$biodiversity)){
          # Combine observed and predictors for the
          df <- cbind(
            data.frame(observed = model$biodiversity[[j]]$observations[['observed']]),
            model$biodiversity[[j]]$predictors
          )
          # Convert to Spatial points
          df <- sf::st_as_sf(df,
                             coords = c('x', 'y'),
                             # data = df[, names(df) %notin% c('x','y')],
                             crs = sf::st_crs(self$get_data('mesh')$crs)
          ) |> as("Spatial")
          assertthat::assert_that(inherits(df, "sf") || inherits(df, "Spatial"))

          # Options for specifying link function of likelihood
          o <- inlabru::bru_options_get()
          # Data type specific. Currently only binomial and poisson supported
          # FIXME: Code below does not work as intended. Worked in earlier versions. To be debugged later!
          # FIXME: Apparently the code works when an SPDE is added, thus cp seems to rely on the matter covs
          # if(model$biodiversity[[j]]$type == 'poipo'){
          #   ips <- self$calc_integration_points(model, mode = 'cp')
          #
          #   # Log gaussian cox process
          #   lh <- inlabru::like(formula = stats::update.formula(model$biodiversity[[j]]$equation, "coordinates ~ ."),
          #                       # include = model$biodiversity[[j]]$predictors_names,
          #                       family = "cp",
          #                       data = df,
          #                       # samplers = methods::as(model$background,"Spatial"),
          #                       # domain = list(coordinates = self$get_data("mesh")),
          #                       ips = ips,
          #                       options = o
          #   )
          #   assertthat::assert_that(sum(lh$response_data$BRU_response_cp>0)>2,
          #                           msg = "Found issues with aligning coordinates within the domain most likely.")
          # # If not poipo but still poisson, prepare data as follows
          # }
          if(model$biodiversity[[j]]$family == "poisson"){
            # Calculate integration points for PPMs and to estimation data.frame
            ips <- self$calc_integration_points(model, mode = 'stack')
            abs_E = ips$E; ips <- ips$ips
            assertthat::assert_that(all(colnames(ips) %in% colnames(df)))
            new <- sp::rbind.SpatialPointsDataFrame(
              df[,c('observed', 'Intercept', model$biodiversity[[j]]$predictors_names)],
              ips[,c('observed', 'Intercept', model$biodiversity[[j]]$predictors_names)])
            # Formulate the likelihood
            lh <- inlabru::like(formula = model$biodiversity[[j]]$equation,
                                family = model$biodiversity[[j]]$family,
                                data = new, # Combine presence and absence information
                                mesh = self$get_data('mesh'),
                                E = c(model$biodiversity[[j]]$expect, abs_E), # Combine Exposure variants
                                # include = include[[i]], # Don't need this as all variables included in equation
                                options = o
            )
          } else if(model$biodiversity[[j]]$family == "binomial"){
            # Set likelihood to cloglog for binomial following Simpson 2016 if multiple likelihoods
            if(length(model$biodiversity)>1){
              o[['control.family']] <- list(link = ifelse(model$biodiversity[[j]]$family=='binomial', 'cloglog', 'default'))
            }

            # Formulate the likelihood
            lh <- inlabru::like(formula = model$biodiversity[[j]]$equation,
                                family = model$biodiversity[[j]]$family,
                                data = df, # Combine presence and absence information
                                mesh = self$get_data('mesh'),
                                Ntrials = model$biodiversity[[j]]$expect,
                                # include = include[[i]], # Don't need this as all variables in equation are included
                                options = o
            )
          }
          # Add to list
          lhl[[j]] <- lh
        }

        # List of likelihoods
        self$set_data("likelihoods", inlabru::like_list(lhl) )

        # --- #
        # Defining the component function
        if(length(model$biodiversity)>1){
          comp <- stats::as.formula(
            paste(' ~ 0 + Intercept(1) ',
                  ifelse(model$biodiversity[[1]]$use_intercept,
                         paste("+",paste0('Intercept_',
                                make.names(tolower(sapply( model$biodiversity, function(x) x$name ))),'_', # Make intercept from name
                                sapply( model$biodiversity, function(x) x$type ),"(1)",collapse = ' + ')
                               ),
                         ""
                          )
                  )
          )
        } else {
          comp <- stats::as.formula(
            paste0( "~ Intercept(1)")
          )
        }

        # Add Offset if set
        if(!is.Waiver(model$offset)){
          ovn <- "spatial_offset"
          comp <- stats::update.formula(comp,
                                 paste(c(' ~ . +', paste0(ovn,'(main = ', ovn, ', model = "offset")')), collapse = " ")
          )
        }

        # --- #
        # Get unified predictors from likelihoods
        pn <- lapply(lhl, function(x) all.vars(x$formula) ) |> (\(.) do.call(c,.))()  |>  unique()
        # pn <- lapply(model$biodiversity, function(x) x$predictors_names ) |> do.call(c,.) |> unique()
        pn <- pn[grep("Intercept|coordinates", pn, invert = TRUE)]
        assertthat::assert_that(length(pn)>0)
        model$predictors_types <- model$predictors_types[which(model$predictors_types$predictors %in% pn),]

        # Add Predictors to component
        for(i in 1:nrow(model$predictors_types)){
          # For numeric
          if(model$predictors_types$type[i] == 'numeric' | model$predictors_types$type[i] == 'integer') {
            # Built component
            if(settings$get('only_linear') == FALSE){
              # if there are less than 50 unique values, create linear variable instead
              if(length(unique(model$predictors[,i])) > 50){
                m <- paste0("rw1","__",model$predictors_types$predictors[i])
              } else m <- "linear"
            } else { m <- 'linear' }

            # Specify priors if set
            if(!is.Waiver(model$priors)){
              # If a prior has been specified
              if(any(model$priors$varnames() == model$predictors_types$predictors[i])){
                vn <- model$priors$varnames()[which(model$priors$varnames() == model$predictors_types$predictors[i])]
                ty <- model$priors$types()[names(vn)]
                if(ty %in% c("gaussian", "normal")){
                  pp <- paste0(c(
                    ', mean.linear = ', model$priors$get(vn)[1],', ',
                    'prec.linear = ', model$priors$get(vn)[2],''
                  ),collapse = "" )
                } else if(ty == "clinear"){
                  pp <- paste0("hyper = list(theta = c(prior = \'clinear\', param = c(",
                               model$priors$get(vn)[1],", ",model$priors$get(vn)[2],")))")
                  m <- "clinear"
                }
              } else {pp <- "" }
            } else { pp <- "" }
            if( m!= "linear" ){
                # Could add RW effects with pc priors. PC priors are on the KL distance (difference between probability distributions), P(sigma >2)=0.05
                # Default is a loggamma prior with mu 1, 5e-05. Better would be 1, 0.5 following Caroll 2015, so we define it like this here
                pp <- ', hyper = list(theta = list(prior = \'loggamma\', param = c(1, 0.5)))'
            }
            comp <- stats::update.formula(comp,
                                   paste(' ~ . +', paste0(model$predictors_types$predictors[i],'(main = ', model$predictors_types$predictors[i],
                                                          pp,', model = "',m,'")'), collapse = " ")
            )
          } else if( model$predictors_types$type[i] == "factor"){
            # factor_full uses the full factor. fact_contrast uses the first level as reference
            # Built component
            comp <- update(comp,
                           paste(c(' ~ . + ', paste0(model$predictors_types$predictors[i],'(main = ', model$predictors_types$predictors[i], ', model = "factor_contrast")')), collapse = " ")
            )
          }
        }

        # Add spatial effect if set
        if("latentspatial" %in% self$list_data() ){
          spde <- self$get_data("latentspatial")
          assertthat::assert_that(inherits(spde, "inla.spde2"))
          if(inherits(spde, "inla.spde") ){
            for(i in 1:length(model$biodiversity)){
              # Add spatial component term
              comp <- stats::update.formula(comp,
                             paste0(c("~ . + "),
                                    paste0("spatial.field", i,
                                           "(main = geometry,",
                                           ifelse( grep('copy', self$get_equation_latent_spatial('spde', vars = i))==1,
                                                     " copy = \'spatial.field1\', fixed = TRUE,",
                                                     ""),
                                           "model = spde)"
                                    )
                             )
              )
            }
          } else {
            # FIXME: Make this more generic so that other latent effects are supported
            stop("Non-SPDE effects not yet implemented")
          }
        }
        # Set component
        self$set_data("components", comp)

        # Set number of threads via set.Options
        inlabru::bru_safe_inla(quietly = TRUE)
        INLA::inla.setOption(num.threads = getOption('ibis.nthread'),
                             blas.num.threads = getOption('ibis.nthread')
                             )

        # Set any other bru options via verbosity of fitting
        inlabru::bru_options_set(bru_verbose = settings$get('verbose'))
        # Newer inlabru versions support quantiles
        if(utils::packageVersion("inlabru") > '2.5.2'){
          inlabru::bru_options_set(quantiles = c(0.05, 0.5, 0.95))
        }
        invisible()
      },
      train = function(self, model, settings) {
        # Check that all inputs are there
        assertthat::assert_that(
          inherits(settings,'Settings'),
          is.list(model),length(model)>1,
          # Check that model id and setting id are identical
          settings$modelid == model$id,
          inherits(self$get_data("likelihoods"), 'list')
        )

        # Get likelihood
        likelihoods <- self$get_data("likelihoods")

        # Get model components
        comp <- self$get_data("components")

        # Get params
        params <- self$get_data("params")

        # Recreate non-linear variables in case they are set
        if(settings$get('only_linear') == FALSE){
          # TODO: Bypass grouping for now until this has been figured out
          m = INLA::inla.models()
          m$latent$rw1$min.diff = NULL
          assign("inla.models", m, INLA::inla.get.inlaEnv())

          for(i in 1:nrow(model$predictors_types)){
            # if there are less than 50 unique values, create linear variable instead
            if(length(unique(model$predictors[,i])) > 50){
              # Create a one-dimensional array
              m <- INLA::inla.mesh.1d(
                seq(min(model$predictors[,i],na.rm = TRUE),
                    max(model$predictors[,i],na.rm = TRUE), length.out = 100),
                degree = 1)
              m <- INLA::inla.spde2.matern(m)
              # Internally assign
              assign(x = paste0("rw1","__",model$predictors_types$predictors[i]),
                     value = m )
              rm(m)
            }
          }
        }

        # Get spatial effect if existent
        if("latentspatial" %in% self$list_data() ){
          spde <- self$get_data("latentspatial")
          assertthat::assert_that(exists("spde"),
                                  inherits(spde, "inla.spde2")
          )
        }
        # Get options
        options <- inlabru::bru_options_get()
        assertthat::assert_that(inlabru::bru_options_check(options))
        # -------- #
        if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Starting fitting.')

        if( settings$get(what='optim_hyperparam')){
          if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Performing incremental variable selection...')

          # Catch all variables with set priors and keep them!
          if(!is.Waiver(model$priors)) keep <- as.character(model$priors$varnames()) else keep <- NULL

          te <- attr(stats::terms.formula(comp), 'term.label')
          test_form <- comp
          # Remove variables that are never removed
          if(!is.null(keep)){
            test_form <- stats::update.formula(test_form, paste0(". ~ . - ",
                                                          paste0(
                                                            grep(pattern = paste0(keep, collapse = '|'),x = te, value = TRUE ),
                                                            collapse = "-"
                                                          ))
                                        )
            te <- te[grep(pattern = paste0(keep,collapse = '|'),x = te, invert = TRUE, fixed = TRUE )]
          }
          te <- te[grep('Intercept',te,ignore.case = T,invert = T)]
          # --- #
          # Iterate through unique combinations of variables backwards
          pb <- progress::progress_bar$new(total = length(te),format = "Backwards eliminating variables... :spin [:elapsedfull]")
          o <- options; o$bru_verbose <- FALSE
          o$bru_max_iter <- 2 # Use only two iteration max for the variable selection
          not_found <- TRUE
          while(not_found) {
            pb$tick()
            # --- #
            # Base Model #
            fit <- try({
              inlabru::bru(components = test_form,
                                      likelihoods,
                                      options = o)
            },silent = TRUE)
            if("error" %in% names(fit)) {not_found <- FALSE;next()}

            results_base <- data.frame(form = deparse1(test_form),
                            converged = fit$ok,
                            waic = fit$waic$waic,
                            dic = fit$dic$dic,
                            mean.deviance = fit$dic$mean.deviance )
            results <- data.frame()

            # Formula terms
            te <- attr(stats::terms.formula(test_form), 'term.label')
            te_int <- te[grep('Intercept',te,ignore.case = T, invert = F)] # capture intercept(s)
            te <- te[grep('Intercept',te,ignore.case = T, invert = T)]
            assertthat::assert_that(length(te) > 0, length(te_int) > 0)

            # Now for each term in variable list
            for(vars in te){
              # New formula
              new_form <- update(test_form, paste0('~ . - ',vars ))
              ll <- likelihoods

              try({fit <- inlabru::bru(components = new_form,
                             ll,
                             options = o)
              },silent = TRUE)
              if("error" %in% names(fit)){
                results <- rbind(results,
                                 data.frame( form = deparse1(new_form),
                                             converged = FALSE,
                                             waic = NA, dic = NA, mean.deviance = NA )
                )
              } else {
                results <- rbind(results, data.frame(form = deparse1(new_form),
                                                     converged = fit$ok,
                                                     waic = fit$waic$waic,
                                                     dic = fit$dic$dic,
                                                     mean.deviance = fit$dic$mean.deviance )
                )
              }
              rm(fit)
            } # End of loop

            if(!is.na(results_base$dic) || nrow(results) > 0) {
              # Now check whether any of the new models are 'better' than the full model
              if(results_base$dic <= min(results$dic, na.rm = TRUE)){
                not_found <- FALSE
                best_found <- results_base$form
              } else {
                # Otherwise continue get best model
                test_form <- stats::as.formula(results$form[which.min(results$dic)])
              }
              rm(results_base, results)
            } else {
              # Check whether formula is empty, if yes, set to not_found to FALSE
              te <- attr(stats::terms.formula(test_form),'term.label')
              if(length(te)<=4){
                not_found <- FALSE
                best_found <- test_form
              }
            }

          } # End of While loop
          # Make sure to add kept variables back
          if(!is.null(keep)){
            te <- attr(stats::terms.formula(comp),'term.label')
            best_found <- stats::update.formula(best_found, paste0(". ~ . + ",
                                                          paste0(
                                                            grep(pattern = paste0(keep, collapse = '|'),x = te, value = TRUE ),
                                                            collapse = "+"
                                                          ))
            )
          }
          # Replace component to be tested with best found
          comp <- stats::as.formula(best_found)
        }

        # --- #
        # Fitting bru model
        try({
          fit_bru <- inlabru::bru(components = comp,
                                  likelihoods,
                                  options = options)
        }, silent = FALSE)
        # --- #

        # Security checks
        if(!exists("fit_bru")){
          stop('Model did not converge. Try to simplify structure and check priors!')
        }
        if(is.null(fit_bru$names.fixed)) stop('Model did not converge. Try to simplify structure and check priors!')

        if(!settings$get('inference_only')){
          # Messenger
          if(getOption('ibis.setupmessages')) myLog('[Estimation]','green','Starting prediction.')

          covs <- model$predictors_object$get_data(df = FALSE)
          covs <- covs[[ which(names(covs) %in% fit_bru$names.fixed) ]]

          # Build coordinates
          suppressWarnings(
            preds <- inla_predpoints(mesh = self$get_data('mesh'),
                                     background = model$background,
                                     cov = covs,
                                     proj_stepsize = self$get_data('proj_stepsize'),
                                     spatial = TRUE
            )
          )
          # Clamp?
          if( settings$get("clamp") ) preds@data <- clamp_predictions(model, preds@data)

          # Set target variables to bias_value for prediction if specified
          if(!is.Waiver(settings$get('bias_variable'))){
            for(i in 1:length(settings$get('bias_variable'))){
              if(settings$get('bias_variable')[i] %notin% names(preds)){
                if(getOption('ibis.setupmessages')) myLog('[Estimation]','red','Did not find bias variable in prediction object!')
                next()
              }
              preds[[settings$get('bias_variable')[i]]] <- settings$get('bias_value')[i]
            }
          }
          # --- #
          # Define formula
          if(params$type == "response"){
            # Transformation to use for prediction scale
            # FIXME: This assumes no custom link function has been set!
            fun <- ifelse(length(model$biodiversity) == 1 && model$biodiversity[[1]]$type == 'poipa', "logistic", "exp")
          } else {
            fun <- "" # Linear predictor
          }

          # Get variables for inlabru
          if(length(model$biodiversity)>1){
            vn <- lapply(model$biodiversity, function(x) x$predictors_names) |> (\(.) do.call(c, .))() |> unique()
            ii <- paste("Intercept",
                        # # If multiple datasets, remove intercept
                        ifelse(length(model$biodiversity)>1,"+ 0", ""),
                        ifelse(model$biodiversity[[1]]$use_intercept,
                               paste("+",paste0('Intercept_',
                                                make.names(tolower(sapply( model$biodiversity, function(x) x$name ))),'_', # Make intercept from name
                                                sapply( model$biodiversity, function(x) x$type ),collapse = ' + ')
                               ),
                               ""
                               )
            )
            # Assert that variables are used in the likelihoods
            assertthat::assert_that(
              all( vn %in% (lapply(likelihoods, function(x) all.vars(x$formula) ) |>  (\(.) do.call(c, .))() |> unique() ) )
            )
          } else {
            # vn <- sapply(model$biodiversity, function(x) x$predictors_names) |> unique()
            vn <- fit_bru$names.fixed[grep('Intercept', fit_bru$names.fixed,invert = TRUE)]
            ii <- "Intercept"
          }
          assertthat::assert_that(all( vn %in% names(preds) ))
          preds <- subset(preds, select = vn )
          # Add offset if set
          if(!is.Waiver(model$offset)){
            ovn <- "spatial_offset"
            ofs <- paste0("", ovn," + ")
          } else { ofs <- ""}

          pfo <- stats::as.formula(
            paste0("~",fun,"( ",ii, " + ", ofs, paste0(vn, collapse = " + "),
                   # Add spatial latent effects
                   ifelse("latentspatial" %in% self$list_data(),
                          paste("+",paste0("spatial.field",1:length(model$biodiversity),collapse = " + ")),
                          ""),
                   ")")
          )
          # --- #
          cores <- if(getOption("ibis.runparallel")) getOption("ibis.nthread") else NULL
          # Make a prediction
          suppressWarnings(
            pred_bru <- inlabru:::predict.bru(
              object = fit_bru,
              num.threads = cores,
              newdata = preds,
              probs = c(0.05, 0.5, 0.95),
              formula = pfo,
              n.samples = 1000 # Pass as parameter?
            )
          )
          pred_bru$cv <- pred_bru$sd / pred_bru$mean
          # Convert to raster if not set already
          if(!is.Raster(pred_bru)) pred_bru <- terra::rast(pred_bru)

          # Get only the predicted variables of interest
          if(utils::packageVersion("inlabru") <= '2.5.2'){
            # Older version where probs are ignored
            prediction <- subset(pred_bru, c("mean","sd","q0.025", "median", "q0.975", "cv"))
            names(prediction) <- c("mean","sd","q0.025", "median", "q0.975", "cv")
          } else {
            prediction <- subset(pred_bru, c("mean","sd","q0.05", "q0.5", "q0.95", "cv") )
            names(prediction) <- c("mean", "sd", "q05", "q50", "q95", "cv")
          }

        } else {
          prediction <- NULL
        }

        # Compute end of computation time
        settings$set('end.time', Sys.time())

        # Definition of INLA Model object ----
        out <- bdproto(
          "INLABRU-Model",
          DistributionModel,
          id = model$id,
          model = model,
          settings = settings,
          fits = list(
            "fit_best" = fit_bru,
            "fit_best_equation" = self$get_data("components"),
            "mesh"     = self$get_data('mesh'),
            "spde"     = self$get_data('latentspatial'),
            "prediction" = prediction
          ),
          # Projection function
          project = function(self, newdata, form = NULL, n.samples = 1000, layer = "mean"){
            assertthat::assert_that('fit_best' %in% names(self$fits),
                                    is.data.frame(newdata) || is.matrix(newdata) || inherits(newdata,'SpatialPixelsDataFrame'),
                                    is.null(form) || is.character(form) || is.formula(form)
            )
            # Get model
            mod <- self$get_data('fit_best')
            model <- self$model
            # Also get settings for bias values
            settings <- self$settings

            # Clamp?
            if( settings$get("clamp") ) newdata <- clamp_predictions(model, newdata)

            # Set target variables to bias_value for prediction if specified
            if(!is.Waiver(settings$get('bias_variable'))){
              for(i in 1:length(settings$get('bias_variable'))){
                if(settings$get('bias_variable')[i] %notin% names(newdata)){
                  if(getOption('ibis.setupmessages')) myLog('[Estimation]','red','Did not find bias variable in prediction object!')
                  next()
                }
                newdata[[settings$get('bias_variable')[i]]] <- settings$get('bias_value')[i]
              }
            }

            # If newdata is not yet a SpatialPixel object, transform
            if(!inherits(newdata,'SpatialPixelsDataFrame')){
              assertthat::assert_that(
                assertthat::has_name(newdata,c('x','y'))
              )
              # Convert predictors to SpatialPixelsDataFrame as required for inlabru
              newdata <- sp::SpatialPointsDataFrame(coords = newdata[,c('x', 'y')],
                                                  data = newdata[, names(newdata) %notin% c('x','y')],
                                                  proj4string = self$get_data('mesh')$crs
              )
              newdata <- subset(newdata, stats::complete.cases(newdata@data)) # Remove missing data
              newdata <- methods::as(newdata, 'SpatialPixelsDataFrame')
            }
            # Check that model variables are in prediction dataset
            assertthat::assert_that(
              all(mod$names.fixed[grep('Intercept', mod$names.fixed,invert = TRUE)] %in% names(newdata))
                  )

            if(is.null(form)){
              # Try and guess backtransformation
              backtransf <- ifelse(mod$bru_info$lhoods[[1]]$family == 'poisson','exp','logistic')

              # Build the formula
              if(length(model$biodiversity)>1){
                vn <- lapply(model$biodiversity, function(x) x$predictors_names) |>  (\(.) do.call(c, .))() |> unique()
                ii <- paste("Intercept + ",paste0('Intercept_',
                                                  make.names(tolower(sapply( model$biodiversity, function(x) x$name ))),'_', # Make intercept from name
                                                  sapply( model$biodiversity, function(x) x$type ),collapse = ' + ')
                )
                # Assert that variables are used in the likelihoods
                assertthat::assert_that(
                  all( vn %in% (lapply(likelihoods, function(x) all.vars(x$formula) ) |>  (\(.) do.call(c, .))() |> unique() ) )
                )
              } else {
                vn <- sapply(model$biodiversity, function(x) x$predictors_names) |> unique()
                # vn <- mod$names.fixed[grep('Intercept', fit_bru$names.fixed,invert = TRUE)]
                assertthat::assert_that(all(vn %in% mod$names.fixed))
                ii <- "Intercept"
              }

              form <- stats::as.formula(
                paste0("~",backtransf,"( ",ii, " + ", paste0(vn, collapse = " + "),
                       ifelse(length(mod$summary.spde2.blc)>0, "+ spatial.field", ""),
                       ")")
              )
            }

            # Perform the projection
            suppressWarnings(
              out <- inlabru:::predict.bru(
                object = mod,
                data = newdata,
                formula = form,
                probs = c(0.05,0.5,0.95),
                n.samples = n.samples
              )
            )
            out$cv <- out$sd / out$mean
            # Get only the predicted variables of interest
            if(utils::packageVersion("inlabru") <= '2.5.2'){
              # Older version where probs are ignored
              out <- terra::rast(
                out[,c("mean","sd","q0.025", "median", "q0.975", "cv")]
              )
              names(out) <- c("mean","sd","q0.025", "median", "q0.975", "cv")
            } else {
              out <- terra::rast(
                out[,c("mean","sd","q0.05", "q0.5", "q0.95", "cv")]
              )
              names(out) <- c("mean", "sd", "q05", "q50", "q95", "cv")
            }
            # Return result
            return(out)
          },
          # Partial response
          partial = function(self, x.var, constant = NULL, variable_length = 100, values = NULL, plot = TRUE, type = "response"){
            # We use inlabru's functionalities to sample from the posterior
            # a given variable. A prediction is made over a generated fitted data.frame
            # Check that provided model exists and variable exist in model
            mod <- self$get_data('fit_best')
            model <- self$model
            df <- model$biodiversity[[1]]$predictors
            assertthat::assert_that(inherits(mod,'bru'),
                                    'model' %in% names(self),
                                    is.character(x.var),
                                    is.numeric(variable_length), variable_length >=1,
                                    is.null(constant) || is.numeric(constant),
                                    is.null(values) || is.numeric(values)
            )

            # Match variable name
            if(!is.null(mod$summary.random)) vn <- names(mod$summary.random) else vn <- ""
            x.var <- match.arg(x.var, c( mod$names.fixed, vn), several.ok = FALSE)

            # Make a prediction via inlabru
            if(any(model$predictors_types$type=="factor")){
              rr <- sapply(df[model$predictors_types$predictors[model$predictors_types$type=="numeric"]],
                           function(x) range(x, na.rm = TRUE)) |> as.data.frame()
            } else {
              rr <- sapply(df, function(x) range(x, na.rm = TRUE)) |> as.data.frame()
            }
            assertthat::assert_that(nrow(rr)>1, ncol(rr)>=1)

            df_partial <- list()
            # Set length out to value length to have equal coverage
            if(!is.null(values)){  variable_length <- length(values) }
            # Add all others as constant
            if(is.null(constant)){
              for(n in names(rr)) df_partial[[n]] <- rep( mean(df[[n]], na.rm = TRUE), variable_length )
            } else {
              for(n in names(rr)) df_partial[[n]] <- rep( constant, variable_length )
            }
            if(!is.null(values)){
              df_partial[[x.var]] <- values
            } else {
              df_partial[[x.var]] <- seq(rr[1,x.var], rr[2,x.var], length.out = variable_length)
            }
            df_partial <- df_partial |> as.data.frame()

            if(any(model$predictors_types$type=="factor")){
              lvl <- levels(model$predictors[[model$predictors_types$predictors[model$predictors_types$type=="factor"]]])
              df_partial[model$predictors_types$predictors[model$predictors_types$type=="factor"]] <-
                factor(lvl[1], levels = lvl)
            }

            ## plot the unique effect of the covariate
            fun <- ifelse(length(model$biodiversity) == 1 && model$biodiversity[[1]]$type == 'poipa', "logistic", "exp")
            pred_cov <- inlabru:::predict.bru(mod,
                                df_partial,
                                stats::as.formula( paste("~ ",fun,"(", paste(mod$names.fixed,collapse = " + ") ,")") ),
                                n.samples = 100,
                                probs = c(0.05,0.5,0.95)
                                )
            pred_cov$cv <- pred_cov$sd / pred_cov$mean

            o <- pred_cov
            names(o)[grep(x.var, names(o))] <- "partial_effect"
            if(utils::packageVersion("inlabru") <= '2.5.2'){
              # Older version where probs are ignored
              o <- subset(o, select = c("partial_effect", "mean", "sd", "median", "q0.025", "q0.975", "cv"))
              names(o) <- c("partial_effect", "mean", "sd", "median", "lower", "upper", "cv")
            } else {
              o <- subset(o, select = c("partial_effect", "mean", "sd", "q0.05", "q0.5", "q0.95", "cv"))
              names(o) <- c("partial_effect", "mean", "sd", "lower", "median", "upper", "cv")
            }

            # Do plot and return result
            if(plot){
              pm <- ggplot2::ggplot(data = o, ggplot2::aes(x = partial_effect, y = median,
                                                     ymin = lower,
                                                     ymax = upper) ) +
                ggplot2::theme_classic() +
                ggplot2::geom_ribbon(fill = "grey90") +
                ggplot2::geom_line() +
                ggplot2::labs(x = x.var, y = "Partial effect")
              print(pm)
            }
            return(o |> as.data.frame() )
          },
          # (S)partial effect
          spartial = function(self, x.var, constant = NULL, plot = TRUE, type = "response"){
            # We use inlabru's functionalities to sample from the posterior
            # a given variable. A prediction is made over a generated fitted data.frame
            # Check that provided model exists and variable exist in model
            mod <- self$get_data('fit_best')
            model <- self$model
            assertthat::assert_that(inherits(mod,'bru'),
                                    'model' %in% names(self),
                                    is.character(x.var),
                                    is.null(constant) || is.numeric(constant)
            )

            # Match variable name
            x.var <- match.arg(x.var, mod$names.fixed, several.ok = FALSE)

            # Convert predictors to SpatialPixelsDataFrame as required for inlabru
            df_partial <- sp::SpatialPointsDataFrame(coords = model$predictors[,c('x', 'y')],
                                                data = model$predictors[, names(model$predictors) %notin% c('x','y')],
                                                proj4string = self$get_data('mesh')$crs
            )
            df_partial <- subset(df_partial, stats::complete.cases(df_partial@data)) # Remove missing data
            df_partial <- methods::as(df_partial, 'SpatialPixelsDataFrame')

            # Add all others as constant
            if(is.null(constant)){
              for(n in names(df_partial)) if(n != x.var) df_partial[[n]] <- suppressWarnings( mean(model$predictors[[n]], na.rm = TRUE) )
            } else {
              for(n in names(df_partial)) if(n != x.var) df_partial[[n]] <- constant
            }
            if(any(model$predictors_types$type=="factor")){
              lvl <- levels(model$predictors[[model$predictors_types$predictors[model$predictors_types$type=="factor"]]])
              df_partial[[model$predictors_types$predictors[model$predictors_types$type=="factor"]]] <-
                factor(lvl[1], levels = lvl)
              # FIXME: Assigning the first level (usually reference) for now. But ideally find a way to skip factors from partial predictions
            }

            fun <- ifelse(length(model$biodiversity) == 1 && model$biodiversity[[1]]$type == 'poipa', "logistic", "exp")
            pred_cov <- inlabru:::predict.bru(mod,
                                              df_partial,
                                              stats::as.formula( paste("~ ",fun,"( Intercept + ", x.var ,")") ),
                                              n.samples = 100,
                                              probs = c(0.05,0.5,0.95)
            )
            pred_cov$cv <- pred_cov$sd / pred_cov$mean

            # Do plot and return result
            if(plot){
              o <- pred_cov
              ggplot2::ggplot() +
                ggplot2::theme_classic(base_size = 18) +
                inlabru::gg(o, ggplot2::aes(fill = mean)) +
                ggplot2::scale_fill_gradientn(colours = ibis_colours$divg_bluegreen) +
                ggplot2::labs(x = "", y = "", title = paste0("Spartial of ", x.var))
            }

            # Depending on the package version return out
            if(utils::packageVersion("inlabru") <= '2.5.2'){
              # Older version where probs are ignored
              return(
                terra::rast(
                  pred_cov[,c("mean","sd","q0.025", "median", "q0.975", "cv")] # Columns need to be adapted if quantiles are changed
                )
              )
            } else {
              return(
                terra::rast(
                  pred_cov[,c("mean","sd","q0.05", "q0.5", "q0.95", "cv")] # Columns need to be adapted if quantiles are changed
                )
              )
            }
          },
          # Residual function
          get_residuals = function(self){
            # Get best object
            obj <- self$get_data("fit_best")
            if(is.Waiver(obj)) return(obj)
            # Get residuals
            rd <- obj$residuals$deviance.residuals
            assertthat::assert_that(length(rd)>0)
            return(rd)
          },
          # Get coefficients
          get_coefficients = function(self){
            # Returns a vector of the coefficients with direction/importance
            cofs <- self$summary()
            cofs <- subset(cofs, select = c("variable", "mean", "sd"))
            names(cofs) <- c("Feature", "Beta", "Sigma")
            # Remove intercept(s)
            int <- grep("Intercept",cofs$Feature,ignore.case = TRUE)
            if(length(int)>0) cofs <- cofs[-int,]
            return(cofs)
          },
          # Function to plot SPDE if existing
          plot_spatial = function(self, spat = NULL, type = "response", what = "spatial.field1", ...){
            # Get mesh, domain and model
            mesh <- self$get_data("mesh")
            domain <- methods::as(self$model$background, "Spatial")
            mod <- self$get_data('fit_best')
            type <- match.arg(type, c("response", "predictor"), several.ok = FALSE)

            assertthat::assert_that(!is.null(mod$model.random),
                                    msg = "No spatial latent was estimated in the model!")

            if(mod$model.random == "SPDE2 model") {
              assertthat::assert_that(inherits(mod,'bru'),
                                      inherits(mesh, 'inla.mesh'),
                                      is.null(spat) || inherits("SpatialPixelsDataFrame"),
                                      'model' %in% names(self),
                                      is.character(what)
              )
              # Check whether random variable exists, otherwise raise warning
              if(!(what %in% names(mod$summary.random))){
                stop(paste0(
                  "Spatial random effect not found. Set 'what' to one of these: ",
                  paste0(names(mod$summary.random),collapse = " | ")
                ))
              }

              # Predict the spatial intensity surface
              if(is.null(spat)){
                spat <- inlabru::pixels(mesh, mask = domain)
              }
              # FIXME: Does not work for other link functions
              if(type == "response") fun <- 'exp' else fun <- ''

              suppressWarnings(
                lambda <- inlabru:::predict.bru(mod,
                                                spat,
                                                stats::as.formula(paste0("~ ",fun,"(",what," + Intercept)"))
                                                )
              )

              # Convert to raster stack
              lambda <- terra::rast(lambda)

              # Also get SPDE posteriors of the matern correlation and coveriance function
              corplot <- inlabru:::plot.prediction(inlabru::spde.posterior(mod, what, what = "matern.correlation")) +
                ggplot2::ggtitle("Matern correlation")
              covplot <- inlabru:::plot.prediction(inlabru::spde.posterior(mod, what, what = "matern.covariance")) +
                ggplot2::ggtitle("Matern covariance")
              inlabru::multiplot(covplot, corplot)

              return(lambda)
            } else {
              message("No SPDE effect found.")
            }
          }
        )
      }
    ))
}
