#' @include utils.R waiver.R bdproto.R bdproto-biodiversitydataset.R
NULL

#' @export
if (!methods::isClass("BiodiversityDistribution")) methods::setOldClass("BiodiversityDistribution")
NULL

#' Biodiversity Distribution master class
#'
#' @description
#' Base [`proto`] class for any biodiversity distribution objects.
#' Serves as container that supplies data and functions to
#' other [`proto`] classes.
#'
#' @details
#' Run [names()] on a [`distribution`] object to show all available functions.
#'
#' @name BiodiversityDistribution-class
#' @aliases BiodiversityDistribution
#' @family bdproto
#' @keywords bdproto
NULL

#' @export
BiodiversityDistribution <- bdproto(
  "BiodiversityDistribution",
  background    = new_waiver(),
  limits        = new_waiver(),
  biodiversity  = bdproto(NULL, BiodiversityDatasetCollection),
  predictors    = new_waiver(),
  priors        = new_waiver(),
  bias          = new_waiver(),
  latentfactors = new_waiver(),
  offset        = new_waiver(),
  log           = new_waiver(),
  engine        = new_waiver(),

  # Self printing function
  print = function(self) {
    # Query information from the distribution object
    ex <- self$show_background_info()
    pn <- ifelse(is.Waiver(self$get_predictor_names()),'None',name_atomic(self$get_predictor_names(), "predictors"))
    of <- ifelse(is.Waiver(self$offset), '', paste0( "\n  offset:         <", name_atomic(self$get_offset()),">" ) )
    pio <- ifelse(is.Waiver(self$priors), '<Default>', paste0('Priors specified (',self$priors$length(), ')') )
    bv <- ifelse(is.Waiver(self$bias), '', paste0( "\n  bias control:   <", self$bias$method, ">" ) )
    en <- ifelse(is.null(self$get_engine()), text_red("<NONE>"), self$get_engine() )

    message(paste0('\033[1m','\033[36m','<', self$name(),'>','\033[39m','\033[22m',
                   ifelse(is.Waiver(self$limits), "\nBackground extent: ", "\nBackground extent (limited): "),
                   "\n     xmin: ", ex[['extent']][1], ", xmax: ", ex[['extent']][2],",",
                   "\n     ymin: ", ex[['extent']][3], ", ymax: ", ex[['extent']][4],
                   "\n   projection: ", ex[['proj']],
                   "\n --------- ",
                   "\n", self$biodiversity$show(),
                   "\n --------- ",
                   "\n  predictors:     ", pn,
                   "\n  priors:         ", pio,
                   "\n  latent:         ", paste(self$get_latent(), collapse = ', '),
                   of,
                   bv,
                   "\n  log:            ", self$get_log(),
                   "\n  engine:         ", en
                   )
            )
  },
  # Print input messages
  show = function(self) {
    self$print()
  },
  # Self description
  name = function(self) {
    "Biodiversity distribution model"
  },
  # Get background stats
  show_background_info = function(self){
    assertthat::assert_that(inherits(self$background,'sf'))
    r <- self$background
    o <- list()
    o[['extent']] <- round( sf::st_bbox(r), 3)
    o[['proj']] <-  sf::st_crs(r)$proj4string
    return(o)
  },
  # Set limits
  set_limits = function(self, x, mcp_buffer = 0, limits_clip = FALSE){
    assertthat::assert_that(is.Raster(x) || inherits(x, "sf"),
                            msg = "Provide a SpatRaster or sf object!")
    # Construct limits object assuming zones
    x <- list(layer = x, limits_method = "zones",
                   "mcp_buffer" = mcp_buffer, "limits_clip" = limits_clip)

    bdproto(NULL, self, limits = x )
  },
  # Get provided limits
  get_limits = function(self){
    if(is.Waiver(self$limits)) return(new_waiver())
    return(self$limits)
  },
  # Remove limits
  rm_limits = function(self){
    bdproto(NULL, self, limits = new_waiver() )
  },
  # Function for querying predictor names if existing
  get_predictor_names = function(self) {
    if(is.Waiver(self$predictors)) return(self$predictors)
    if(inherits(self$predictors, "PredictorDataset")) {
      self$predictors$get_names()
    } else {
      stop("feature data is of an unrecognized class")
    }
  },
  # Adding latent factors
  set_latent = function(self, type, method = NULL, separate_spde = FALSE){
    assertthat::assert_that(is.character(type),
                            type %in% c('<Spatial>','<Temporal>','<Spatial-temporal>'),
                            is.character(method),
                            is.logical(separate_spde))
    # Assign argument if existing
    if(!is.null(method)){
      type <- paste0('<Spatial | ',method,'>')
      attr(type, 'method') <- method
      attr(type, 'separate_spde') <- separate_spde
    }
    if(!is.Waiver(self$latentfactors)){
      bdproto(NULL, self, latentfactors = type )
      # bdproto(NULL, self, latentfactors = unique(c(self$latentfactors,type)) )
    } else {
      bdproto(NULL, self, latentfactors = type )
    }
  },
  # Get latent factors
  get_latent = function(self){
    if(is.Waiver(self$latentfactors)) return('None')
    self$latentfactors
  },
  # Remove latent
  rm_latent = function(self){
    bdproto(NULL, self, latentfactors = new_waiver() )
  },
  # Get priors
  get_priors = function(self){
    return( self$priors )
  },
  # Set new priors
  set_priors = function(self, x ){
    assertthat::assert_that(inherits(x, 'PriorList'),
                            msg = 'An object created through `priors` has to be provided.')
    # Check if a priorlist is set. If yes, then combine the new one with existing priors
    if(is.Waiver(self$priors)){
      bdproto(NULL, self, priors = x )
    } else {
      # Get old prior list and combine with the new one
      pl <- self$priors
      pl$combine( x )
      bdproto(NULL, self, priors = pl )
    }
  },
  # Set biodiversity function
  set_biodiversity = function(self, id, p){
    assertthat::assert_that(inherits(self$biodiversity,'BiodiversityDatasetCollection'),
                            is.Id(id),
                            inherits(p, "BiodiversityDataset")
                            )
    # Get biodiversity dataset collection
    bdcol <- bdproto(NULL, self)
    # Set the object
    bdcol$biodiversity$set_data(id, p)
    return(bdcol)
  },
  # Set predictors
  set_predictors = function(self, x){
    assertthat::assert_that(inherits(x, "PredictorDataset"))
    bdproto(NULL, self, predictors = x)
  },
  # Set Engine
  set_engine = function(self, x) {
    assertthat::assert_that(inherits(x, "Engine"))
    if(!is.Waiver(self$engine)) warning("Overwriting previously defined engine.")
    bdproto(NULL, self, engine = x)
  },
  # Get Engine
  get_engine = function(self){
    if(is.Waiver(self$engine)) return(NULL)
    self$engine$show()
  },
  # Remove engine
  rm_engine = function(self){
    if(!is.Waiver(self$engine)){
      bdproto(NULL, self, engine = new_waiver())
    }
  },
  # Get prior variables
  get_prior_variables = function(self){
    if(is.Waiver(self$priors)) return(NULL)
    self$priors$varnames()
  },
  # Set offset
  set_offset = function(self, x){
    assertthat::assert_that(is.Raster(x))
    bdproto(NULL, self, offset = x )
  },
  # Get offset (print name)
  get_offset = function(self){
    if(is.Waiver(self$offset)) return( self$offset )
    names( self$offset )
  },
  # Remove offsets
  rm_offset = function(self, what = NULL){
    if(is.null(what)){
      bdproto(NULL, self, offset = new_waiver() )
    } else {
      of <- self$offset
      of <- terra::subset(of, -what)
      bdproto(NULL, self, offset = of )
    }
  },
  # Plot offset
  plot_offsets = function(self){
    if(is.Waiver(self$offset)) return( self$offset )
    if(terra::nlyr(self$offset)>1){
      of <- sum(self$offset, na.rm = TRUE)
      of <- terra::mask(of, self$background)
    } else {of <- self$offset}
    terra::plot(of, col = ibis_colours$viridis_orig, main = "Combined offset")
  },
  # Offset type
  get_offset_type = function(self){
    if(is.Waiver(self$offset)) return( self$offset )
    # Get attributes
    at <- list()
    at[['distance_function']] <- attr(self$offset, 'distance_function')
    at[['distance_max']] <- attr(self$offset, 'distance_function')
    if(!is.null(attr(self$offset,"logistic_coefficients")))
      at[['logistic_coefficients']] <- attr(self$offset, "logistic_coefficients")
    return(at)
  },
  # set_biascontrol
  set_biascontrol = function(self, x, method, value){
    assertthat::assert_that(missing(x) || is.Raster(x),
                            all(is.numeric(value)))
    if(missing(x)) {
      assertthat::assert_that(method == "proximity",
                              msg = paste0("Supply a layer for method ", method))
      x <- NULL
    }
    bdproto(NULL, self, bias = list(layer = x, method = method, bias_value = value) )
  },
  # Get bias control (print name)
  get_biascontrol = function(self){
    if(is.Waiver(self$bias)) return( self$bias )
    names( self$bias )
  },
  # Remove bias controls
  rm_biascontrol = function(self){
    bdproto(NULL, self, bias = new_waiver() )
  },
  # Plot bias variable
  plot_bias = function(self){
    if(is.Waiver(self$bias)) return( self$bias )
    terra::plot(self$bias$layer, col = ibis_colours$viridis_plasma, main = "Bias variable")
  },
  # Get log
  get_log = function(self){
    if(is.Waiver(self$log)){
      return('<Console>')
    } else {
      # Print filename
      self$log$get_filename()
    }
  },
  # Set log
  set_log = function(self, x){
    assertthat::assert_that(inherits(x, "Log"))
    bdproto(NULL, self, log = x)
  },
  # Get extent
  get_extent = function(self){
    # Calculate the extent from the background
    if(!is.Waiver(self$background)) terra::ext(self$background) else NULL
  },
  # Get dimensions of extent
  get_extent_dimensions = function(self){
    # Calculate the dimensions of the background
    if(!is.Waiver(self$background)) extent_dimensions(self$background) else NULL
  },
  # Get projection
  get_projection = function(self){
    assertthat::assert_that(inherits(self$background,'sf'))
    sf::st_crs(self$background)
  },
  # Get resolution
  get_resolution = function(self){
    if(!is.Waiver(self$predictors)){
      self$predictors$get_resolution()
    }
  },
  # Remove predictors
  rm_predictors = function(self, names){
    if(is.Waiver(self$predictors) || is.null(self$predictors)) return(NULL)
    if(missing(names)){
      names <- self$get_predictor_names() # Assume all names
    }
    assertthat::assert_that(
      is.character(names) || assertthat::is.scalar(names) || is.vector(names)
    )
    # Get predictor collection
    prcol <- bdproto(NULL, self)
    # Set the object
    prcol$predictors$rm_data(names)
    if(base::length(prcol$get_predictor_names())==0) prcol$predictors <- new_waiver()
    return(prcol)
  },
  # Remove priors
  rm_priors = function(self, names = NULL){
    assertthat::assert_that(is.null(names) || is.vector(names) || is.character(names))
    if(is.Waiver(self$priors)) {return( bdproto(NULL, self) )}
    # Delete selected priors
    if(is.null(names)){
      priors <- new_waiver()
    } else {
      priors <- self$priors
      ids <- priors$ids()[which( priors$varnames() %in% names)]
      for(id in ids) priors$rm(id)
    }
    bdproto(NULL, self, priors = priors )
  },
  # Show number of biodiversity records
  show_biodiversity_length = function(self){
    sum( self$biodiversity$length() )
  },
  # Get Equations of biodiversity records
  show_biodiversity_equations = function(self){
    self$biodiversity$show_equations()
  },
  # Get equations of biodiversity records
  get_biodiversity_equations = function(self){
    self$biodiversity$get_equations()
  },
  # Get biodiversity types
  get_biodiversity_types = function(self){
    self$biodiversity$get_types()
  },
  # Dummy function for plot and summary
  plot = function(self){
    message("No generic plotting implemented!")
  },
  summary = function(self){
    message("No generic summary function implemented! Try print.")
  }
)
