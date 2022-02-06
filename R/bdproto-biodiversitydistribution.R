#' @include utils.R waiver.R bdproto.R bdproto-biodiversitydataset.R
NULL

#' @export
if (!methods::isClass("BiodiversityDistribution")) methods::setOldClass("BiodiversityDistribution")
NULL

#' Biodiversity Distribution master class
#'
#' Base [`proto`] class for any biodiversity distribution objects.
#' Serves as container that supplies data and functions to
#' other [`proto`] classes.
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
  latentfactors = new_waiver(),
  offset        = new_waiver(),
  log           = new_waiver(),
  engine        = new_waiver(),

  # Self printing function
  print = function(self) {
    # TODO: Prettify below
    # Query information from the distribution object
    ex =  self$show_background_info()
    pn = ifelse(is.Waiver(self$get_predictor_names()),'None',name_atomic(self$get_predictor_names(), "predictors"))
    of = ifelse(is.Waiver(self$offset), '', paste0( "\n  offset:         <", name_atomic(self$get_offset()),">" ) )
    pio = ifelse(is.Waiver(self$priors), '<Default>', paste0('Priors specified (',self$priors$length(), ')') )

    message(paste0('\033[1m','\033[36m','<',self$name(),'>','\033[39m','\033[22m',
                   ifelse(is.Waiver(self$limits),"\nBackground extent: ","\nBackground extent (limited): "),
                   "\n     xmin: ", ex[['extent']][1], ", xmax: ", ex[['extent']][2],",",
                   "\n     ymin: ", ex[['extent']][3], ", ymax: ", ex[['extent']][4],
                   "\n   projection: ", ex[['proj']],
                   "\n --------- ",
                   "\n",self$biodiversity$show(),
                   "\n --------- ",
                   "\n  predictors:     ", pn,
                   "\n  priors:         ", pio,
                   "\n  latent:         ", paste(self$get_latent(),collapse = ', '),
                   of,
                   "\n  log:            ", self$get_log(),
                   "\n  engine:         ", self$get_engine()
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
    o[['proj']] <-  raster::projection(r)
    return(o)
  },
  # Get provided limits
  get_limits = function(self){
    if(is.Waiver(self$limits)) return(NULL)
    return(self$limits)
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
  # Get priors
  get_priors = function(self){
    return( self$priors )
  },
  # Set new priors
  set_priors = function(self, x ){
    assertthat::assert_that(inherits(x, 'PriorList'),msg = 'An object created through `priors` has to be provided.')
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
    if (!is.Waiver(self$engine)) warning("Overwriting previously defined engine.")
    bdproto(NULL, self, engine = x)
  },
  # Get Engine
  get_engine = function(self){
    if(is.Waiver(self$engine)) return('None')
    self$engine$show()
  },
  # Get prior variables
  get_prior_variables = function(self){
    if(is.Waiver(self$priors)) return(NULL)
    self$priors$varnames()
  },
  # Set offset
  # FIXME: For logical consistency could define a new bdproto object
  set_offset = function(self, x){
    assertthat::assert_that(is.Raster(x))
    bdproto(NULL, self, offset = x )
  },
  # Get offset (print name)
  get_offset = function(self){
    if(is.Waiver(self$offset)) return( self$offset() )
    names( self$offset )
  },
  # Plot offset
  plot_offset = function(self){
    if(is.Waiver(self$offset)) return( self$offset() )
    if(raster::nlayers(self$offset)>1){
      of <- sum(self$offset,na.rm = TRUE)
      of <- raster::mask(of, self$background)
    } else {of <- self$offset}
    raster::plot(of, col = ibis_colours$viridis_orig, main = "Combined offset")
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
    if(!is.Waiver(self$background)) raster::extent(self$background) else NULL
  },
  # Get dimensions of extent
  get_extent_dimensions = function(self){
    # Calculate the dimensions of the background
    if(!is.Waiver(self$background)) extent_dimensions(self$background) else NULL
  },
  # Remove predictors
  rm_predictors = function(self, names){
    assertthat::assert_that(
      is.character(names) || assertthat::is.scalar(names) || is.vector(names)
    )
    # Get predictor collection
    prcol <- bdproto(NULL, self)
    # Set the object
    prcol$predictors$rm_data(names)
    return(prcol)
  },
  # Remove priors
  rm_priors = function(self, names = NULL){
    assertthat::assert_that(is.null(names) || is.vector(names) || is.character(names))
    if(is.Waiver(self$priors)) {return(NULL)}
    # Delete selected priors
    if(is.null(names)){
      self$priors <- new_waiver()
    } else {
      priors <- self$priors
      ids <- priors$ids()[which( priors$varnames() %in% names)]
      for(id in ids) priors$rm(id)
      self$priors <- priors
    }
    invisible()
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
  }
)
