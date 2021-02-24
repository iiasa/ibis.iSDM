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
NULL

#' @export
BiodiversityDistribution <- bdproto(
  "BiodiversityDistribution",
  background    = new_waiver(),
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

    ex =  self$show_background_info()
    pn = ifelse(is.Waiver(self$get_predictor_names()),'None',name_atomic(self$get_predictor_names(), "predictors"))
    of = ifelse(is.Waiver(self$offset), '', paste0( "\n  offset:         <", self$get_offset(),">" ) )
    message(paste0('\033[1m','\033[36m','<',self$name(),'>','\033[39m','\033[22m',
                   "\nBackground extent: ",
                   "\n     xmin: ", ex[['extent']][1], ", xmax: ", ex[['extent']][2],",",
                   "\n     ymin: ", ex[['extent']][3], ", ymax: ", ex[['extent']][4],
                   "\n   projection: ", ex[['proj']],
                   "\n --------- ",
                   "\n",self$biodiversity$show(),
                   "\n --------- ",
                   "\n  predictors:     ", pn,
                   "\n  priors:         ", "Not yet implemented",
                   "\n  latent factors: ", paste(self$get_latent(),collapse = ', '),
                   of,
                   "\n  log:            ", "Not yet implemented",
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
    o[['extent']] <- round( st_bbox(r), 3)
    o[['proj']] <-  projection(r)
    return(o)
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
  set_latent = function(self, type, spatial_model = NULL){
    assertthat::assert_that(is.character(type),
                            type %in% c('<Spatial>','<Temporal>','<Spatial-temporal>'))
    # Assign argument if existing
    if(!is.null(spatial_model)){
      type <- paste0('<Spatial | ',spatial_model,'>')
      attr(type, 'spatial_model') <- spatial_model
    }
    if(!is.Waiver(self$latentfactors)){
      self$latentfactors <- unique(c(self$latentfactors,type))
    } else { self$latentfactors <- type }
  },
  # Get latent factors
  get_latent = function(self){
    if(is.Waiver(self$latentfactors)) return('None')
    self$latentfactors
  },
  # Set predictors
  set_predictors = function(self, x){
    assertthat::assert_that(inherits(x, "PredictorDataset"))
    self$predictors <- x
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
  # Set offset
  # FIXME: For logical consistency could define a new bdproto object
  set_offset = function(self, x){
    assertthat::assert_that(inherits(x, "Raster"))
    if (!is.Waiver(self$offset)) warning("Overwriting previously defined offset.")
    self$offset <- x
  },
  # Get offset (print name)
  get_offset = function(self){
    if(is.Waiver(self$offset)) return('None')
    names( self$offset )
  },
  # Remove predictors
  rm_predictors = function(self, names){
    assertthat::assert_that(
      is.character(names) || assertthat::is.scalar(names) || is.vector(names)
    )
    self$predictors$rm_data(names)
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
  }
)
