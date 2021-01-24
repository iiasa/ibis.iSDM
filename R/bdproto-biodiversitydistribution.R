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
  equation      = new_waiver(),
  background    = new_waiver(),
  biodiversity  = bdproto(NULL, BiodiversityDatasetCollection),
  predictors    = new_waiver(),
  priors        = new_waiver(),
  latentfactors = new_waiver(),
  log           = new_waiver(),
  engine        = new_waiver(),

  # Self printing function
  print = function(self) {
    # r <- vapply(list(self$objective, self$targets), function(x) {
    #   if (is.Waiver(x))
    #     return("none")
    #   return(x$name())
    # }, character(1))
    #
    # d <- vapply(list(self$solver, self$decisions, self$portfolio), function(x) {
    #   if (is.Waiver(x))
    #     return("default")
    #   return(x$name())
    # }, character(1))
    # FIXME: Prettify below

    ex =  self$show_background()
    pn = ifelse(is.Waiver(self$predictor_names()),'None',name_atomic(self$predictor_names(), "predictors"))
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
  show_background = function(self){
    assertthat::assert_that(inherits(self$background,'Raster'))
    r <- self$background
    o <- list()
    o[['extent']] <- round( as.vector(raster::extent(r)), 3)
    o[['proj']] <-  projection(r)
    return(o)
  },
  # Function to print the equation
  show_equation = function(self){
    if(!is.Waiver(equation) && !is.null(equation))
      message(equation)
    else message('None set. Default equation used (response ~ .)')
  },
  # Function for querying predictor names if existing
  predictor_names = function(self) {
    if(is.Waiver(self$predictors)) return(self$predictors)
    if(inherits(self$predictors, "PredictorDataset")) {
      self$predictors$get_names()
    } else {
      stop("feature data is of an unrecognized class")
    }
  },
  # Adding latent factors
  set_latent = function(self, type){
    assertthat::assert_that(is.character(type),
                            type %in% c('<Spatial>','<Temporal>','<Spatial-temporal>'))
    if(!is.Waiver(self$latentfactors)){
      self$latentfactors <- unique(c(self$latentfactors,type))
    } else { self$latentfactors <- type}
  },
  # Get latent factors
  get_latent = function(self){
    if(is.Waiver(self$latentfactors)) return('None')
    self$latentfactors
  },
  # Set Engine
  set_engine = function(self, x) {
    assertthat::assert_that(inherits(x, "Engine"))
    if (!is.Waiver(self$engine)) warning("overwriting previously defined solver")
    bdproto(NULL, self, engine = x)
  },
  # Get Engine
  get_engine = function(self){
    if(is.Waiver(self$engine)) return('None')
    self$engine$show()
  }
)
