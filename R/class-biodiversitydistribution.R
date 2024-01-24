#' @include waiver.R class-biodiversitydataset.R
NULL

if (!methods::isClass("BiodiversityDistribution")) methods::setOldClass("BiodiversityDistribution")

#' Biodiversity Distribution master class
#'
#' @description Base [`R6`] class for any biodiversity distribution objects.
#' Serves as container that supplies data and functions to other [`R6`]
#' classes. Generally stores all objects and parameters added to a model.
#'
#' @details Run [`names()`] on a [`distribution`] object to show all available
#' functions.
#' @examples
#' # Query available functions and entries
#' background <- terra::rast(system.file('extdata/europegrid_50km.tif',
#' package='ibis.iSDM',mustWork = TRUE))
#' # Define model
#' x <- distribution(background)
#' names(x)
#'
#' @keywords classes
#'
#' @name BiodiversityDistribution-class
NULL

#' @rdname BiodiversityDistribution-class
#' @export
BiodiversityDistribution <- R6::R6Class(
  "BiodiversityDistribution",
  public = list(
    #' @field background A [`SpatRaster`] or [`sf`] object delineating the modelling extent.
    #' @field limits An optional [`sf`] object on potential extrapolation limits
    #' @field biodiversity A [`BiodiversityDatasetCollection-class`] object.
    #' @field predictors A [`PredictorDataset-class`] object.
    #' @field priors An optional [`PriorList`] object.
    #' @field control An optional Control object.
    #' @field latentfactors A [`character`] on whether latentfactors are used.
    #' @field offset A [`character`] on whether methods are used.
    #' @field log An optional [`Log-class`] object.
    #' @field engine A [`Engine-class`] object.
    background    = new_waiver(),
    limits        = new_waiver(),
    biodiversity  = new_waiver(),
    predictors    = new_waiver(),
    priors        = new_waiver(),
    control       = new_waiver(),
    latentfactors = new_waiver(),
    offset        = new_waiver(),
    log           = new_waiver(),
    engine        = new_waiver(),

    #' @description
    #' Initializes the object and creates an BiodiversityDataset by default.
    #' @param background A [`SpatRaster`] or [`sf`] object delineating the modelling extent.
    #' @param limits An optional [`sf`] object on potential extrapolation limits
    #' @param biodiversity A [`BiodiversityDatasetCollection-class`] object.
    #' @param ... Any other objects
    #' @return NULL
    initialize = function(background, limits, biodiversity, ...){
      assertthat::assert_that(
        is.Raster(background) ||inherits(background, "sf"),
        is.null(limits) || is.list(limits)
      )
      self$background <- background
      self$limits <- limits
      self$biodiversity <- biodiversity
      # Get Dots and save too
      dots <- list(...)
      for(el in names(dots)){
        self[[el]] <- dots[[el]]
      }
    },

    #' @description
    #' Looks for and returns the properties of all contained objects.
    #' @return A message on screen
    print = function() {
      # Query information from the distribution object
      ex <- self$show_background_info()
      pn <- ifelse(is.Waiver(self$get_predictor_names()),'None',
                   name_atomic(self$get_predictor_names(), "predictors"))
      of <- ifelse(is.Waiver(self$offset), '',
                   paste0( "\n  offset:         <", name_atomic(self$get_offset()),">" ) )
      pio <- ifelse(is.Waiver(self$priors),
                    '<Default>', paste0('Priors specified (',self$priors$length(), ')') )
      bv <- ifelse(is.Waiver(self$control), '',
                   paste0( "\n  control:        <", name_atomic(
                     paste0( self$control$type, " - ", self$control$method)
                   ), ">" ) )
      li <- ifelse(is.Waiver(self$limits), '',
                   paste0( "\n  limits:         <",paste0( self$limits$limits_method,collapse = ", "), ">" ))
      en <- ifelse(is.null(self$get_engine()),
                   text_red("<NONE>"), self$get_engine() )

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
                     li,
                     "\n  log:            ", self$get_log(),
                     "\n  engine:         ", en
      )
      )
    },

    #' @description
    #' An alias for print
    #' @return A message on screen
    show = function() {
      self$print()
    },

    #' @description
    #' Returns self-describing name
    #' @return A [`character`] with the name
    name = function() {
      "Biodiversity distribution model"
    },

    #' @description
    #' Summarizes extent and projection from set background
    #' @return A [`character`] with the name
    show_background_info = function(){
      assertthat::assert_that(inherits(self$background,'sf'))
      r <- self$background
      o <- list()
      o[['extent']] <- round( sf::st_bbox(r), 3)
      o[['proj']] <-  sf::st_crs(r)$proj4string
      return(o)
    },

    #' @description
    #' Specify new limits to the background
    #' @param x A [`list`] object with method and limit type.
    #' @seealso [add_control_extrapolation()]
    #' @return This object.
    set_limits = function(x){
      # Specify list
      assertthat::assert_that(is.list(x),
                              msg = "Provide a prepared list for the limits!")
      assertthat::assert_that(
        utils::hasName(x, "layer"), utils::hasName(x, "limits_method")
      )
      self$limits <- x
      return(self)
    },

    #' @description
    #' Get provided limits if set or a waiver
    #' @return A [`list`] or waiver.
    get_limits = function(){
      if(is.Waiver(self$limits)) return(new_waiver())
      return(self$limits)
    },

    #' @description
    #' Remove limits if set.
    #' @return This object.
    rm_limits = function(){
      self$limits <- new_waiver()
      return(self)
    },

    #' @description
    #' Function for querying predictor names if existing
    #' @return A [`character`] vector.
    get_predictor_names = function() {
      if(is.Waiver(self$predictors)) return(self$predictors)
      if(inherits(self$predictors, "PredictorDataset")) {
        self$predictors$get_names()
      } else {
        stop("feature data is of an unrecognized class")
      }
    },

    #' @description
    #' Adding latent factors to the object.
    #' @param type A [`character`] with the given type.
    #' @param method A [`character`] with a method.
    #' @param separate_spde A [`logical`] flag whether duplicate of SPDE effects are to be created.
    #' @seealso [add_latent_spatial()]
    #' @return This object.
    set_latent = function(type, method = NULL, separate_spde = FALSE){
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
      self$latentfactors <- type
      return(self)
    },

    #' @description
    #' Get latent factors if found in object.
    #' @return A [`character`] with those objects.
    get_latent = function(){
      if(is.Waiver(self$latentfactors)) return('None')
      self$latentfactors
    },

    #' @description
    #' Remove latent factors if found in object.
    #' @return This object.
    rm_latent = function(){
      self$latentfactors <- new_waiver()
      return(self)
    },

    #' @description
    #' Get prior object if found in object.
    #' @return This object.
    get_priors = function(){
      return( self$priors )
    },

    #' @description
    #' Specify new prior object. Overwrites existing ones
    #' @param x A [`PriorList-class`] object.
    #' @seealso [add_priors()]
    #' @return This object.
    set_priors = function(x){
      assertthat::assert_that(inherits(x, 'PriorList'),
                              msg = 'An object created through `priors` has to be provided.')
      # Check if a priorlist is set. If yes, then combine the new one with existing priors
      if(is.Waiver(self$priors)){
        self$priors <- x
      } else {
        # Get old prior list and combine with the new one
        pl <- self$priors
        pl$combine( x )
        self$priors <- pl
      }
      return(self)
    },

    #' @description
    #' Adds a new biodiversity object to the existing empty collection.
    #' @param id A [`character`] or id defining this object.
    #' @param p A [`BiodiversityDataset-class`] object.
    #' @seealso [add_biodiversity_poipa()], [add_biodiversity_poipo()], [add_biodiversity_polpa()], [add_biodiversity_polpo()]
    #' @return This object.
    set_biodiversity = function(id, p){
      assertthat::assert_that(inherits(self$biodiversity,'BiodiversityDatasetCollection'),
                              is.Id(id) || is.character(id),
                              inherits(p, "BiodiversityDataset")
      )
      # Get biodiversity dataset collection
      bdcol <- self
      # Set the object
      bdcol$biodiversity$set_data(id, p)
      return(bdcol)
    },

    #' @description
    #' Set a new Predictor object to this object.
    #' @param x A [`PredictorDataset-class`] with predictors for this object.
    #' @seealso [add_predictors()]
    #' @return This object.
    set_predictors = function(x){
      assertthat::assert_that(inherits(x, "PredictorDataset"))
      self$predictors <- x
      return(self)
    },

    #' @description
    #' Set a new Engine object to this object.
    #' @param x A [`Engine-class`] for this object.
    #' @return This object.
    set_engine = function(x) {
      assertthat::assert_that(inherits(x, "Engine"))
      if(!is.Waiver(self$engine)) warning("Overwriting previously defined engine.")
      self$engine <- x
      return(self)
    },

    #' @description
    #' Gets the name of the current engine if set.
    #' @return A [`character`] with the engine name
    get_engine = function(){
      if(is.Waiver(self$engine)) return(NULL)
      self$engine$show()
    },

    #' @description
    #' Removes the current engine if set.
    #' @return This object
    rm_engine = function(){
      if(!is.Waiver(self$engine)){
        self$engine <- new_waiver()
      }
      return( self )
    },

    #' @description
    #' Get prior variables
    #' @return A [`character`] with the variable names for which priors have been added.
    get_prior_variables = function(){
      if(is.Waiver(self$priors)) return(NULL)
      self$priors$varnames()
    },

    #' @description
    #' Specify new offsets.
    #' @param x A new [`SpatRaster`] object to be used as offset.
    #' @seealso [add_offset()]
    #' @return This object.
    set_offset = function(x){
      assertthat::assert_that(is.Raster(x))
      self$offset <- x
      return(self)
    },

    #' @description
    #' Get offset (print name)
    #' @return A [`character`] with all the offsets in here.
    get_offset = function(){
      if(is.Waiver(self$offset)) return( self$offset )
      names( self$offset )
    },

    #' @description
    #' Remove offsets if found.
    #' @param what Optional [`character`] of specific offsets to remove.
    #' @return This object.
    rm_offset = function(what = NULL){
      assertthat::assert_that(
        is.null(what) || is.character(what)
      )
      if(is.null(what)){
        self$offset <- new_waiver()
      } else {
        of <- self$offset
        of <- terra::subset(of, -what)
        self$offset <- of
      }
      return(self)
    },

    #' @description
    #' Plot offset if found.
    #' @return A graphical element.
    plot_offsets = function(){
      if(is.Waiver(self$offset)) return( self$offset )
      if(terra::nlyr(self$offset)>1){
        of <- sum(self$offset, na.rm = TRUE)
        of <- terra::mask(of, self$background)
      } else {of <- self$offset}
      terra::plot(of, col = ibis_colours$viridis_orig, main = "Combined offset")
    },

    #' @description
    #' Get offset parameters if found
    #' @return A [`list`] with the offset parameters if found.
    get_offset_type = function(){
      if(is.Waiver(self$offset)) return( self$offset )
      # Get attributes
      at <- list()
      at[['distance_function']] <- attr(self$offset, 'distance_function')
      at[['distance_max']] <- attr(self$offset, 'distance_function')
      if(!is.null(attr(self$offset,"logistic_coefficients")))
        at[['logistic_coefficients']] <- attr(self$offset, "logistic_coefficients")
      return(at)
    },

    #' @description
    #' Set new bias control
    #' @param type A [`character`] with the type of control object.
    #' @param x A new bias control object. Expecting a [`SpatRaster`] object.
    #' @param method The method used to create the object.
    #' @param value A bias value as [`numeric`].
    #' @return This object.
    set_control = function(type = "bias", x, method, value){
      assertthat::assert_that(missing(x) || is.Raster(x),
                              all(is.numeric(value)))
      # Check type of control
      type <- match.arg(type, c("bias"), several.ok = FALSE)
      if(type == "bias"){
        if(missing(x)) {
          assertthat::assert_that(method == "proximity",
                                  msg = paste0("Supply a layer for method ", method))
          x <- NULL
        }
        self$control <- list(type = type, layer = x,
                             method = method, bias_value = value)
        return(self)
      }
    },

    #' @description
    #' Get bias control (print name)
    #' @param type A [`character`] with the type of control object.
    #' @return A [`character`] with the bias object if found.
    get_control = function(type = "bias"){
      # Check type of control
      type <- match.arg(type, c("bias"), several.ok = FALSE)
      control <- self$control
      if(is.Waiver(control)) return( control )
      if(control$type == "bias" && type == "bias") return( control )
    },

    #' @description
    #' Remove bias controls if found.
    #' @return This object.
    rm_control = function(){
      self$control <- new_waiver()
      return(self)
    },

    #' @description
    #' Plot bias variable if set.
    #' @return A graphical element.
    plot_bias = function(){
      if(is.Waiver(self$control)) return( self$control )
      control <- self$control
      if(control$type == "bias"){
        terra::plot(control$layer,
                    col = ibis_colours$viridis_plasma,
                    main = "Bias variable")
      }
    },

    #' @description
    #' Returns the output filename of the current log object if set.
    #' @return A [`character`] where the output is returned.
    get_log = function(){
      if(is.Waiver(self$log)){
        return('<Console>')
      } else {
        # Print filename
        self$log$get_filename()
      }
    },

    #' @description
    #' Set a new log object
    #' @param x A [`Log-class`] object.
    #' @return This object
    set_log = function(x){
      assertthat::assert_that(inherits(x, "Log"))
      self$log <- x
      return(self)
    },

    #' @description
    #' Get extent
    #' @return Background extent or NULL.
    get_extent = function(){
      # Calculate the extent from the background
      if(!is.Waiver(self$background)) terra::ext(self$background) else NULL
    },

    #' @description
    #' Get dimensions of extent from background
    #' @return Background extent or NULL.
    get_extent_dimensions = function(){
      # Calculate the dimensions of the background
      if(!is.Waiver(self$background)) extent_dimensions(self$background) else NULL
    },

    #' @description
    #' Get projection from the background in crs format.
    #' @return A [`character`] of the projection
    get_projection = function(){
      assertthat::assert_that(inherits(self$background,'sf'))
      sf::st_crs(self$background)
    },

    #' @description
    #' Return resolution of the background object.
    #' @return A [`vector`] with the resolution.
    get_resolution = function(){
      if(!is.Waiver(self$predictors)){
        self$predictors$get_resolution()
      }
    },

    #' @description
    #' Remove predictiors. Either all of them or specific ones.
    #' @param names A [`character`] with the predictors to be removed.
    #' @return This object.
    rm_predictors = function(names){
      if(is.Waiver(self$predictors) || is.null(self$predictors)) return(NULL)
      if(missing(names)){
        names <- self$get_predictor_names() # Assume all names
      }
      assertthat::assert_that(
        is.character(names) || assertthat::is.scalar(names) || is.vector(names)
      )
      # Get predictor collection
      prcol <- self
      # Set the object
      prcol$predictors$rm_data(names)
      if(base::length(prcol$get_predictor_names())==0) prcol$predictors <- new_waiver()
      return(prcol)
    },

    #' @description
    #' Remove priors. Either all of them or specific ones.
    #' @param names A [`character`] with the priors to be removed.
    #' @return This object.
    rm_priors = function(names = NULL){
      assertthat::assert_that(is.null(names) || is.vector(names) || is.character(names))
      if(is.Waiver(self$priors)) return( self )

      # Delete selected priors
      if(is.null(names)){
        priors <- new_waiver()
      } else {
        priors <- self$priors
        ids <- priors$ids()[which( priors$varnames() %in% names)]
        for(id in ids) priors$rm(id)
      }
      self$priors <- priors
      return(self)
    },

    #' @description
    #' Show number of biodiversity records
    #' @return A [`numeric`] with sum of biodiversity records
    show_biodiversity_length = function(){
      sum( self$biodiversity$length() )
    },

    #' @description
    #' Show Equations of biodiversity records
    #' @return A message on screen.
    show_biodiversity_equations = function(){
      self$biodiversity$show_equations()
    },

    #' @description
    #' Get equations of biodiversity records
    #' @return A [`list`] vector.
    get_biodiversity_equations = function(){
      self$biodiversity$get_equations()
    },

    #' @description
    #' Query all biodiversity types in this object
    #' @return A [`character`] vector.
    get_biodiversity_types = function(){
      self$biodiversity$get_types()
    },

    #' @description
    #' Plots the content of this class.
    #' @note
    #' Not implemented yet.
    #' @return A message.
    plot = function(){
      message("No generic plotting implemented!")
    },

    #' @description
    #' Summary function for this object.
    #' @note
    #' Not implemented yet.
    #' @return A message.
    summary = function(){
      message("No generic summary function implemented! Try print.")
    }
  ),

  # Private counters
  private = list(
    finalize = function() {
    }
  ),
  # Don't lock objects
  lock_objects = FALSE
)
