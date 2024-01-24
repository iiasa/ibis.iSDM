#' @include waiver.R
NULL

if (!methods::isClass("BiodiversityDatasetCollection")) methods::setOldClass("BiodiversityDatasetCollection")
if (!methods::isClass("BiodiversityDataset")) methods::setOldClass("BiodiversityDataset")

#---- BiodiversityDatasetCollection ----
#' BiodiversityDatasetCollection super class description
#'
#' @description
#' Acts a container for a specified set of BiodiversityDataset contained within.
#' Functions are provided to summarize across the BiodiversityDataset-class objects.
#'
#' @keywords classes
#'
#' @name BiodiversityDatasetCollection-class
NULL

#' @rdname BiodiversityDatasetCollection-class
#' @export
BiodiversityDatasetCollection <- R6::R6Class(
  "BiodiversityDatasetCollection",
  public = list(
    #' @field data A [`list`] of [`BiodiversityDataset-class`] objects.
    #' @field name The default name of this collection as [`character`].
    data = NULL,
    name = NULL,

    #' @description
    #' Initializes the object and creates an empty list
    #' @return NULL
    initialize = function(){
      self$data <- list()
      self$name <- 'Biodiversity data'
    },

    #' @description
    #' Print the names and properties of all Biodiversity datasets contained within
    #' @param format A [`logical`] flag on whether a message should be printed.
    #' @return A message on screen
    print = function(format = TRUE) {
      ty = if(is.Waiver(self$get_types())){
        ty = '\n   \033[31mNone\033[39m'} else { ty = paste0("\n   ",self$get_types()) }
      if(self$length()>0){ obs = paste0(" <",self$get_observations()," records>") } else obs = ''
      # FIXME: Prettify
      m <- paste0(self$name,':', paste0(ty,obs,collapse = ''))
      if(format) message(m) else return(m)
    },

    #' @description
    #' Aliases that calls print.
    #' @return A message on screen
    show = function(){
      self$print(format = FALSE)
    },

    #' @description
    #' Types of all biodiversity datasets included in this
    #' @param short A [`logical`] flag whether types should be in short format.
    #' @return A [`character`] vector.
    get_types = function(short = FALSE){
      if (base::length(self$data) > 0)
      {
        return(sapply(self$data, function(z) z$get_type( short ) ))
      } else return(new_waiver())
    },

    #' @description
    #' Get names and format them if necessary
    #' @param format A [`logical`] flag whether names are to be formatted
    #' @return A [`character`] vector.
    get_names = function(format = FALSE){
      x <- lapply(self$data, function(z) z$name)
      if(format) x <- make.names( tolower(x) )
      x
    },

    #' @description
    #' Add a new Biodiversity dataset to this collection.
    #' @param x A [`character`] with the name or id of this dataset.
    #' @param value A BiodiversityDataset
    #' @return Invisible
    set_data = function(x, value){
      assertthat::assert_that(assertthat::is.string(x),
                              inherits(value, "BiodiversityDataset"))
      self$data[[x]] <- value
      invisible(self)
    },

    #' @description
    #' Get a specific Biodiversity dataset by id
    #' @param id A [`character`] with a given id for the dataset.
    #' @return Returns a BiodiversityDataset.
    get_data_object = function(id) {
      assertthat::assert_that(is.Id(id) || is.character(id) )
      if (!id %in% names(self$data))
        return(new_waiver())
      return(self$data[[id]])
    },

    #' @description
    #' Get all biodiversity observations from a given dataset.
    #' @param id A [`character`] with a given id for the dataset.
    #' @return Returns all data from a set BiodiversityDataset.
    get_data = function(id){
      assertthat::assert_that(is.Id(id) || is.character(id))
      o <- self$get_data_object(id)
      o$get_data()
    },

    #' @description
    #' Get coordinates for a given biodiversity dataset. Else return a wkt object
    #' @param id A [`character`] with a given id for the dataset.
    #' @return All coordinates from a given object in `data.frame`.
    get_coordinates = function(id){
      assertthat::assert_that(is.Id(id) || is.character(id))
      # Get data
      o <- self$get_data(id)
      o <- guess_sf(o)
      # Add lowercase coordinates for consistency
      o$x <- sf::st_coordinates(o)[,1]
      o$y <- sf::st_coordinates(o)[,2]

      # Return coordinates
      if(utils::hasName(o,'geom')) sf::st_coordinates(o) else o[,c('x','y')]
    },

    #' @description
    #' Convenience function to mask all input datasets.
    #' @param mask A \code{SpatRaster} or `sf` object.
    #' @param inverse A `logical` flag if the inverse should be masked instead.
    #' @return Invisible
    mask = function(mask, inverse = FALSE){
      # Check whether prediction has been created
      biob <- self$data
      if(!is.Waiver(biob)){
        # If mask is SpatRaster, convert to polygons
        if(!inherits(mask, 'sf')){
          mask <- terra::as.polygons(mask) |> sf::st_as_sf()
        }
        # Now for each element in biob, mask
        for(id in names(biob)){
          biob[[id]]$mask(mask, inverse = inverse)
        }
        invisible(self)
      }
    },

    #' @description
    #' Remove a specific biodiversity dataset by id
    #' @param id A [`character`] with a given id for the dataset.
    #' @return Invisible
    rm_data = function(id) {
      assertthat::assert_that(is.Id(id) || is.character(id),
                              id %in% names(self$data) )
      self$data[[id]] <- NULL
      invisible(self)
    },

    #' @description
    #' Number of Biodiversity Datasets in connection
    #' @return A [`numeric`] with the number of datasets.
    length = function() {
      base::length(self$data)
    },

    #' @description
    #' Get number of observations of all datasets
    #' @return A [`numeric`] with the number of observations across datasets.
    get_observations = function() {
      x <- sapply(self$data, function(z) z$get_observations())
      x
    },

    #' @description
    #' Get equations from all datasets
    #' @return A [`list`] vector with all equations across datasets.
    get_equations = function(){
      x <- lapply(self$data, function(z) z$get_equation())
      x
    },

    #' @description
    #' Get families from datasets.
    #' @return A [`list`] vector with all families across datasets.
    get_families = function(){
      x <- lapply(self$data, function(z) z$get_family())
      x
    },

    #' @description
    #' Get custom link functions
    #' @return A [`list`] vector with all link functions across datasets.
    get_links = function(){
      x <- lapply(self$data, function(z) z$get_link() )
      x
    },

    #' @description
    #' Get fields with observation columns
    #' @return A [`list`] vector with the names of observation columns.
    get_columns_occ = function(){
      x <- lapply(self$data, function(z) z$get_column_occ())
      x
    },

    #' @description
    #' Get the weights across datasets.
    #' @return A [`list`] vector with the weights if set per dataset.
    get_weights = function(){
      x <- lapply(self$data, function(z) z$get_weight())
      x
    },

    #' @description
    #' Get ids of all assets in the collection.
    #' @return A [`list`] vector with the ids of all datasets.
    get_ids = function(){
      x <- lapply(self$data, function(z) z$get_id())
      x
    },

    #' @description
    #' Search for a specific biodiversity dataset with type
    #' @param type A [`character`] for a given data type.
    #' @return A [`character`] with the id(s) of datasets with the given type.
    get_id_byType = function(type){
      assertthat::assert_that(is.character(type), !missing(type))
      # Check whether type is correctly set
      if(type %notin% c("Point - Presence only","Point - Presence absence",
                        'Polygon - Presence only','Polygon - Presence absence')){
        type <- switch(type,
                       'poipo' = "Point - Presence only",
                       'poipa' = "Point - Presence absence",
                       'polpo' = 'Polygon - Presence only',
                       'polpa' = 'Polygon - Presence absence'
        )
      }
      if(is.null(type)) stop('Dataset type not found!')
      w <- which(self$get_types() %in% type)
      self$get_types()[w]
    },

    #' @description
    #' Get id by name
    #' @param name A [`character`] for a given name.
    #' @return A [`character`] with the id(s) of datasets with the given name.
    get_id_byName = function(name){
      assertthat::assert_that(is.character(name), !missing(name))
      # Get id(s) of dataset with given name
      r <- lapply(self$data, function(z) z$name)
      if(name %in% r){
        r[which(r==name)]
      } else character(0)
    },

    #' @description
    #' Show equations of all datasets
    #' @param msg A [`logical`] on whether to use print a message instead.
    #' @return Shows equations on screen or as [`character`].
    show_equations = function(msg = TRUE) {
      # Get equations
      x <- self$get_equations()

      # new names
      # n <- c(
      #     poipo = 'Point - Presence only',
      #     poipa = 'Point - Presence absence',
      #     polpo = 'Polygon - Presence only',
      #     polpa = 'Polygon - Presence absence'
      #   )
      # names(x) <- as.vector( n[match(names(x), names(n))] )
      # Prettify
      o <- paste0(names(x),":\n ",x,collapse = '\n')
      if(msg) message(o) else o
    },

    #' @description
    #' Plot the whole collection
    #' @note
    #' This can likely be beautified further.
    #' @return Invisible
    plot = function(){
      # Get observed columns
      cols <- self$get_columns_occ()

      par.ori <- graphics::par(no.readonly = TRUE)
      # Base plot
      g <- ggplot2::ggplot() + ggplot2::geom_sf() + ggplot2::labs( title = self$name)
      # Adding the other elements
      for(dataset in names(cols)){

        if('Polygon - Presence only' == self$get_types()[dataset] ) g <- g + ggplot2::geom_sf(data = st_as_sf(self$get_data(dataset))[cols[[dataset]]], fill = 'lightblue', alpha = .35 )

        if('Polygon - Presence absence' == self$get_types()[dataset] ){
          dd <- st_as_sf(self$get_data(dataset))[cols[[dataset]]]
          dd[[cols[[dataset]]]] <- factor(dd[[cols[[dataset]]]])
          g <- g + ggplot2::geom_sf(data = dd, fill = 'lightgreen', alpha = .35 )
        }

        if('Point - Presence only' == self$get_types()[dataset] ) g <- g + ggplot2::geom_sf(data = st_as_sf(self$get_data(dataset))[cols[[dataset]]],
                                                                                            colour = 'grey20', alpha = .5 )

        if('Point - Presence absence' == self$get_types()[dataset] ){
          dd <- st_as_sf(self$get_data(dataset))[cols[[dataset]]]
          dd[[cols[[dataset]]]] <- factor(dd[[cols[[dataset]]]])
          dd$observed <- dd[[cols[[dataset]]]]
          g <- g + ggplot2::geom_sf(data = dd, ggplot2::aes(colour = observed) )
        }
      }
      g
    }
  ),
  private = list(
    finalize = function() {

    }
  ),
  # No lock, thus allow members to be added.
  lock_objects = FALSE
)

#---- BiodiversityDataset ----
#' BiodiversityDataset prototype description
#'
#' @keywords classes
#'
#' @name BiodiversityDataset-class
NULL

#' @rdname BiodiversityDataset-class
#' @export
BiodiversityDataset <- R6::R6Class(
  "BiodiversityDataset",
  public = list(
    #' @field name The default name of this dataset as [`character`].
    #' @field id A [`character`] with the unique id for this dataset.
    #' @field equation A [`formula`] object containing the equation of how this dataset is modelled.
    #' @field family The family used for this dataset as [`character`].
    #' @field link The link function used for this data as [`character`].
    #' @field type A [`character`] with the type as [`character`].
    #' @field weight A [`numeric`] containing custom weights per observation for this dataset.
    #' @field field_occurrence A [`character`] with the name of the column name containing observations.
    #' @field data Contains the observational data in [`sf`] format.
    #' @field use_intercept A [`logical`] flag on whether intercepts are included for this dataset.
    #' @field pseudoabsence_settings Optionally provided pseudoabsence settings.
    name             = character(0),
    id               = character(0),
    equation         = new_waiver(),
    family           = character(0),
    link             = new_waiver(),
    type             = new_waiver(),
    weight           = new_waiver(),
    field_occurrence = character(0),
    data             = new_waiver(),
    use_intercept    = logical(0),
    pseudoabsence_settings = new_waiver(),

    #' @description
    #' Initializes the object and creates an empty list
    #' @param name The default name of this dataset as [`character`].
    #' @param id A [`character`] with the unique id for this dataset.
    #' @param equation A [`formula`] object containing the equation of how this dataset is modelled.
    #' @param family The family used for this dataset as [`character`].
    #' @param link The link function used for this data as [`character`].
    #' @param type A [`character`] with the type as [`character`].
    #' @param weight A [`numeric`] containing custom weights per observation for this dataset.
    #' @param field_occurrence A [`character`] with the name of the column name containing observations.
    #' @param data Contains the observational data in [`sf`] format.
    #' @param use_intercept A [`logical`] flag on whether intercepts are included for this dataset.
    #' @param pseudoabsence_settings Optionally provided pseudoabsence settings.
    #' @return NULL
    initialize = function(name, id, equation, family, link, type, weight,
                          field_occurrence, data, use_intercept,
                          pseudoabsence_settings){
      assertthat::assert_that(
        is.character(name) || is.null(name),
        is.character(id) || is.Id(id),
        (is.character(equation) || is.formula(equation)) || is.Waiver(equation),
        is.character(family) || is.function(family),
        is.character(link) || is.null(link),
        is.character(type),
        is.numeric(weight) || is.null(weight),
        is.character(field_occurrence),
        is.logical(use_intercept)
      )
      # Match type
      type <- match.arg(type, c("poipo", "poipa", "polpo", "polpa"), several.ok = FALSE)

      # If pseudo-absence missing, set default options
      if(missing(pseudoabsence_settings)){
        pseudoabsence_settings <- pseudoabs_settings()
      } else if(is.null(pseudoabsence_settings)) {
        pseudoabsence_settings <- pseudoabs_settings()
      }

      self$name <- name
      self$id <- id
      self$equation <- equation
      self$family <- family
      self$link <- link
      self$type <- type
      self$weight <- weight
      self$field_occurrence <- field_occurrence
      self$data <- data
      self$use_intercept <- use_intercept
      self$pseudoabsence_settings <- pseudoabsence_settings
    },

    #' @description
    #' Print the names and properties of all Biodiversity datasets contained within
    #' @return A message on screen
    print = function(){
      message(paste0('Biodiversity data:',
                     '\n Name:  ',self$name,
                     '\n Type:  ',self$get_type()
      ))
    },

    #' @description
    #' Set new equation and writes it into \code{formula}
    #' @param x A new [`formula`] object.
    #' @return Invisible
    set_equation = function(x){
      assertthat::assert_that(inherits(x, "formula"))
      self$formula <- x
      invisible(self)
    },

    #' @description
    #' Get equation
    #' @return A placeholder or [`formula`] object.
    get_equation = function(){
      if(is.Waiver(self$equation)) return('<Default>')
      self$equation
    },

    #' @description
    #' Function to print the equation
    #' @return A message on screen.
    show_equation = function(){
      if(!is.Waiver(equation) && !is.null(equation))
        message(equation)
      else message('None set. Default equation used (response ~ .)')
    },

    #' @description
    #' Get Id within the dataset
    #' @return A [`character`] with the id.
    get_id = function(){
      self$id
    },

    #' @description
    #' Get type of the dataset.
    #' @param short A [`logical`] flag if this should be formatted in shortform.
    #' @return A [`character`] with the type
    get_type = function(short = FALSE){
      if(short){
        self$type
      } else {
        switch (self$type,
                poipo = 'Point - Presence only',
                poipa = 'Point - Presence absence',
                polpo = 'Polygon - Presence only',
                polpa = 'Polygon - Presence absence'
        )
      }
    },

    #' @description
    #' Get field with occurrence information
    #' @return A [`character`] with the occurence field
    get_column_occ = function(){
      self$field_occurrence
    },

    #' @description
    #' Get family
    #' @return A [`character`] with the family for the dataset
    get_family = function(){
      assertthat::assert_that(is.character(self$family))
      self$family
    },

    #' @description
    #' Get custom link function
    #' @return A [`character`] with the family for the dataset
    get_link = function(){
      self$link
    },

    #' @description
    #' Get data from the object
    #' @return A [`sf`] object with the data
    get_data = function(){
      self$data
    },

    #' @description
    #' Get weight
    #' @return A [`numeric`] with the weights within the dataset.
    get_weight = function(){
      self$weight
    },

    #' @description
    #' Print input messages
    #' @return A message on screen.
    show = function(){
      self$print()
    },

    #' @description
    #' Collect info statistics about number of observations
    #' @return A [`numeric`] with the number of observations.
    get_observations = function() {
      nrow(self$data)
    },

    #' @description
    #' Convenience function to mask all input datasets.
    #' @param mask A \code{SpatRaster} or `sf` object.
    #' @param inverse A `logical` flag if the inverse should be masked instead.
    #' @param ... Any other parameters passed on to mask
    #' @return Invisible
    mask = function(mask, inverse = FALSE, ...){
      # Check whether prediction has been created
      biob <- self$data
      if(!is.Waiver(biob)){
        # If mask is SpatRaster, convert to polygons
        if(!inherits(mask, 'sf')){
          mask <- terra::as.polygons(mask) |> sf::st_as_sf()
        }
        # Now for each element in biob, mask
        biob <-
          suppressMessages(
            suppressWarnings(
              sf::st_crop(guess_sf(biob), mask)
            )
          )
        self$data <- biob
        invisible(self)
      }
    }
  ),

  # Any private entries
  private = list(
    finalize = function() {
    }
  ),
  # No lock, thus allow members to be added.
  lock_objects = FALSE
)
