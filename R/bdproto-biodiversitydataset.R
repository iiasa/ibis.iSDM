#' @include utils.R waiver.R bdproto.R
NULL

#' @export
if (!methods::isClass("BiodiversityDatasetCollection")) methods::setOldClass("BiodiversityDatasetCollection")
if (!methods::isClass("BiodiversityDataset")) methods::setOldClass("BiodiversityDataset")
NULL

#' BiodiversityDatasetCollection super prototype description
#'
#' Acts a container for BiodiversityDataset within
#'
#' @name BiodiversityDatasetCollection-class
#' @keywords bdproto
#' @family bdproto
#' @aliases BiodiversityDatasetCollection
NULL

#' @export
BiodiversityDatasetCollection <- bdproto(
  "BiodiversityDatasetCollection",
  data = list(),
  # Print the names of all Biodiversity datasets
  print = function(self) {
    message(self$show())
  },
  # Print a summary of all datasets
  show = function(self){
    ty = if(is.Waiver(self$get_types())){
      ty = '\n   \033[31mNone\033[39m'} else { ty = paste0("\n   ",self$get_types()) }
    if(self$length()>0){ obs = paste0(" <",self$get_observations()," records>") } else obs = ''
    # FIXME: Prettify
    paste0(self$name(),':',
                   paste0(ty,obs,collapse = '')
           )
  },
  # Name of this object
  name = function(self){
    'Biodiversity data'
  },
  # Types of all biodiversity datasets
  get_types = function(self, short = FALSE){
    if (base::length(self$data) > 0)
    {
      return(sapply(self$data, function(z) z$get_type( short ) ))
    } else return(new_waiver())
  },
  # Get names. Format if necessary
  get_names = function(self, format = FALSE){
    x <- lapply(self$data, function(z) z$name)
    if(format) x <- make.names( tolower(x) )
    x
  },
  # Add a new Biodiversity dataset to this collection
  set_data = function(self, x, value){
    assertthat::assert_that(assertthat::is.string(x),
                            inherits(value, "BiodiversityDataset"))
    self$data[[x]] <- value
    invisible()
  },
  # Get a specific Biodiversity dataset by name
  get_data_object = function(self, id) {
    assertthat::assert_that(is.Id(id) || is.character(id) )
    if (!id %in% names(self$data))
      return(new_waiver())
    return(self$data[[id]])
  },
  # Get biodiversity observations
  get_data = function(self, id){
    assertthat::assert_that(is.Id(id) || is.character(id))
    o <- self$get_data_object(id)
    o$get_data()
  },
  # Get coordinates for a given biodiversity dataset. Else return a wkt object
  get_coordinates = function(self, id){
    assertthat::assert_that(is.Id(id) || is.character(id))
    # Get data
    o <- self$get_data(id)
    o <- guess_sf(o)
    # Add lowercase coordinates for consistency
    o$x <- sf::st_coordinates(o)[,1]
    o$y <- sf::st_coordinates(o)[,2]

    # Return coordinates
    if(hasName(o,'geom')) sf::st_coordinates(o) else o[,c('x','y')]
  },
  # Remove a specific biodiversity dataset by id
  rm_data = function(self, id) {
    assertthat::assert_that(is.Id(id) || is.character(id),
                            id %in% names(self$data) )
    self$data[[id]] <- NULL
    invisible()
  },
  # Number of Biodiversity Datasets in connection
  length = function(self) {
    base::length(self$data)
  },
  # Get number of observations of all datasets
  get_observations = function(self) {
    x <- sapply(self$data, function(z) z$get_observations())
    x
  },
  # Get equations
  get_equations = function(self){
    x <- lapply(self$data, function(z) z$get_equation())
    x
  },
  # Get family
  get_families = function(self){
    x <- lapply(self$data, function(z) z$get_family())
    x
  },
  # Get custom link functions
  get_links = function(self){
    x <- lapply(self$data, function(z) z$get_link() )
    x
  },
  # Get fields with observation columns
  get_columns_occ = function(self){
    x <- lapply(self$data, function(z) z$get_column_occ())
    x
  },
  # Get weights
  get_weights = function(self){
    x <- lapply(self$data, function(z) z$get_weight())
    x
  },
  # Get ids of all assets in collection
  get_ids = function(self){
    x <- lapply(self$data, function(z) z$id)
    x
  },
  # Search for a specific biodiversity dataset with type
  get_id_byType = function(self, type){
    assertthat::assert_that(is.character(type), !missing(type))
    # Check whether type is correctly set
    if(type %notin% c("Point - Presence only","Point - Presence absence",'Polygon - Presence only','Polygon - Presence absence')){
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
  # Get id by name
  get_id_byName = function(self, name){
    assertthat::assert_that(is.character(name), !missing(name))
    # Get id(s) of dataset with given name
    r <- lapply(self$data, function(z) z$name)
    if(name %in% r){
      r[which(r==name)]
    } else character(0)
  },
  # Show equations of all datasets
  show_equations = function(self, msg = TRUE) {
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
  # Plot the whole collection
  plot = function(self){
    # FIXME: Can quite likely be beautified
    # Get observed columns
    cols <- self$get_columns_occ()

    par.ori <- par(no.readonly = TRUE)
    # Base plot
    g <- ggplot2::ggplot() + ggplot2::geom_sf() + ggplot2::labs( title = self$name())
    # Adding the other elements
    for(dataset in names(cols)){

      if('Polygon - Presence only' == self$get_types()[dataset] ) g <- g + ggplot2::geom_sf(data = st_as_sf(self$get_data(dataset))[cols[[dataset]]], fill = 'lightblue', alpha = .35 )

      if('Polygon - Presence absence' == self$get_types()[dataset] ){
        dd <- st_as_sf(self$get_data(dataset))[cols[[dataset]]]
        dd[[cols[[dataset]]]] <- factor(dd[[cols[[dataset]]]])
        g <- g + ggplot2::geom_sf(data = dd, fill = 'lightgreen', alpha = .35 )
      }

      if('Point - Presence only' == self$get_types()[dataset] ) g <- g + ggplot2::geom_sf(data = st_as_sf(self$get_data(dataset))[cols[[dataset]]], colour = 'grey20', alpha = .5 )

      if('Point - Presence absence' == self$get_types()[dataset] ){
        dd <- st_as_sf(self$get_data(dataset))[cols[[dataset]]]
        dd[[cols[[dataset]]]] <- factor(dd[[cols[[dataset]]]])
        dd$observed <- dd[[cols[[dataset]]]]
        g <- g + ggplot2::geom_sf(data = dd, ggplot2::aes(colour = observed) )
      }
    }
    g
  }

)


#' BiodiversityDataset prototype description
#'
#' @name BiodiversityDataset-class
#' @aliases BiodiversityDataset
#' @keywords bdproto
NULL

#' @export
BiodiversityDataset <- bdproto(
  "BiodiversityDataset",
  name             = character(0),
  id               = character(0),
  equation         = new_waiver(),
  family           = character(0),
  link             = new_waiver(),
  type             = new_waiver(),
  weight           = new_waiver(),
  field_occurrence = character(0),
  data             = new_waiver(),
  # Set new equation
  set_equation = function(self, x){
    assertthat::assert_that(inherits(x, "formula"))
    self$formula <- x
  },
  # Get equation
  get_equation = function(self){
    if(is.Waiver(self$equation)) return('<Default>')
    self$equation
  },
  # Function to print the equation
  show_equation = function(self){
    if(!is.Waiver(equation) && !is.null(equation))
      message(equation)
    else message('None set. Default equation used (response ~ .)')
  },
  # Printing function
  print = function(self){
    message(paste0('Biodiversity data:',
            '\n Name:  ',self$name,
            '\n Type:  ',self$get_type()
    ))
  },
  # Return name
  name = function(self){
    self$name
  },
  # Get Id
  id = function(self){
    self$id
  },
  # Get type
  get_type = function(self, short = FALSE){
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
  # Get field with occurrence information
  get_column_occ = function(self){
    self$field_occurrence
  },
  # Get family
  get_family = function(self){
    assertthat::assert_that(is.character(self$family))
    self$family
  },
  # Get custom link function
  get_link = function(self){
    self$link
  },
  # Get data
  get_data = function(self){
    self$data
  },
  # Get weight
  get_weight = function(self){
    self$weight
  },
  # Print input messages
  show = function(self) {
    self$print()
  },
  # Collect info statistics
  get_observations = function(self) {
    nrow(self$data)
  }
)
