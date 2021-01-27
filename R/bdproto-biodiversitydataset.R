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
#' @aliases BiodiversityDatasetCollection
NULL

#' @export
BiodiversityDatasetCollection <- bdproto(
  "BiodiversityDatasetCollection",
  data = list(), # FIXME: Possibly replace this in the future by a UUID / object-variable sytem
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
  get_types = function(self){
    if (base::length(self$data) > 0)
    {
      return(sapply(self$data, function(z) z$get_type() ))
    } else return(new_waiver())
  },
  # Add a new Biodiversity dataset to this collection
  set_data = function(self, x, value){
    assertthat::assert_that(assertthat::is.string(x),
                            inherits(value, "BiodiversityDataset"))
    self$data[[x]] <- value
    invisible()
  },
  # Get a specific Biodiversity dataset by name
  get_data_object = function(self, x) {
    assertthat::assert_that(assertthat::is.string(x))
    if (!x %in% names(self$data))
      return(new_waiver())
    return(self$data[[x]])
  },
  # Get biodiversity observations
  get_data = function(self, x){
    assertthat::assert_that(assertthat::is.string(x))
    o <- self$get_data_object(x)
    o$get_data()
  },
  # Get coordinates for a given biodiversity dataset. Else return a wkt object
  get_coordinates = function(self, x){
    assertthat::assert_that(assertthat::is.string(x))
    if(x %in% c('poipo','poipa') ){
      o <- self$get_data(x)
      o[,c('X','Y')]
    } else {
      # TODO: WKT to be implemented
      return(new_waiver())
    }
  },
  # Remove a specific biodiversity dataset by name
  rm_data = function(self, x) {
    assertthat::assert_that(assertthat::is.string(x),
                            x %in% names(self$data) )
    self$data[[x]] <- NULL
    invisible()
  },
  # Number of Biodiversity Datasets in connection
  length = function(self) {
    base::length(self$data)
  },
  # Get observations of all datasets
  get_observations = function(self) {
    x <- sapply(self$data, function(z) z$get_observations())
    x
  },
  # Get equations
  get_equations = function(self){
    x <- lapply(self$data, function(z) z$get_equation())
    x
  },
  # Show equations of all datasets
  show_equations = function(self, msg = TRUE) {
    x <- self$get_equations()
    #if(length(x)== 0) return( new_waiver() )
    # new names
    n <- c(poipo = 'Point - Presence only',poipa = 'Point - Presence absence',
      polpo = 'Polygon - Presence only',polpa = 'Polygon - Presence absence')
    names(x) <- as.vector( n[match(names(x), names(n))] )
    # Prettify
    o <- paste0(names(x),":\n ",x,collapse = '\n')
    if(msg) message(o) else o
  }
)


#' BiodiversityDataset prototype description
#'
#' @name BiodiversityDataset-class
#' @aliases BiodiversityDataset
NULL

#' @export
BiodiversityDataset <- bdproto(
  "BiodiversityDataset",
  name         = character(0),
  id           = character(0),
  equation     = new_waiver(),
  type         = new_waiver(),
  data         = new_waiver(),
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
  get_type = function(self){
    switch (self$type,
      poipo = 'Point - Presence only',
      poipa = 'Point - Presence absence',
      polpo = 'Polygon - Presence only',
      polpa = 'Polygon - Presence absence'
    )
  },
  # Get data
  get_data = function(self){
    self$data
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
