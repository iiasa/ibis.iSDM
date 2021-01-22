#' @include utils.R bdproto.R bdproto-biodiversitydataset.R bdproto-biodiversitydistribution.R
NULL

#' Add biodiversity point dataset to a distribution object (presence-only)
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param po A [`data.frame`], [`sf`] or [`Spatial`]) object of presence-only point occurrences.
#' @param field_occurrence A [`numeric`] or [`character`] location of biodiversity point records.
#' @param ... Other parameters passed down
#'
#' @details Say something about presence only biodiversity records in \pkg{ibis}
#' @section Notes:
#' @references
#'
#' @examples
#' \dontrun{
#'  TBD
#' }
#' @name add_biodiversity_poipo
NULL

#' @name add_biodiversity_poipo
#' @rdname add_biodiversity_poipo
#' @exportMethod add_biodiversity_poipo
#' @export
methods::setGeneric(
  "add_biodiversity_poipo",
  signature = methods::signature("x", "poipo"),
  function(x, poipo, name = NULL, field_occurrence = 'Observed',...) standardGeneric("add_biodiversity_poipo"))

# TODO: Support supplement of other object types, such as data.frame, sp, etc...

#' @name add_biodiversity_poipo
#' @rdname add_biodiversity_poipo
#' @usage \S4method{add_biodiversity_poipo}{BiodiversityDistribution,sf}(x, poipo)
methods::setMethod(
  "add_biodiversity_poipo",
  methods::signature(x = "BiodiversityDistribution", poipo = "sf"),
  function(x, poipo, name = NULL, field_occurrence = 'Observed', ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            inherits(poipo,'Spatial') || inherits(poipo,'sf') || inherits(poipo,'data.frame') || inherits(poipo,'tibble'),
                            assertthat::is.scalar(field_occurrence), assertthat::has_name(poipo,field_occurrence)
                            )
    assertthat::assert_that(length(unique(poipo[[field_occurrence]])) <= 2,
                            msg = "More 2 unique values. Specify a column ")

    # Assess whether poipo data already exists in the distribution object
    if(!is.Waiver( x$biodiversity$get_data('poipo') )) message('Overwriting existing poipo data.')

    # Finally set the data to the BiodiversityDistribution object
    x$biodiversity$set_data(
      'poipo',
      bdproto(NULL, BiodiversityDataset,
              name = ifelse(is.null(name), 'Species: ',name),
              id = new_id(),
              type = 'poipo',
              data = format_biodiversity_data(poipo,field_occurrence)
      )
    )
    return(x)
  }
)

#' Add biodiversity polygon dataset to a distribution object (presence-only)
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param po A [`sf`] or [`Spatial`]) object of presence-only occurrences.
#' @param field_occurrence A [`numeric`] or [`character`] location of biodiversity records.
#' @param ... Other parameters passed down
#'
#' @details Say something about presence only biodiversity records in \pkg{ibis} and
#' integration of ranges
#' @section Notes:
#' @references
#'
#' @examples
#' \dontrun{
#'  TBD
#' }
#' @name add_biodiversity_polpo
NULL

#' @name add_biodiversity_polpo
#' @rdname add_biodiversity_polpo
#' @exportMethod add_biodiversity_polpo
#' @export
methods::setGeneric(
  "add_biodiversity_polpo",
  signature = methods::signature("x", "polpo"),
  function(x, polpo, name = NULL, field_occurrence = 'Observed',...) standardGeneric("add_biodiversity_polpo"))

# TODO: Support supplement of other object types, such as data.frame, sp, etc...

#' @name add_biodiversity_polpo
#' @rdname add_biodiversity_polpo
#' @usage \S4method{add_biodiversity_polpo}{BiodiversityDistribution,sf}(x, polpo)
methods::setMethod(
  "add_biodiversity_polpo",
  methods::signature(x = "BiodiversityDistribution", polpo = "sf"),
  function(x, polpo, name = NULL, field_occurrence = 'Observed', ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            inherits(polpo,'Spatial') || inherits(polpo,'sf') || inherits(polpo,'data.frame') || inherits(polpo,'tibble'),
                            assertthat::is.scalar(field_occurrence), assertthat::has_name(polpo,field_occurrence)
    )
    assertthat::assert_that(length(unique(polpo[[field_occurrence]])) <= 2,
                            msg = "More 2 unique values. Specify a column ")

    # Assess whether poipo data already exists in the distribution object
    if(!is.Waiver( x$biodiversity$get_data('polpo') )) message('Overwriting existing polpo data.')

    # Finally set the data to the BiodiversityDistribution object
    x$biodiversity$set_data(
      'polpo',
      bdproto(NULL, BiodiversityDataset,
              name = ifelse(is.null(name), 'Species: ',name),
              id = new_id(),
              type = 'polpo',
              data = format_biodiversity_data(polpo,field_occurrence)
      )
    )
    return(x)
  }
)

#' Format biodiversity dataset to standardized format
#'
#' @param x A [`data.frame`], [`sf`] or [`Spatial`]) object of biodiversity information
#' @param field_occurrence A [`numeric`] or [`character`] location of biodiversity records.
#' @param field_space A [`vector`] on the column names (Default: 'X','Y')
#' @param ... Other parameters passed down
#'
#' @import sf
#' @name format_biodiversity_data
#' @rdname format_biodiversity_data

format_biodiversity_data <- function(x, field_occurrence, field_space = c('X','Y'),...){
  # Final checks
  if(inherits(x,'sf')) assertthat::assert_that(unique(sf::st_geometry_type(x)) %in% c('POINT','MULTIPOINT',
                                                                                      'POLYGON','MULTIPOLYGON'))

  # If data.frame or tibble, check whether coordinates are in there
  if((!inherits(x,'sf') && is.data.frame(x)) || tibble::is_tibble(x)){
    # Check whether field_space columns are present

    # FIXME: Implement a lookup of commonly used column names
    # if(!all(assertthat::has_name(x, field_space))){
    #   # Spatial column suggestions
    #   spat_cols <- c('x','y','long','lat','long','longitude','latitude')
    #   assertthat::assert_that(
    #     any(assertthat::has_name(x,spat_cols)) || any( has_name(x, sapply(spat_cols, capitalize_text)) ),
    #     msg = 'No spatial column found in the dataset. Specify manually or set to [x] and [y].'
    #   )
    # }
    # TODO: Make sure this function is comprehensive and captures all possible spatial formats
    stop('To be coded')
    assertthat::assert_that( all(assertthat::has_name(x, field_space)),
                             msg ='No spatial column found in the dataset. Specify manually or set to [X] and [Y].')

    # Select and format
    out <- subset(x, select = c(field_space, field_occurrence) ) %>%
      as_tibble()
  } else {
    if(inherits(x, 'Spatial')) x <- sf::st_as_sf(x) # First convert to sf
    #if(inherits(x,'sf')) coords <- sf::st_coordinates(x) %>% tibble::as_tibble()

    if(unique(sf::st_geometry_type(x)) %in% c('POINT','MULTIPOINT')){
      # Take target column and append coordinates to it
      out <- cbind(subset(x, select = field_occurrence),
                   sf::st_coordinates(x)) %>% tibble::as_tibble()
    } else if(unique(sf::st_geometry_type(x)) %in% c('POLYGON','MULTIPOLYGON')){
      # Return target column and spatial object as such
      out <- subset(x, select = field_occurrence)
    }
  }
  out
}

