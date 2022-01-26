#' @include utils.R bdproto.R bdproto-biodiversitydataset.R bdproto-biodiversitydistribution.R
NULL

#' Add biodiversity point dataset to a distribution object (presence-only)
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param poipo A [`data.frame`], [`sf`] or [`Spatial`]) object of presence-only point occurrences.
#' @param name The name of the biodiversity dataset used as internal identifier
#' @param field_occurrence A [`numeric`] or [`character`] location of biodiversity point records.
#' @param formula A [`character`] or [`formula`] object to be passed. Default is to use all covariates (if specified)
#' @param family A [`character`] stating the family to be used (Default: Poisson)
#' @param separate_intercept A [`boolean`] value stating whether a separate intercept is to be added in
#' shared likelihood models for engines [engine_inla], [engine_inlabru] and [engine_stan].
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
  function(x, poipo, name = NULL, field_occurrence = "Observed", formula = NULL, family = "poisson", separate_intercept = TRUE, ...) standardGeneric("add_biodiversity_poipo"))

# TODO: Support supplement of other object types, such as data.frame, sp, etc...

#' @name add_biodiversity_poipo
#' @rdname add_biodiversity_poipo
#' @usage \S4method{add_biodiversity_poipo}{BiodiversityDistribution,sf}(x, poipo)
methods::setMethod(
  "add_biodiversity_poipo",
  methods::signature(x = "BiodiversityDistribution", poipo = "sf"),
  function(x, poipo, name = NULL, field_occurrence = "Observed", formula = NULL, family = "poisson", separate_intercept = TRUE, ...) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            inherits(poipo, "Spatial") || inherits(poipo, "sf") || inherits(poipo, "data.frame") || inherits(poipo, "tibble"),
                            assertthat::is.scalar(field_occurrence), assertthat::has_name(poipo, field_occurrence),
                            inherits(formula, "formula") || is.null(formula) || is.character(formula),
                            is.character(family),
                            is.logical(separate_intercept)
                            )
    assertthat::assert_that(length(unique(poipo[[field_occurrence]])) <= 2,
                            msg = "More 2 unique values. Specify a column.")

    # Messager
    if(getOption("ibis.setupmessages")) myLog("[Setup]","green","Adding poipo dataset...")

    # Get only those records that fall onto the background data
    suppressMessages( poipo <- point_in_polygon(poly = x$background, points = poipo) )

    # Convert formula if necessary
    formula <- to_formula(formula)

    # Create a new id for this dataset
    id <- new_id()

    # Finally set the data to the BiodiversityDistribution object
    x$set_biodiversity(
      id,
      bdproto(NULL, BiodiversityDataset,
              name = ifelse(is.null(name), "Species: ",name),
              id = id,
              equation = formula,
              family = family,
              type = "poipo",
              field_occurrence = field_occurrence,
              data = format_biodiversity_data(poipo,field_occurrence),
              use_intercept = separate_intercept
      )
    )
  }
)

#' Add biodiversity point dataset to a distribution object (presence-absence)
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param poipa A [`data.frame`], [`sf`] or [`Spatial`]) object of presence-absence point occurrences.
#' @param name The name of the biodiversity dataset used as internal identifier
#' @param field_occurrence A [`numeric`] or [`character`] location of biodiversity point records indicating presence/absence
#' @param formula A [`character`] or [`formula`] object to be passed. Default is to use all covariates (if specified)
#' @param family A [`character`] stating the family to be used (Default: binomial)
#' @param separate_intercept A [`boolean`] value stating whether a separate intercept is to be added in
#' shared likelihood models for engines [engine_inla], [engine_inlabru] and [engine_stan].
#' @param ... Other parameters passed down
#'
#' @details Say something about presence-absence biodiversity records in \pkg{ibis}
#' @section Notes:
#' @references
#'
#' @examples
#' \dontrun{
#'  TBD
#' }
#' @name add_biodiversity_poipa
NULL

#' @name add_biodiversity_poipa
#' @rdname add_biodiversity_poipa
#' @exportMethod add_biodiversity_poipa
#' @export
methods::setGeneric(
  "add_biodiversity_poipa",
  signature = methods::signature("x", "poipa"),
  function(x, poipa, name = NULL, field_occurrence = "Observed", formula = NULL, family = "binomial", separate_intercept = TRUE, ...) standardGeneric("add_biodiversity_poipa"))

#' @name add_biodiversity_poipa
#' @rdname add_biodiversity_poipa
#' @usage \S4method{add_biodiversity_poipa}{BiodiversityDistribution,sf}(x, poipa)
methods::setMethod(
  "add_biodiversity_poipa",
  methods::signature(x = "BiodiversityDistribution", poipa = "sf"),
  function(x, poipa, name = NULL, field_occurrence = "Observed", formula = NULL, family = "binomial", separate_intercept = TRUE,  ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            inherits(poipa, "Spatial") || inherits(poipa, "sf") || inherits(poipa, "data.frame") || inherits(poipa, "tibble"),
                            assertthat::is.scalar(field_occurrence), assertthat::has_name(poipa, field_occurrence),
                            inherits(formula, "formula") || is.null(formula) || is.character(formula),
                            is.character(family),
                            is.logical(separate_intercept)
    )
    assertthat::assert_that(length(unique(poipa[[field_occurrence]])) == 2,
                            msg = "Presence-Absence requires at exactly 2 unique values.")

    # Messager
    if(getOption("ibis.setupmessages")) myLog("[Setup]","green","Adding poipa dataset...")

    # Get only those records that fall onto the background data
    suppressMessages( poipa <- point_in_polygon(poly = x$background, points = poipa) )

    # Record presence absence to 0 and 1
    if(is.character(poipa[[field_occurrence]]) ){
      # TODO:
      stop("Guessing conversion to be coded")
    }
    poipa[[field_occurrence]] <- as.numeric(poipa[[field_occurrence]]) # Convert to numeric for ocassional factor formatting

    # Convert formula if necessary
    formula = to_formula(formula)

    # Create a new id for this dataset
    id = new_id()

    # Finally set the data to the BiodiversityDistribution object
    x$set_biodiversity(
      id,
      bdproto(NULL, BiodiversityDataset,
              name = ifelse(is.null(name), "Species: ",name),
              id = id,
              equation = formula,
              family = family,
              type = "poipa",
              field_occurrence = field_occurrence,
              data = format_biodiversity_data(poipa,field_occurrence),
              use_intercept = separate_intercept
      )
    )
  }
)


#' Add biodiversity polygon dataset to a distribution object (presence-only)
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param polpo A [`sf`] or [`Spatial`]) polygon object of presence-only occurrences.
#' @param name The name of the biodiversity dataset used as internal identifier
#' @param field_occurrence A [`numeric`] or [`character`] location of biodiversity point records.
#' @param formula A [`character`] or [`formula`] object to be passed. Default is to use all covariates (if specified)
#' @param family A [`character`] stating the family to be used (Default: Poisson)
#' @param simulate Simulate poipa points within its boundaries. Result are passed to [`add_biodiversity_poipa`] (Default: FALSE)
#' @param simulate_points A [`numeric`] number of points to be created by simulation
#' @param simulate_weights A [`Raster`] layer describing an eventual preference for simulation (Default: NULL)
#' @param simulate_strategy A [`character`] stating the strategy for sampling. Can be set to either
#' 'random' or 'regular', the latter requiring a raster supplied in the [simulate_weights]
#' parameter.
#' @param separate_intercept A [`boolean`] value stating whether a separate intercept is to be added in
#' shared likelihood models for engines [engine_inla], [engine_inlabru] and [engine_stan].
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
  function(x, polpo, name = NULL, field_occurrence = "Observed", formula = NULL, family = "poisson",
           simulate = FALSE, simulate_points = 100, simulate_weights = NULL, simulate_strategy = "random", separate_intercept = TRUE, ...) standardGeneric("add_biodiversity_polpo"))

#' @name add_biodiversity_polpo
#' @rdname add_biodiversity_polpo
#' @usage \S4method{add_biodiversity_polpo}{BiodiversityDistribution,sf}(x, polpo)
methods::setMethod(
  "add_biodiversity_polpo",
  methods::signature(x = "BiodiversityDistribution", polpo = "sf"),
  function(x, polpo, name = NULL, field_occurrence = "Observed", formula = NULL, family = "poisson",
           simulate = FALSE, simulate_points = 100, simulate_weights = NULL, simulate_strategy = "random", separate_intercept = TRUE, ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            inherits(polpo, "Spatial") || inherits(polpo, "sf") || inherits(polpo, "data.frame") || inherits(polpo, "tibble"),
                            assertthat::is.scalar(field_occurrence), assertthat::has_name(polpo, field_occurrence),
                            inherits(formula, "formula") || is.null(formula) || is.character(formula),
                            is.character(family),
                            assertthat::is.flag(simulate), is.numeric(simulate_points),
                            is.null(simulate_weights) || inherits(simulate_weights, "Raster"),
                            is.logical(separate_intercept)
    )

    # Check type and ensure that is a polygon
    assertthat::assert_that(all( unique( st_geometry_type(polpo) ) %in% c("POLYGON","MULTIPOLYGON") ),
                            msg = "This method works for spatial data of type polygon only.")

    assertthat::assert_that(length(unique(polpo[[field_occurrence]])) <= 2,
                            msg = "More 2 unique values. Specify a column ")

    # Messager
    if(getOption("ibis.setupmessages")) myLog("[Setup]","green","Adding polpo dataset...")

    # Simulate presence absence points rather than using the range directly
    if(simulate){
      if(family == "poisson") warning("Simulated points created. Binomial distribution is advised.")
      # Simulation strategy
      simulate_strategy <- match.arg(simulate_strategy, c('random', 'regular'), several.ok = FALSE)

      if(!is.null(simulate_weights)){
        # Crop to target range
        simulate_weights <- raster::crop(simulate_weights, polpo)
        # Normalize the weight layer if is not a factorized, else set everything to 1
        if(is.null(levels(simulate_weights))) simulate_weights <- predictor_transform(simulate_weights, "norm") else simulate_weights[simulate_weights>0] <- 1

        # Weighted sampling on background raster, the greater the value, the more likely sampled points
        ptscell <- sample(which(!is.na(simulate_weights[])),
                          size = simulate_points,
                          prob = simulate_weights[which(!is.na(simulate_weights[]))],
                          replace = TRUE)
        poipa_on <- as.data.frame(raster::xyFromCell(simulate_weights, ptscell))
        poipa_on <- sf::st_as_sf(poipa_on, coords = c("x","y"),crs = sf::st_crs(simulate_weights))

      } else {
      # Simply sample presence points within as determined
        suppressMessages(
          poipa_on <- sf::st_as_sf(
            sf::st_sample(x = polpo, size = simulate_points, type = simulate_strategy)
          )
        )
      }

      names(poipa_on) <- "geometry"; st_geometry(poipa_on) <- "geometry"
      poipa_on[[field_occurrence]] <- 1
      poipa_on$x <- st_coordinates(poipa_on)[,1];poipa_on$y <- st_coordinates(poipa_on)[,2]
      # Get absence data
      # FIXME: Quick fix. Ideally cookie cut the range out instead.
      suppressMessages(
        poipa_off <- sf::st_as_sf(
          sf::st_sample(x = x$background, size = simulate_points*2, type = "random")
        )
      )
      names(poipa_off) <- "geometry"; st_geometry(poipa_off) <- "geometry"
      # Remove points on the range
      suppressMessages(
        wi <- sf::st_within(poipa_off, polpo, sparse = FALSE)
      )
      poipa_off <- poipa_off[!apply(wi, 1, any),]
      poipa_off[[field_occurrence]] <- 0
      poipa_off$x <- st_coordinates(poipa_off)[,1];poipa_off$y <- st_coordinates(poipa_off)[,2]

      poipa <- rbind(poipa_on,poipa_off)

      # Add simulated poipa object instead
      add_biodiversity_poipa(x, poipa = poipa, name = paste0(name, "_simulated"),
                             field_occurrence = field_occurrence, formula = formula, family = family, ... )
    } else {
      # Convert formula if necessary
      formula = to_formula(formula)

      # Create a new id for this dataset
      id = new_id()

      # Finally set the data to the BiodiversityDistribution object
      x$set_biodiversity(
        id,
        bdproto(NULL, BiodiversityDataset,
                name = ifelse(is.null(name), "Species: ", name),
                id = id,
                equation = formula,
                family = family,
                type = "polpo",
                field_occurrence = field_occurrence,
                data = format_biodiversity_data(polpo,field_occurrence),
                use_intercept = separate_intercept
        )
      )
    }
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
  if(!(inherits(x,'sf') )){
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

