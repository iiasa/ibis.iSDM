#' @include class-biodiversitydataset.R class-biodiversitydistribution.R
NULL

#' Add biodiversity point dataset to a distribution object (presence-only)
#'
#' @description This function adds a presence-only biodiversity dataset to a
#'   distribution object.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param poipo A [`data.frame`] or [`sf`] object of presence-only point occurrences.
#' @param name The name of the biodiversity dataset used as internal identifier.
#' @param field_occurrence A [`numeric`] or [`character`] location of biodiversity point records.
#' @param formula A [`character`] or [`formula`] object to be passed. Default is
#' to use all covariates (if specified).
#' @param family A [`character`] stating the family to be used (Default: \code{'Poisson'}).
#' @param link A [`character`] to overwrite the default link function (Default: \code{NULL}).
#' @param weight A [`numeric`] value acting as a multiplier with regards to any
#' weights used in the modelling. Larger weights indicate higher weighting
#' relative to any other datasets. By default set to \code{1} if only one
#' dataset is added. A [`vector`] is also supported but must be of the same
#' length as \code{"poipo"}.
#' **Note: Weights are reformated to the inverse for models with area offsets (e.g. 5 is converted to 1/5).**
#' @param separate_intercept A [`logical`] value stating whether a separate
#' intercept is to be added in shared likelihood models for engines
#' [engine_inla], [engine_inlabru] and [engine_stan]. Otherwise ignored.
#' @param docheck [`logical`] on whether additional checks should be performed
#' (e.g. intersection tests) (Default: \code{TRUE}).
#' @param pseudoabsence_settings Either \code{NULL} or a
#' [`pseudoabs_settings()`] created settings object.
#' @param ... Other parameters passed down to the object. Normally not used
#' unless described in details.
#'
#' @details This function allows to add presence-only biodiversity records to a
#' [distribution] \pkg{ibis.iSDM} Presence-only data are usually modelled
#' through an inferential model (see Guisan and Zimmerman, 2000) that relate
#' their occurrence in relation to environmental covariates to a selected
#' sample of 'background' points. The most common approach for estimation and
#' the one supported by this type of dataset are poisson-process models (PPM)
#' in which presence-only points are fitted through a down-weighted Poisson
#' regression. See Renner et al. 2015 for an overview.
#'
#' @returns Adds biodiversity data to [distribution] object.
#'
#' @references
#' * Guisan A. and Zimmerman N. 2000. Predictive habitat distribution models in ecology.
#' Ecol. Model. 135: 147–186.
#' * Renner, I. W., J. Elith, A. Baddeley, W. Fithian, T. Hastie, S. J. Phillips,
#' G. Popovic, and D. I. Warton. 2015. Point process models for presence-only analysis.
#' Methods in Ecology and Evolution 6:366–379.
#'
#' @family add_biodiversity
#' @keywords biodiversity
#'
#' @examples
#' # Load background
#' background <- terra::rast(system.file('extdata/europegrid_50km.tif',
#' package='ibis.iSDM',mustWork = TRUE))
#' # Load virtual species
#' virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg',
#' package='ibis.iSDM',mustWork = TRUE),'points',quiet = TRUE)
#' # Define model
#' x <- distribution(background) |>
#' add_biodiversity_poipo(virtual_points, field_occurrence = "Observed")
#'
#' @import sf
#'
#' @name add_biodiversity_poipo
NULL

#' @rdname add_biodiversity_poipo
#' @export
methods::setGeneric(
  "add_biodiversity_poipo",
  signature = methods::signature("x", "poipo"),
  function(x, poipo, name = NULL, field_occurrence = "observed", formula = NULL, family = "poisson", link = NULL,
           weight = 1, separate_intercept = TRUE, docheck = TRUE, pseudoabsence_settings = NULL, ...) {standardGeneric("add_biodiversity_poipo") })

#' @rdname add_biodiversity_poipo
methods::setMethod(
  "add_biodiversity_poipo",
  methods::signature(x = "BiodiversityDistribution", poipo = "sf"),
  function(x, poipo, name = NULL, field_occurrence = "observed", formula = NULL, family = "poisson", link = NULL,
           weight = 1, separate_intercept = TRUE, docheck = TRUE, pseudoabsence_settings = NULL, ...) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            inherits(poipo, "Spatial") || inherits(poipo, "sf") || inherits(poipo, "data.frame") || inherits(poipo, "tibble"),
                            assertthat::is.scalar(field_occurrence), assertthat::has_name(poipo, field_occurrence),
                            inherits(formula, "formula") || is.null(formula) || is.character(formula),
                            is.character(family), is.null(link) || is.character(link),
                            is.logical(separate_intercept),
                            is.logical(docheck),
                            is.numeric(weight) && all(weight > 0)
                            )

    # Messenger
    if(getOption("ibis.setupmessages", default = TRUE)) myLog("[Setup]","green","Adding poipo dataset...")

    # Transform to background for analysis
    if(sf::st_crs(x$background) != sf::st_crs(poipo)){
      poipo <- poipo |> sf::st_transform(crs = sf::st_crs(x$background))
    }

    if(docheck){
      # Get only those records that fall onto the background data
      suppressMessages( poipo <- point_in_polygon(poly = x$background, points = poipo) )
    }

    # Convert occurrence field to numeric
    if(is.factor(poipo[[field_occurrence]])) {
      poipo[[field_occurrence]] <- as.numeric(as.character(poipo[[field_occurrence]]))
    }

    # Check whether there are any absence point, if so stop with error
    if(any(poipo[[field_occurrence]] == 0)){
      if(getOption("ibis.setupmessages", default = TRUE)) myLog("[Setup]","yellow",
                                                "Absence points found. Potentially this data needs to be added as presence-absence instead?")
    }

    # Convert formula if necessary
    formula <- to_formula(formula)

    # Create a new id for this dataset
    id <- new_id()

    # Check that weights are correctly set
    if(length(weight)>1){
      assertthat::assert_that(length(weight) == nrow(poipo))
    } else {
      weight <- rep(weight, nrow(poipo))
    }

    # Create a new biodiversity dataset
    bd <- BiodiversityDataset$new(
      name = ifelse(is.null(name), substring(id[[1]], 1, 8), name),
      id = id,
      equation = formula,
      family = family,
      link = link,
      type = "poipo",
      weight = weight,
      field_occurrence = field_occurrence,
      data = format_biodiversity_data(poipo, field_occurrence),
      use_intercept = separate_intercept,
      pseudoabsence_settings = pseudoabsence_settings
    )

    # Finally set the data to the BiodiversityDistribution object
    y <- x$clone(deep = TRUE)
    y$set_biodiversity(id, bd)
  }
)

#' Add biodiversity point dataset to a distribution object (presence-absence).
#'
#' @description This function adds a presence-absence biodiversity dataset to a
#' distribution object. Opposed to presence-only data, presence-absence
#' biodiversity records usually originate from structured biodiversity surveys
#' where the absence of a species in a given region was specifically assessed.
#'
#' If it is the analysts choice it is also possible to format presence-only
#' biodiversity data into a presence-absence form, by adding pseudo-absence
#' through [`add_pseudoabsence`]. See the help file for more information.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param poipa A [`data.frame`] or [`sf`] object of presence-absence point
#' occurrences.
#' @param name The name of the biodiversity dataset used as internal identifier.
#' @param field_occurrence A [`numeric`] or [`character`] location of
#' biodiversity point records indicating presence/absence. By default set to
#' \code{"observed"} and an error will be thrown if a [`numeric`] column with
#'that name does not exist.
#' @param formula A [`character`] or [`formula`] object to be passed. Default (\code{NULL})
#' is to use all covariates.
#' @param family A [`character`] stating the family to be used  (Default: \code{'binomial'}).
#' @param link A [`character`] to overwrite the default link function (Default: \code{NULL}).
#' @param weight A [`numeric`] value acting as a multiplier with regards to any
#' weights used in the modelling. Larger weights indicate higher weighting
#' relative to any other datasets. By default set to \code{1} if only one
#' dataset is added. A [`vector`] is also supported but must be of the same
#' length as parameter \code{"poipa"}.
#' @param separate_intercept A [`logical`] value stating whether a separate
#' intercept is to be added in. shared likelihood models for engines
#' [engine_inla], [engine_inlabru] and [engine_stan].
#' @param docheck [`logical`] on whether additional checks should be performed
#' (e.g. intersection tests) (Default: \code{TRUE}).
#' @param ... Other parameters passed down.
#'
#' @details By default, the logit link function is used in a logistic regression
#' setting unless the specific engine does not support generalised linear
#' regressions (e.g. [engine_bart]).
#'
#' @returns Adds biodiversity data to [distribution] object.
#'
#' @references
#' * Renner, I. W., J. Elith, A. Baddeley, W. Fithian, T. Hastie, S. J. Phillips, G.
#' Popovic, and D. I. Warton. 2015. Point process models for presence-only analysis. Methods in Ecology and Evolution 6:366–379.
#' * Guisan A. and Zimmerman N. 2000. Predictive habitat distribution models in
#' ecology. Ecol. Model. 135: 147–186.
#'
#' @family add_biodiversity
#' @keywords biodiversity
#'
#' @examples
#' \dontrun{
#' # Define model
#' x <- distribution(background) |> add_biodiversity_poipa(virtual_species)
#' }
#'
#' @name add_biodiversity_poipa
NULL

#' @rdname add_biodiversity_poipa
#' @export
methods::setGeneric(
  "add_biodiversity_poipa",
  signature = methods::signature("x", "poipa"),
  function(x, poipa, name = NULL, field_occurrence = "observed", formula = NULL, family = "binomial", link = NULL,
           weight = 1, separate_intercept = TRUE, docheck = TRUE, ...) standardGeneric("add_biodiversity_poipa"))

#' @rdname add_biodiversity_poipa
methods::setMethod(
  "add_biodiversity_poipa",
  methods::signature(x = "BiodiversityDistribution", poipa = "sf"),
  function(x, poipa, name = NULL, field_occurrence = "observed", formula = NULL, family = "binomial", link = NULL,
           weight = 1, separate_intercept = TRUE, docheck = TRUE, ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            inherits(poipa, "Spatial") || inherits(poipa, "sf") || inherits(poipa, "data.frame") || inherits(poipa, "tibble"),
                            assertthat::is.scalar(field_occurrence), assertthat::has_name(poipa, field_occurrence),
                            inherits(formula, "formula") || is.null(formula) || is.character(formula),
                            is.character(family),
                            is.null(link) || is.character(link),
                            is.logical(separate_intercept),
                            is.numeric(weight) && all(weight > 0),
                            is.logical(docheck)
    )
    # Raise a warning if there NA data
    assertthat::assert_that(!anyNA(poipa[[field_occurrence]]) || !is.null(poipa[[field_occurrence]]),
                            msg = "NA observations or NULL occurrence field?")

    # Check for that presence and absence are there.
    assertthat::assert_that(length(unique(poipa[[field_occurrence]])) == 2,
                            msg = "Presence-Absence requires at exactly 2 unique values.")

    # Messenger
    if(getOption("ibis.setupmessages", default = TRUE)) myLog("[Setup]","green","Adding poipa dataset...")

    # Transform to background for analysis
    if(sf::st_crs(x$background) != sf::st_crs(poipa)){
      poipa <- poipa |> sf::st_transform(crs = sf::st_crs(x$background))
    }

    if(docheck){
      # Get only those records that fall onto the background data
      suppressMessages( poipa <- point_in_polygon(poly = x$background, points = poipa) )
    }

    # Record presence absence to 0 and 1
    if(is.factor(poipa[[field_occurrence]])) {
      poipa[[field_occurrence]] <- as.numeric(as.character(poipa[[field_occurrence]]))
    }

    # Convert formula if necessary
    formula = to_formula(formula)

    # Create a new id for this dataset
    id = new_id()

    # Check that weights are correctly set
    if(length(weight)>1) {
      assertthat::assert_that(length(weight) == nrow(poipa))
    } else {
      weight <- rep(weight, nrow(poipa))
    }

    # Define the biodiversity object
    bd <- BiodiversityDataset$new(
      name = ifelse(is.null(name), substring(id[[1]], 1, 8), name),
      id = id,
      equation = formula,
      family = family,
      link = link,
      type = "poipa",
      weight = weight,
      field_occurrence = field_occurrence,
      data = format_biodiversity_data(poipa, field_occurrence),
      use_intercept = separate_intercept
    )

    # Finally set the data to the BiodiversityDistribution object
    y <- x$clone(deep = TRUE)
    y$set_biodiversity(id, bd)
  }
)

#' Add biodiversity polygon dataset to a distribution object (presence-only)
#'
#' @description This function can be used to add a [`sf`] polygon dataset to an
#' existing distribution object. Presence-only polygon data is treated
#' differential than point data in some engines particular through the way that
#' points are generated.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param polpo A [`sf`] polygon object of presence-only occurrences.
#' @param name The name of the biodiversity dataset used as internal identifier.
#' @param field_occurrence A [`numeric`] or [`character`] location of
#' biodiversity point records.
#' @param formula A [`character`] or [`formula`] object to be passed. Default is
#' to use all covariates (if specified).
#' @param family A [`character`] stating the family to be used (Default: \code{poisson}).
#' @param link A [`character`] to overwrite the default link function (Default: \code{NULL}).
#' @param weight A [`numeric`] value acting as a multiplier with regards to any
#' weights used in the modelling. Larger weights indicate higher weighting
#' relative to any other datasets. By default set to \code{1} if only one
#' dataset is added. A [`vector`] is also supported but must be of the same
#' length as \code{"polpo"}.
#' @param simulate Simulate poipo points within its boundaries. Result are
#' passed to [`add_biodiversity_poipo`] (Default: \code{FALSE}).
#' @param simulate_points A [`numeric`] number of points to be created by
#' simulation (Default: \code{100}).
#' @param simulate_bias A [`SpatRaster`] layer describing an eventual preference
#' for simulation (Default: \code{NULL}).
#' @param simulate_strategy A [`character`] stating the strategy for sampling.
#' Can be set to either. \code{'random'} or \code{'regular'}, the latter
#' requiring a raster supplied in the \code{'simulate_weights'} parameter.
#' @param separate_intercept A [`logical`] value stating whether a separate
#' intercept is to be added in shared likelihood models for engines
#' [engine_inla], [engine_inlabru] and [engine_stan].
#' @param docheck [`logical`] on whether additional checks should be performed
#' (e.g. intersection tests) (Default: \code{TRUE}).
#' @param pseudoabsence_settings Either \code{NULL} or a [`pseudoabs_settings()`] created settings object.
#' @param ... Other parameters passed down.
#'
#' @details The default approach for polygon data is to sample presence-only
#' points across the region of the polygons. This function thus adds as a
#' wrapper to [`add_biodiversity_poipo()`] as presence-only points are created
#' by the model. If no points are simulated directly (Default) then the
#' polygon is processed by [`train()`] by creating regular point data over the
#' supplied predictors.
#'
#' Use [`add_biodiversity_polpa()`] to create binomial distributed
#' inside-outside points for the given polygon!
#'
#' For an integration of range data as predictor or offset, see
#'
#' [`add_predictor_range()`] and [`add_offset_range()`] instead.
#'
#' @returns Adds biodiversity data to [distribution] object.
#'
#' @family add_biodiversity
#' @keywords biodiversity
#'
#' @examples
#' \dontrun{
#'  x <- distribution(mod) |>
#'    add_biodiversity_polpo(protectedArea)
#' }
#' @name add_biodiversity_polpo
NULL

#' @rdname add_biodiversity_polpo
#' @export
methods::setGeneric(
  "add_biodiversity_polpo",
  signature = methods::signature("x", "polpo"),
  function(x, polpo, name = NULL, field_occurrence = "observed", formula = NULL, family = "poisson", link = NULL,
           weight = 1, simulate = FALSE, simulate_points = 100, simulate_bias = NULL, simulate_strategy = "random",
           separate_intercept = TRUE, docheck = TRUE, pseudoabsence_settings = NULL, ...) standardGeneric("add_biodiversity_polpo"))

#' @rdname add_biodiversity_polpo
methods::setMethod(
  "add_biodiversity_polpo",
  methods::signature(x = "BiodiversityDistribution", polpo = "sf"),
  function(x, polpo, name = NULL, field_occurrence = "observed", formula = NULL,family = "poisson", link = NULL,
           weight = 1, simulate = FALSE, simulate_points = 100, simulate_bias = NULL, simulate_strategy = "random",
           separate_intercept = TRUE, docheck = TRUE, pseudoabsence_settings = NULL, ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            inherits(polpo, "Spatial") || inherits(polpo, "sf") || inherits(polpo, "data.frame") || inherits(polpo, "tibble"),
                            assertthat::is.scalar(field_occurrence), assertthat::has_name(polpo, field_occurrence),
                            inherits(formula, "formula") || is.null(formula) || is.character(formula),
                            is.character(family),
                            is.null(link) || is.character(link),
                            assertthat::is.flag(simulate), is.numeric(simulate_points),
                            is.null(simulate_bias) || inherits(simulate_bias, "SpatRaster"),
                            is.logical(separate_intercept),
                            is.numeric(weight) && all(weight > 0)
    )

    # Check type and ensure that is a polygon
    assertthat::assert_that(all( unique( sf::st_geometry_type(polpo) ) %in% c("POLYGON","MULTIPOLYGON") ),
                            msg = "This method works for spatial data of type polygon only.")

    assertthat::assert_that(length(unique(polpo[[field_occurrence]])) <= 2,
                            msg = "More 2 unique values. Specify a column ")

    # Messenger
    if(getOption("ibis.setupmessages", default = TRUE)) myLog("[Setup]","green","Adding polpo dataset...")

    # Transform to background for analysis
    if(sf::st_crs(x$background) != sf::st_crs(polpo)){
      polpo <- polpo |> sf::st_transform(crs = sf::st_crs(x$background))
    }

    # Record presence absence to 0 and 1
    if(is.factor(polpo[[field_occurrence]])) {
      polpo[[field_occurrence]] <- as.numeric(as.character(polpo[[field_occurrence]]))
    }

    # Simulate presence absence points rather than using the range directly
    if(simulate){
      # Simulation strategy
      simulate_strategy <- match.arg(simulate_strategy, c('random', 'regular'), several.ok = FALSE)

      if(!is.null(simulate_bias)){
        # Crop to target range
        simulate_bias <- terra::crop(simulate_bias, polpo)
        # Normalize the weight layer if is not a factorized, else set everything to 1
        if(is.null(levels(simulate_bias))) simulate_bias <- predictor_transform(simulate_bias, "norm") else simulate_bias[simulate_bias>0] <- 1

        # Weighted sampling on background raster, the greater the value, the
        # more likely sampled points
        ptscell <- sample(which(!is.na(simulate_bias[])),
                          size = simulate_points,
                          prob = simulate_bias[which(!is.na(simulate_bias[]))],
                          replace = TRUE)
        poipo_on <- terra::as.data.frame(terra::xyFromCell(simulate_bias, ptscell))
        poipo_on <- sf::st_as_sf(poipo_on, coords = c("x","y"), crs = sf::st_crs(simulate_bias))

      } else {
      # Simply sample presence points within as determined
        suppressMessages(
          poipo_on <- sf::st_as_sf(
            sf::st_sample(x = polpo, size = simulate_points, type = simulate_strategy)
          )
        )
      }

      names(poipo_on) <- "geometry"; sf::st_geometry(poipo_on) <- "geometry"
      poipo_on[[field_occurrence]] <- 1
      poipo_on$x <- sf::st_coordinates(poipo_on)[,1];poipo_on$y <- sf::st_coordinates(poipo_on)[,2]
      poipo <- poipo_on

      # Check that weights are correctly set
      if(length(weight)>1) {
        assertthat::assert_that(length(weight) == nrow(poipo))
      } else {
        weight <- rep(weight, nrow(poipo))
      }

      # Add simulated poipo object instead
      add_biodiversity_poipo(x, poipo = poipo, name = paste0(name, "_simulated"),
                             field_occurrence = field_occurrence, formula = formula, family = family, link = link,
                             weight = weight, ... )
    } else {
      # Convert formula if necessary
      formula = to_formula(formula)

      # Create a new id for this dataset
      id = new_id()

      # Check that weights are correctly set
      if(length(weight)>1) {
        assertthat::assert_that(length(weight) == nrow(polpo))
      } else {
        weight <- rep(weight, nrow(polpo))
      }

      # Define the biodiversity object
      bd <- BiodiversityDataset$new(
        name = ifelse(is.null(name), substring(id[[1]], 1, 8), name),
        id = id,
        equation = formula,
        family = family,
        link = link,
        type = "polpo",
        weight = weight,
        field_occurrence = field_occurrence,
        data = format_biodiversity_data(polpo,field_occurrence),
        use_intercept = separate_intercept,
        pseudoabsence_settings = pseudoabsence_settings
      )

      # Finally set the data to the BiodiversityDistribution object
      y <- x$clone(deep = TRUE)
      y$set_biodiversity(id, bd)
    }
  }
)

#' Add biodiversity polygon dataset to a distribution object (presence-absence)
#'
#' @description This function can be used to add a [`sf`] polygon dataset to an
#' existing distribution object. Presence-absence polygon data assumes that each
#' area within the polygon can be treated as 'presence' for the species, while
#' each area outside the polygon is where the species is absent.
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param polpa A [`sf`] polygon object of presence-absence occurrences.
#' @param name The name of the biodiversity dataset used as internal identifier.
#' @param field_occurrence A [`numeric`] or [`character`] location of biodiversity point records.
#' @param formula A [`character`] or [`formula`] object to be passed.
#' Default (\code{NULL}) is to use all covariates .
#' @param family A [`character`] stating the family to be used (Default: \code{binomial}).
#' @param link A [`character`] to overwrite the default link function (Default: \code{NULL}).
#' @param weight A [`numeric`] value acting as a multiplier with regards to any
#' weights used in the modelling. Larger weights indicate higher weighting
#' relative to any other datasets. By default set to \code{1} if only one
#' dataset is added. A [`vector`] is also supported but must be of the same
#' length as \code{"polpa"}.
#' @param simulate Simulate poipa points within its boundaries. Result are
#' passed to [`add_biodiversity_poipa`] (Default: \code{FALSE}).
#' @param simulate_points A [`numeric`] number of points to be created by
#' simulation.
#' @param simulate_bias A [`SpatRaster`] layer describing an eventual preference
#' for simulation (Default: \code{NULL}).
#' @param simulate_strategy A [`character`] stating the strategy for sampling.
#' Can be set to either. \code{'random'} or \code{'regular'}, the latter
#' requiring a raster supplied in the \code{'simulate_weights'} parameter.
#' @param separate_intercept A [`logical`] value stating whether a separate
#' intercept is to be added in shared likelihood models for engines
#' [engine_inla], [engine_inlabru] and [engine_stan].
#' @param docheck [`logical`] on whether additional checks should be performed
#' (e.g. intersection tests) (Default: \code{TRUE}).
#' @param pseudoabsence_settings Either \code{NULL} or a
#' [`pseudoabs_settings()`] created settings object.
#' @param ... Other parameters passed down.
#'
#' @details The default approach for polygon data is to sample presence-absence
#' points across the region of the polygons. This function thus adds as a
#' wrapper to [`add_biodiversity_poipa()`] as presence-only points are created
#' by the model. Note if the polygon is used directly in the modelling the
#' link between covariates and polygonal data is established by regular
#' sampling of points within the polygon and is thus equivalent to simulating
#' the points directly.
#'
#' For an integration of range data as predictor or offset, see [`add_predictor_range()`]
#' and [`add_offset_range()`] instead.
#'
#' @returns Adds biodiversity data to [distribution] object.
#'
#' @family add_biodiversity
#' @keywords biodiversity
#'
#' @examples
#' \dontrun{
#'  x <- distribution(background) |>
#'    add_biodiversity_polpa(protectedArea)
#' }
#' @name add_biodiversity_polpa
NULL

#' @rdname add_biodiversity_polpa
#' @export
methods::setGeneric(
  "add_biodiversity_polpa",
  signature = methods::signature("x", "polpa"),
  function(x, polpa, name = NULL, field_occurrence = "observed", formula = NULL, family = "binomial", link = NULL,
           weight = 1, simulate = FALSE, simulate_points = 100, simulate_bias = NULL, simulate_strategy = "random",
           separate_intercept = TRUE, docheck = TRUE, pseudoabsence_settings = NULL, ...) standardGeneric("add_biodiversity_polpa"))

#' @rdname add_biodiversity_polpa
methods::setMethod(
  "add_biodiversity_polpa",
  methods::signature(x = "BiodiversityDistribution", polpa = "sf"),
  function(x, polpa, name = NULL, field_occurrence = "observed", formula = NULL, family = "binomial", link = NULL,
           weight = 1, simulate = FALSE, simulate_points = 100, simulate_bias = NULL, simulate_strategy = "random",
           separate_intercept = TRUE, docheck = TRUE, pseudoabsence_settings = NULL, ... ) {
    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            inherits(polpa, "Spatial") || inherits(polpa, "sf") || inherits(polpa, "data.frame") || inherits(polpa, "tibble"),
                            assertthat::is.scalar(field_occurrence), assertthat::has_name(polpa, field_occurrence),
                            inherits(formula, "formula") || is.null(formula) || is.character(formula),
                            is.character(family),
                            is.null(link) || is.character(link),
                            assertthat::is.flag(simulate), is.numeric(simulate_points),
                            is.null(simulate_bias) || inherits(simulate_bias, "SpatRaster"),
                            is.numeric(weight) && all(weight > 0),
                            is.logical(separate_intercept)
    )

    # Check type and ensure that is a polygon
    assertthat::assert_that(all( unique( sf::st_geometry_type(polpa) ) %in% c("POLYGON","MULTIPOLYGON") ),
                            msg = "This method works for spatial data of type polygon only.")

    assertthat::assert_that(length(unique(polpa[[field_occurrence]])) <= 2,
                            msg = "More 2 unique values. Specify a column.")

    # Transform to background for analysis
    if(sf::st_crs(x$background) != sf::st_crs(polpa)){
      polpa <- polpa |> sf::st_transform(crs = sf::st_crs(x$background))
    }

    # Record presence absence to 0 and 1
    if(is.factor(polpa[[field_occurrence]])) {
      polpa[[field_occurrence]] <- as.numeric(as.character(polpa[[field_occurrence]]))
    }

    # Simulate presence absence points rather than using the range directly
    if(simulate){
      # Simulation strategy
      simulate_strategy <- match.arg(simulate_strategy, c('random', 'regular'), several.ok = FALSE)

      if(!is.null(simulate_bias)){
        # Crop to target range
        simulate_bias <- terra::crop(simulate_bias, polpa)
        # Normalize the weight layer if is not a factorized, else set everything to 1
        if(is.null(levels(simulate_bias))) simulate_bias <- predictor_transform(simulate_bias, "norm") else simulate_bias[simulate_bias>0] <- 1

        # Weighted sampling on background raster, the greater the value, the more likely sampled points
        ptscell <- sample(which(!is.na(simulate_bias[])),
                          size = simulate_points,
                          prob = simulate_bias[which(!is.na(simulate_bias[]))],
                          replace = TRUE)
        poipa_on <- terra::as.data.frame(terra::xyFromCell(simulate_bias, ptscell))
        poipa_on <- sf::st_as_sf(poipa_on, coords = c("x","y"), crs = sf::st_crs(simulate_bias))

      } else {
        # Simply sample presence points within as determined
        suppressMessages(
          poipa_on <- sf::st_as_sf(
            sf::st_sample(x = polpa, size = simulate_points, type = simulate_strategy)
          )
        )
      }

      names(poipa_on) <- "geometry"; sf::st_geometry(poipa_on) <- "geometry"
      poipa_on[[field_occurrence]] <- 1
      poipa_on$x <- sf::st_coordinates(poipa_on)[,1]; poipa_on$y <- sf::st_coordinates(poipa_on)[,2]

      # Get absence data for poipa data, masking out the range of first
      if(is.Raster(x$background)){
        bg_masked <- terra::mask(x$background, polpa, inverse = TRUE)
      } else bg_masked <- x$background

      suppressMessages(
        poipa_off <- sf::st_as_sf(
          sf::st_sample(x = bg_masked, size = simulate_points, type = "random")
        )
      )
      names(poipa_off) <- "geometry"; sf::st_geometry(poipa_off) <- "geometry"

      if(!is.Raster(x$background)){
        # Remove points on the range
        suppressMessages(
          wi <- sf::st_within(poipa_off, polpa, sparse = FALSE)
        )
        poipa_off <- poipa_off[!apply(wi, 1, any),]
      }

      # Remove points on the range
      poipa_off[[field_occurrence]] <- 0
      poipa_off$x <- sf::st_coordinates(poipa_off)[,1];poipa_off$y <- sf::st_coordinates(poipa_off)[,2]

      poipa <- rbind(poipa_on,poipa_off)

      # Check that weights are correctly set
      if(length(weight)>1) {
        assertthat::assert_that(length(weight) == nrow(poipa))
      } else {
        weight <- rep(weight, nrow(poipa))
      }

      # Add simulated poipa object instead
      add_biodiversity_poipa(x, poipa = poipa, name = paste0(name, "_simulated"),
                             field_occurrence = field_occurrence, formula = formula, family = family, link = link,
                             weight = weight, ... )
    } else {

      # Messenger
      if(getOption("ibis.setupmessages", default = TRUE)) myLog("[Setup]","green","Adding polpa dataset...")

      # If no points are simulated, ensure that the polygon has objects with at
      # least 2 factor levels
      assertthat::assert_that(
        length(unique(polpa[[field_occurrence]]))==2
      )
      # Convert formula if necessary
      formula = to_formula(formula)

      # Create a new id for this dataset
      id = new_id()

      # Check that weights are correctly set
      if(length(weight)>1) {
        assertthat::assert_that(length(weight) == nrow(polpa))
      } else {
        weight <- rep(weight, nrow(polpa))
      }

      # Define the biodiversity object
      bd <- BiodiversityDataset$new(
        name = ifelse(is.null(name), substring(id[[1]], 1, 8), name),
        id = id,
        equation = formula,
        family = family,
        link = link,
        type = "polpa",
        weight = weight,
        field_occurrence = field_occurrence,
        data = format_biodiversity_data(polpa, field_occurrence),
        use_intercept = separate_intercept,
        pseudoabsence_settings = pseudoabsence_settings
      )

      # Finally set the data to the BiodiversityDistribution object
      y <- x$clone(deep = TRUE)
      y$set_biodiversity(id, bd)
    }
  }
)
# Removal function ----

#' Remove specific BiodiversityDataset from a [distribution] object
#'
#' @description Remove a particular dataset (or all) from an [distribution] object with
#' a [`BiodiversityDatasetCollection-class`].
#'
#' @param x [distribution()] (i.e. [`BiodiversityDistribution-class`]) object.
#' @param name A [`character`] with the name of the biodiversity dataset.
#' @param id A [`character`] with the id of the biodiversity dataset.
#'
#' @examples
#' \dontrun{
#' distribution(background) |>
#'  add_biodiversity_poipa(species, "Duckus communus")
#'  rm_biodiversity(names = "Duckus communus")
#' }
#'
#' @name rm_biodiversity
NULL

#' @rdname rm_biodiversity
#' @export
methods::setGeneric(
  "rm_biodiversity",
  signature = methods::signature("x"),
  function(x, name, id) standardGeneric("rm_biodiversity"))

#' @rdname rm_biodiversity
methods::setMethod(
  "rm_biodiversity",
  methods::signature(x = "BiodiversityDistribution"),
  function(x, name, id) {
    assertthat::assert_that(
      inherits(x, "BiodiversityDistribution")
    )
    # Make a deep copy
    y <- x$clone(deep = TRUE)

    # If both are missing, get all ids
    if(missing(name) && missing(id)){
      y$biodiversity <- BiodiversityDatasetCollection$new()
      return(y)
    } else {
      # Get id of name
      if(!missing(name)){
        n <- x$get_biodiversity_names()
        if((name %in% n)) id <- names(n)[which(name %in% n)]
      }
      # If id is set
      if(!missing(id)){
        n <- x$get_biodiversity_ids()
        if(!(id %in% n)) id <- NULL
      }
      # Is there anything to remove?
      assertthat::assert_that(is.character(id) || is.Id(id),
                              msg = "Provided biodiversity dataset not found in object!")

      # Remove dataset
      y$biodiversity <- y$biodiversity$rm_data(id)
      return(y)
    }
  }
)

#' Format biodiversity dataset to standardized format
#'
#' @param x A [`data.frame`] or [`sf`] object of biodiversity information.
#' @param field_occurrence A [`numeric`] or [`character`] location of biodiversity records.
#' @param field_space A [`vector`] on the column names (Default: \code{'x'}, \code{'y'}).
#' @param ... Other parameters passed down.
#'
#' @noRd
#'
#' @keywords internal
format_biodiversity_data <- function(x, field_occurrence, field_space = c("x","y"), ... ){
  # Final checks
  if(inherits(x,'sf')) assertthat::assert_that(unique(sf::st_geometry_type(x)) %in% c("POINT","MULTIPOINT",
                                                                                      "POLYGON","MULTIPOLYGON"))

  # If data.frame or tibble, check whether coordinates are in there
  if(!(inherits(x,'sf') )){
    # Check whether field_space columns are present
    if(all(assertthat::has_name(x, field_space))){
      # Convert data.frame to spatial format
      x <- sf::st_as_sf(x, coords = field_space)
    } else {
      x <- guess_sf(x)
      # Spatial column suggestions
      assertthat::assert_that(
        inherits(x, "sf"),
        msg = 'Data could not be converted to spatial format. Specify manually or set to [x] and [y].'
      )
      # Add field_space columns if not already existing
      if(!all(assertthat::has_name(x, field_space))){
        x[[field_space[1]]] <- sf::st_coordinates(x)[,1]
        x[[field_space[2]]] <- sf::st_coordinates(x)[,2]
      }
    }
    assertthat::assert_that( all(assertthat::has_name(x, field_space)),
                             msg ='No spatial column found in the dataset. Specify manually or set to [x] and [y].')
    # Select and format
    out <- subset(x, select = c(field_space, field_occurrence) ) |>
      tibble::as_tibble()
  } else {
    if(inherits(x, "Spatial")) x <- sf::st_as_sf(x) # First convert to sf

    if(unique(sf::st_geometry_type(x)) %in% c("POINT","MULTIPOINT")){
      # Take target column and append coordinates to it
      out <- cbind(subset(x, select = field_occurrence),
                   sf::st_coordinates(x)) |> tibble::as_tibble()
    } else if(unique(sf::st_geometry_type(x)) %in% c('POLYGON','MULTIPOLYGON')){
      # Return target column and spatial object as such
      out <- subset(x, select = field_occurrence)
    }
  }
  out
}

