#' @include utils.R bdproto.R bdproto-biodiversitydistribution.R bdproto-predictors.R bdproto-biodiversityscenario.R
NULL

#' Add GLOBIOM-DownScaleR derived predictors to a Biodiversity distribution
#' object
#'
#' @description This is a customized function to format and add downscaled
#' land-use shares from the [Global Biosphere Management Model
#' (GLOBIOM)](https://iiasa.github.io/GLOBIOM/) to a [distribution] or
#' [BiodiversityScenario] in ibis.iSDM. GLOBIOM is a partial-equilibrium model
#' developed at IIASA and represents land-use sectors with a rich set of
#' environmental and socio-economic parameters, where for instance the
#' agricultural and forestry sector are estimated through dedicated
#' process-based models. GLOBIOM outputs are spatial explicit and usually at a
#' half-degree resolution globally. For finer grain analyses GLOBIOM outputs can
#' be produced in a downscaled format with a customized statistical [downscaling
#' module](https://github.com/iiasa/DownScale).
#'
#' The purpose of this script is to format the GLOBIOM outputs of *DownScale*
#' for the use in the ibis.iSDM package.
#'
#' @param x A [`BiodiversityDistribution-class`] or [`BiodiversityScenario-class`] object.
#' @param fname A [`character`] pointing to a netCDF with the GLOBIOM data.
#' @param names A [`vector`] of character names describing the environmental
#' stack in case they should be renamed (Default: \code{NULL}).
#' @param transform A [`vector`] stating whether predictors should be preprocessed
#' in any way (Options: \code{'none'},\code{'pca'}, \code{'scale'}, \code{'norm'})
#' @param derivates A Boolean check whether derivate features should be considered
#' (Options: \code{'none'}, \code{'thresh'}, \code{'hinge'}, \code{'quad'}) )
#' @param derivate_knots A single [`numeric`] or [`vector`] giving the number of
#' knots for derivate creation if relevant (Default: \code{4}).
#' @param int_variables A [`vector`] with length greater or equal than \code{2}
#' specifying the covariates (Default: \code{NULL}).
#' @param bgmask Check whether the environmental data should be masked with the
#' background layer (Default: \code{TRUE})
#' @param harmonize_na A [`logical`] value indicating of whether NA values should
#' be harmonized among predictors (Default: \code{FALSE})
#' @param priors A [`PriorList-class`] object. Default is set to \code{NULL}
#' which uses default prior assumptions.
#' @param ... Other parameters passed down
#'
#' @details See [`add_predictors()`] for additional parameters and
#' customizations. For more (manual) control the function for formatting the
#' GLOBIOM data can also be called directly via `formatGLOBIOM()`.
#'
#' @seealso [add_predictors]
#'
#' @examples
#' \dontrun{
#'  obj <- distribution(background) |>
#'         add_predictors_globiom(fname = "", transform = 'none')
#'  obj
#' }
#'
#' @name add_predictors_globiom
NULL

#' @rdname add_predictors_globiom
#' @export
methods::setGeneric(
  "add_predictors_globiom",
  signature = methods::signature("x", "fname"),
  function(x, fname, names = NULL, transform = 'none', derivates = 'none', derivate_knots = 4, int_variables = NULL,
           bgmask = TRUE, harmonize_na = FALSE,
           priors = NULL, ...) standardGeneric("add_predictors_globiom"))

#' @rdname add_predictors_globiom
methods::setMethod(
  "add_predictors_globiom",
  methods::signature(x = "BiodiversityDistribution", fname = "character"),
  function(x, fname, names = NULL, transform = 'none', derivates = 'none', derivate_knots = 4, int_variables = NULL,
           bgmask = TRUE, harmonize_na = FALSE, priors = NULL, ... ) {
    # Try and match transform and derivatives arguments
    transform <- match.arg(transform, c('none','pca', 'scale', 'norm', 'windsor'), several.ok = TRUE)
    derivates <- match.arg(derivates, c('none','thresh', 'hinge', 'quadratic', 'bin'), several.ok = TRUE)

    # Check that file exists and has the correct endings
    assertthat::assert_that(is.character(fname),
                            file.exists(fname),
                            assertthat::is.readable(fname),
                            assertthat::has_extension(fname, "nc"),
                            msg = "The provided path to GLOBIOM land-use shares could not be found or is not readable!"
    )

    assertthat::assert_that(inherits(x, "BiodiversityDistribution"),
                            is.null(names) || assertthat::is.scalar(names) || is.vector(names),
                            is.null(priors) || inherits(priors,'PriorList'),
                            is.vector(derivate_knots) || is.null(derivate_knots),
                            is.null(int_variables) || is.vector(int_variables)
    )

    # Messenger
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Formatting GLOBIOM inputs for species distribution modelling.')

    # Get and format the GLOBIOM data
    env <- formatGLOBIOM(fname = fname,
                         oftype = "raster",
                         period = "reference",
                         template = x$background
    )

    if(is.list(env)) env <- env[[1]] # Take the first reference entry
    assertthat::assert_that(is.Raster(env),
                            terra::nlyr(env)>0)

    if(!is.null(names)) {
      assertthat::assert_that(terra::nlyr(env)==length(names),
                              all(is.character(names)),
                              msg = 'Provided names not of same length as environmental data.')
      # Set names of env
      names(env) <- names
    }

    # Check that all names allowed
    problematic_names <- grep("offset|w|weight|spatial_offset|Intercept|spatial.field", names(env),fixed = TRUE)
    if( length(problematic_names)>0 ){
      stop(paste0("Some predictor names are not allowed as they might interfere with model fitting:", paste0(names(env)[problematic_names],collapse = " | ")))
    }

    # If priors have been set, save them in the distribution object
    if(!is.null(priors)) {
      assertthat::assert_that( all( priors$varnames() %in% names(env) ) )
      x <- x$set_priors(priors)
    }
    # Harmonize NA values
    if(harmonize_na){
      if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Harmonizing missing values...')
      env <- predictor_homogenize_na(env, fill = FALSE)
    }

    # Standardization and scaling
    if('none' %notin% transform){
      if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Transforming predictors...')
      for(tt in transform) env <- predictor_transform(env, option = tt)
    }

    # Calculate derivates if set
    if('none' %notin% derivates){
      if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Creating predictor derivates...')
      new_env <- terra::rast()
      for(dd in derivates) {
        suppressWarnings(
          new_env <- c(new_env, predictor_derivate(env, option = dd, nknots = derivate_knots, int_variables = int_variables) )
        )
      }
      # Add to env
      env <- c(env, new_env)
    }

    # Generally not relevant for GLOBIOM unless created as derivate
    attr(env, 'has_factors') <- FALSE

    # Assign an attribute to this object to keep track of it
    attr(env,'transform') <- transform

    # Mask predictors with existing background layer
    if(bgmask){
      env <- terra::mask(env, mask = x$background)
      env <- terra::rast(env)
    }

    # Check whether predictors already exist, if so overwrite
    if(!is.Waiver(x$predictors)) myLog('[Setup]','yellow','Overwriting existing predictors.')

    # Sanitize names if specified
    if(getOption('ibis.cleannames')) names(env) <- sanitize_names(names(env))

    # Finally set the data to the BiodiversityDistribution object
    x$set_predictors(
      bdproto(NULL, PredictorDataset,
              id = new_id(),
              data = env,
              ...
      )
    )
  }
)

#' @rdname add_predictors_globiom
methods::setMethod(
  "add_predictors_globiom",
  methods::signature(x = "BiodiversityScenario", fname = "character"),
  function(x, fname, names = NULL, transform = 'none', derivates = 'none', derivate_knots = 4, int_variables = NULL,
           harmonize_na = FALSE, ... ) {
    # Try and match transform and derivatives arguments
    transform <- match.arg(transform, c('none','pca', 'scale', 'norm', 'windsor') , several.ok = TRUE)
    derivates <- match.arg(derivates, c('none','thresh', 'hinge', 'quadratic', 'bin') , several.ok = TRUE)

    # Check that file exists and has the correct endings
    assertthat::assert_that(is.character(fname),
                            file.exists(fname),
                            assertthat::is.readable(fname),
                            assertthat::has_extension(fname, "nc"),
                            msg = "The provided path to GLOBIOM land-use shares could not be found or is not readable!"
    )
    assertthat::assert_that(inherits(x, "BiodiversityScenario"),
                            is.null(names) || assertthat::is.scalar(names) || is.vector(names),
                            is.logical(harmonize_na),
                            is.vector(derivate_knots) || is.null(derivate_knots),
                            is.null(int_variables) || is.vector(int_variables)
    )

    # Get model object
    obj <- x$get_model()

    # Messenger
    if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Adding GLOBIOM predictors to scenario object...')

    # Get and format the GLOBIOM data
    env <- formatGLOBIOM(fname = fname,
                         oftype = "stars",
                         period = "projection",
                         template = obj$model$background
    )
    assertthat::assert_that( inherits(env, "stars") )

    # Rename attributes if names is specified
    if(!is.null(names)){
      assertthat::assert_that(length(names) == length(env))
      names(env) <- names
    }

    # Harmonize NA values
    if(harmonize_na){
      stop('Missing data harmonization for stars not yet implemented!') #TODO
      if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Harmonizing missing values...')
      env <- predictor_homogenize_na(env, fill = FALSE)
    }

    # Standardization and scaling
    if('none' %notin% transform){
      if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Transforming predictors...')
      for(tt in transform) env <- predictor_transform(env, option = tt)
    }

    # # Calculate derivates if set
    if('none' %notin% derivates){
      # Get variable names
      varn <- obj$get_coefficients()[['Feature']]
      # Are there any derivates present in the coefficients?
      if(any( length( grep("hinge__|bin__|quad__|thresh__", varn ) ) > 0 )){
        if(getOption('ibis.setupmessages')) myLog('[Setup]','green','Creating predictor derivates...')
        for(dd in derivates){
          if(any(grep(dd, varn))){
            env <- predictor_derivate(env, option = dd, nknots = derivate_knots, int_variables = int_variables, deriv = varn)
          } else {
            if(getOption('ibis.setupmessages')) myLog('[Setup]','red', paste0(derivates,' derivates should be created, but not found among coefficients!'))
          }
        }
      } else {
        if(getOption('ibis.setupmessages')) myLog('[Setup]','red','No derivates found among coefficients. None created for projection!')
      }
    }

    # Get and format Time period
    env_dim <- stars::st_dimensions(env)
    timeperiod <- stars::st_get_dimension_values(env, "time", center = TRUE)
    if(is.numeric(timeperiod)){
      # Format to Posix. Assuming years only
      timeperiod <- as.POSIXct(paste0(timeperiod,"-01-01"))
    }
    if(anyNA(timeperiod)) stop('Third dimension is not a time value!')

    # Check whether predictors already exist, if so overwrite
    if(!is.Waiver(x$predictors)) myLog('[Setup]','yellow','Overwriting existing predictors.')

    # Sanitize names if specified
    if(getOption('ibis.cleannames')) names(env) <- sanitize_names(names(env))

    # Finally set the data to the BiodiversityScenario object
    x$set_predictors(
      bdproto(NULL, PredictorDataset,
              id = new_id(),
              data = env,
              timeperiod = timeperiod,
              ...
      )
    )
  }
)

#' Function to format a prepared GLOBIOM netCDF file for use in Ibis.iSDM
#'
#' @description This function expects a downscaled GLOBIOM output as created in
#' the BIOCLIMA project. Likely of little use for anyone outside IIASA.
#'
#' @param fname A filename in [`character`] pointing to a GLOBIOM output in netCDF format.
#' @param oftype A [`character`] denoting the output type (Default: \code{'raster'}).
#' @param ignore A [`vector`] of variables to be ignored (Default: \code{NULL}).
#' @param period A [`character`] limiting the period to be returned from the
#' formatted data. Options include \code{"reference"} for the first entry, \code{"projection"}
#' for all entries but the first, and \code{"all"} for all entries (Default: \code{"reference"}).
#' @param template An optional [`SpatRaster`] object towards which projects
#' should be transformed.
#' @param shares_to_area A [`logical`] on whether shares should be corrected to
#' areas (if identified).
#' @param use_gdalutils (Deprecated) [`logical`] on to use gdalutils hack-around.
#' @param verbose [`logical`] on whether to be chatty.
#'
#' @return A [`SpatRaster`] stack with the formatted GLOBIOM predictors.
#'
#' @keywords internal, utils
#'
#' @examples
#' \dontrun{
#' # Expects a filename pointing to a netCDF file.
#' covariates <- formatGLOBIOM(fname)
#' }
#'
#' @export
formatGLOBIOM <- function(fname, oftype = "raster", ignore = NULL,
                          period = "all", template = NULL, shares_to_area = FALSE,
                          use_gdalutils = FALSE,
                          verbose = getOption("ibis.setupmessages")){
  assertthat::assert_that(
    file.exists(fname),
    assertthat::has_extension(fname, "nc"),
    is.character(oftype),
    is.null(ignore) || is.character(ignore),
    is.character(period),
    is.character(fname),
    is.logical(shares_to_area),
    is.logical(use_gdalutils),
    is.logical(verbose)
  )
  period <- match.arg(period, c("reference", "projection", "all"), several.ok = FALSE)
  check_package("stars")
  check_package("dplyr")
  check_package("cubelyr")
  check_package("ncdf4")

  # Try and load in the GLOBIOM file to get the attributes
  fatt <- ncdf4::nc_open(fname)
  if(verbose) myLog('[Setup]','green',"Found ", fatt$ndims, " dimensions and ", fatt$nvars, " variables")

  # Get all dimension names and variable names
  dims <- names(fatt$dim)
  vars <- names(fatt$var)
  if(!is.null(ignore)) assertthat::assert_that( all( ignore %in% vars ) )

  attrs <- list() # For storing the attributes
  sc <- vector() # For storing the scenario files
  sc_area <- new_waiver() # For storing any area information if set

  # Now open the netcdf file with stars
  if( length( grep("netcdf", stars::detect.driver(fname), ignore.case = TRUE) )>0 ){
    if(verbose){
      myLog('[Predictor]','green',"Loading in predictor file...")
      pb <- progress::progress_bar$new(total = length(vars),
                                       format = "Loading :variable (:spin) [:bar] :percent")
    }

    for(v in vars) {
      if(verbose) pb$tick(tokens = list(variable = v))
      if(!is.null(ignore)) if(ignore == v) next()

      # Get and save the attributes of each variable
      attrs[[v]] <- ncdf4::ncatt_get(fatt, varid = v, verbose = FALSE)

      # Load in the variable
      suppressWarnings(
        suppressMessages(
          ff <- stars::read_ncdf(fname,
                                 var = v,
                                 proxy = FALSE,
                                 make_time = TRUE, # Make time on 'time' band
                                 make_units = FALSE # To avoid unnecessary errors due to unknown units
          )
        )
      )

      # Sometimes variables don't seem to have a time dimension
      if(!"time" %in% names(stars::st_dimensions(ff))) {
        if(shares_to_area && length(grep("area",names(ff)))>0){
          # Check that the unit is a unit
          if(fatt$var[[v]]$units %in% c("km2","ha","m2")){
            sc_area <- ff
          }
        } else {
          next()
        }
      }

      # Crop to background extent if set
      # if(!is.null(template)){
      # FIXME: Currently this code, while working clips too much of Europe.
      # Likely need to
      # bbox <- sf::st_bbox(template) |> sf::st_as_sfc() |>
      #   sf::st_transform(crs = sf::st_crs(ff))
      # suppressMessages(
      # ff <- ff |> stars:::st_crop.stars(bbox)
      # )
      # }

      # Record dimensions for later
      full_dis <- stars::st_dimensions(ff)

      # Get dimensions other that x,y and time and split
      # Commonly used column names
      check = c("x","X","lon","longitude", "y", "Y", "lat", "latitude", "time", "Time", "year", "Year")
      chk <- which(!names(stars::st_dimensions(ff)) %in% check)

      if(length(chk)>0){
        for(i in chk){
          col_class <- names(stars::st_dimensions(ff))[i]
          # FIXME: Dirty hack to remove forest zoning
          if(length( grep("zone",col_class,ignore.case = T) )>0) next()

          # And class units as description from over
          class_units <- fatt$dim[[col_class]]$units
          class_units <-  class_units |>
            strsplit(";") |>
            # Remove emptyspace and special symbols
            sapply(function(y)  gsub("[^0-9A-Za-z///' ]", "" , y, ignore.case = TRUE) ) |>
            sapply(function(y)  gsub(" ", "" , y, ignore.case = TRUE) )
          # Convert to vector and make names
          class_units <- paste0(
            v, "__",
            make.names(unlist(class_units)) |> as.vector()
          )

          ff <- ff |> stars:::split.stars(col_class) |> stats::setNames(nm = class_units)

          # FIXME: Dirty hack to deal with the forest zone dimension
          # If there are more dimensions than 3, aggregate over them
          if( length(stars::st_dimensions(ff)) >3){
            # Aggregate spatial-temporally
            ff <- stars::st_apply(ff, c("longitude", "latitude", "time"), sum, na.rm = TRUE)
          }
        }
      }

      # Finally aggregate
      if(!is.null(template) && is.Raster(template)){
        # FIXME: MJ 14/11/2022 - The code below is buggy, resulting in odd
        # curvilinear extrapolations for Europe Hacky approach now is to convert
        # to raster, crop, project and then convert back. Only use if gdalUtils
        # is installed
        if(("gdalUtilities" %in% utils::installed.packages()[,1])&&use_gdalutils){
          ff <- hack_project_stars(ff, template, use_gdalutils)
        } else {
          # Make background
          bg <- stars::st_as_stars(template)

          # # Get resolution
          res <- stars::st_res(bg)
          assertthat::assert_that(!anyNA(res))

          # # And warp by projecting and resampling
          ff <- ff |> stars::st_warp(bg, crs = sf::st_crs(bg),
                                     cellsize = res, method = "near") |>
            sf::st_transform(crs = sf::st_crs(template))
        }
        # Overwrite full dimensions
        full_dis <- stars::st_dimensions(ff)
      }
      # Now append to vector
      sc <- c(sc, ff)
      rm(ff)
    }
    invisible(gc())
    assertthat::assert_that(length(names(full_dis))>=3)

    # Format sc object as stars and set dimensions again
    sc <- stars::st_as_stars(sc)
    assertthat::assert_that(length(sc)>0)
    full_dis <- full_dis[c(
      grep("x|longitude",names(full_dis), ignore.case = TRUE,value = TRUE),
      grep("y|latitude",names(full_dis), ignore.case = TRUE,value = TRUE),
      grep("year|time",names(full_dis), ignore.case = TRUE,value = TRUE)
    )] # Order assumed to be correct
    assertthat::assert_that(length(names(full_dis))==3)
    stars::st_dimensions(sc) <- full_dis # Target dimensions

  } else { stop("Fileformat not recognized!")}

  # Get time dimension (without applying offset) so at the centre
  times <- stars::st_get_dimension_values(sc, "time", center = TRUE)

  # Make checks on length of times and if equal to one, drop. check.
  if(length(times)==1){
    if(period == "projection") stop("Found only a single time slot. Projections not possible.")
    if(verbose) myLog('[Setup]','yellow','Found only a single time point in file. Dropping time dimension.')
    # Drop the time dimension
    sc <- stars:::adrop.stars(sc, drop = which(names(stars::st_dimensions(sc)) == "time") )
  }

  # Formate times unit and convert to posix if not already set
  if(is.numeric(times) && length(times) > 1){
    # Assume year and paste0 as properly POSIX formatted
    times <- as.POSIXct( paste0(times, "-01-01") )
    sc <- stars::st_set_dimensions(sc, "time", times)
  }

  # Depending on the period, slice the input data
  if(period == "reference"){
    # Get the first entry and filter
    if(length(times)>1){
      # In case times got removed
      times_first <- stars::st_get_dimension_values(sc, "time")[1]
      sc <- sc |> stars:::filter.stars(time == times_first)
      times <- times_first;rm(times_first)
    }
  } else if(period == "projection"){
    # Remove the first time entry instead, only using the last entries
    times_allbutfirst <- stars::st_get_dimension_values(sc, "time")[-1]
    sc <- sc |> stars:::filter.stars(time %in% times_allbutfirst)
    times <- times_allbutfirst; rm(times_allbutfirst)
  }
  assertthat::assert_that(length(times)>0,
                          length(sc)>=1)

  # Create raster template if set
  if(!is.null(template)){
    # Check that template is a raster, otherwise rasterize for GLOBIOM use
    if(inherits(template, "sf")){
      o <- sc |> stars:::slice.stars("time" , 1) |> terra::rast()
      template <- terra::rasterize(template, o, field = 1)
      rm(o)
    }
  }

  # Correct shares to area if set
  if(shares_to_area && inherits(sc_area,"stars")){
    # Transform and warp the shares
    sc_area <- stars::st_warp(sc_area, stars::st_as_stars(template), crs = sf::st_crs(sc),method = "near")
    # grep those layers with the name share
    shares <- grep(pattern = "share|fraction|proportion", names(sc),value = TRUE)
    sc[shares] <- sc[shares] * sc_area
  }

  # Now format outputs depending on type, either returning the raster or the stars object
  if(oftype == "raster"){
    # Output type raster, use function from utils_scenario
    out <- stars_to_raster(sc, which = NULL, template = template)
    return(out)
  } else {
    return( sc )
    }
}
