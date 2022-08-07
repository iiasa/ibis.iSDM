#' Generic function to write spatial outputs
#'
#' @description
#' The `write_output` function is a generic wrapper to writing any output files (e.g. projections) created with
#' the [`ibis.iSDM-package`]. It is possible to write outputs of fitted [`DistributionModel`],
#' [`BiodiversityScenario`] or individual [`Raster`] or [`stars`] objects. In case a [`data.frame`]
#' is supplied, the output is written as csv file.
#' **For creating summaries of distribution and scenario parameters and performance, see `write_summary()`**
#' @note
#' By default output files will be overwritten if already existing!
#' @param mod Provided [`DistributionModel`], [`BiodiversityScenario`], [`Raster`] or [`stars`] object.
#' @param fname A [`character`] depicting an output filename.
#' @param dt A [`character`] for the output datatype. Following the [`raster::dataType()`] options (Default: \code{'FLT4S'}).
#' @param verbose [`logical`] indicating whether messages should be shown. Overwrites `getOption("ibis.setupmessages")` (Default: \code{TRUE}).
#' @param ... Any other arguements passed on the individual functions.
#' @returns No R-output is created. A file is written to the target direction.
#' @examples \dontrun{
#' x <- distribution(background) %>%
#'  add_biodiversity_poipo(virtual_points, field_occurrence = 'Observed', name = 'Virtual points') %>%
#'  add_predictors(pred_current, transform = 'scale',derivates = 'none') %>%
#'  engine_xgboost(nrounds = 2000) %>%  train(varsel = FALSE, only_linear = TRUE)
#' write_output(x, "testmodel.tif")
#' }

#' @name write_output
#' @aliases write_output
#' @keywords utils
#' @exportMethod write_output
#' @export
NULL
methods::setGeneric("write_output",
                    signature = methods::signature("mod"),
                    function(mod, fname, dt = "FLT4S", verbose = getOption("ibis.setupmessages"), ...) standardGeneric("write_output"))

#' @name write_output
#' @rdname write_output
#' @usage \S4method{write_output}{ANY, character, character, logical}(mod, fname, dt, verbose)
methods::setMethod(
  "write_output",
  methods::signature("ANY"),
  function(mod, fname, dt = "FLT4S", verbose = getOption("ibis.setupmessages"), ...){
    assertthat::assert_that(
      !missing(mod),
      is.character(fname),
      is.character(dt),
      is.logical(verbose)
    )

    if(verbose && getOption('ibis.setupmessages')) myLog('[Output]','green','Saving output(s)...')

  # This function will only capture the distribution model object and will save them separately
  if(any(class(mod) %in% getOption("ibis.engines")) ){
    # FIXME: If errors occur, check harmonization of saving among engines.
    mod$save(fname = fname)
  } else if(is.Raster(mod)){
      if(raster::extension(fname) %in% c('.tif', '.TIF')) {
        writeGeoTiff(file = mod, fname = fname, dt = dt)
      } else if(raster::extension(fname) %in% c('.nc', '.NC', '.ncdf', '.NCDF')){
        writeNetCDF(file = mode, fname = fname, varName = names(mod), dt = dt)
      } else {
        stop("Output type could not be determined. Currently only geoTIFF and netCDF are supported.")
      }
  } else if(is.data.frame(mod)){
    write.csv(x = mod,file = fname,...)
  } else {
    # Check that a save function exists for object
    assertthat::assert_that("save" %in%names(mod),
                            msg = "No method to save the output could be found!")
    # Try a generic save
    mod$save(fname = fname)
  }
  invisible()
  }
)

#' @name write_output
#' @rdname write_output
#' @usage \S4method{write_output}{BiodiversityScenario, character, character, logical}(mod, fname, dt, verbose)
methods::setMethod(
  "write_output",
  methods::signature(mod = "BiodiversityScenario"),
  function(mod, fname, dt = "FLT4S", verbose = getOption("ibis.setupmessages"), ...) {
    assertthat::assert_that(
      !missing(mod),
      is.character(fname),
      is.character(dt),
      is.logical(verbose)
    )
    # Get outputs
    mod$save(fname = fname, type = tools::file_ext(fname), dt = dt)
    invisible()
  }
)

#' @name write_output
#' @rdname write_output
#' @usage \S4method{write_output}{RasterLayer, character, character, logical}(mod, fname, dt, verbose)
methods::setMethod(
  "write_output",
  methods::signature(mod = "RasterLayer"),
  function(mod, fname, dt = "FLT4S", verbose = getOption("ibis.setupmessages"), ...) {
    assertthat::assert_that(
      !missing(mod),
      is.Raster(mod),
      is.character(fname),
      is.character(dt),
      is.logical(verbose)
    )

    # Write output depending on type
    if(raster::extension(fname) %in% c('.tif', '.TIF')) {
      writeGeoTiff(file = mod, fname = fname, dt = dt)
    } else if(raster::extension(fname) %in% c('.nc', '.NC', '.ncdf', '.NCDF')){
      writeNetCDF(file = mode, fname = fname, varName = names(mod), dt = dt)
    } else {
      stop("Output type could not be determined. Currently only geoTIFF and netCDF are supported.")
    }
    invisible()
  }
)

#' @name write_output
#' @rdname write_output
#' @usage \S4method{write_output}{RasterStack, character, character, logical}(mod, fname, dt, verbose)
methods::setMethod(
  "write_output",
  methods::signature(mod = "RasterStack"),
  function(mod, fname, dt = "FLT4S", verbose = getOption("ibis.setupmessages"),...) {
    assertthat::assert_that(
      !missing(mod),
      is.Raster(mod),
      is.character(fname),
      is.character(dt),
      is.logical(verbose)
    )

    # Write output depending on type
    if(raster::extension(fname) %in% c('.tif', '.TIF')) {
      writeGeoTiff(file = mod, fname = fname, dt = dt)
    } else if(raster::extension(fname) %in% c('.nc', '.NC', '.ncdf', '.NCDF')){
      writeNetCDF(file = mode, fname = fname, varName = names(mod), dt = dt)
    } else {
      stop("Output type could not be determined. Currently only geoTIFF and netCDF are supported.")
    }
    invisible()
  }
)

#' @name write_output
#' @rdname write_output
#' @usage \S4method{write_output}{data.frame, character, character, logical}(mod, fname, dt, verbose)
methods::setMethod(
  "write_output",
  methods::signature(mod = "data.frame"),
  function(mod, fname, dt = "FLT4S", verbose = getOption("ibis.setupmessages"),...) {
    assertthat::assert_that(
      !missing(mod),
      is.data.frame(mod),
      is.character(fname),
      is.logical(verbose)
    )

    # data.frames will be written by default as csv files for consistency
    fname <- paste0( tools::file_path_sans_ext(fname), ".csv")
    write.csv(x = mod, file = fname, ...)
    invisible()
  }
)

#' Saves a raster file in Geotiff format
#'
#' @description Functions that acts as a wrapper to [raster::writeRaster].
#' @param file A [`raster`] object to be saved.
#' @param fname A [`character`] stating the output destination.
#' @param dt The datatype to be written (Default: *Float64*)
#' @param varNA The nodata value to be used (Default: \code{-9999}).
#' @keywords utils, internal
#' @noRd
writeGeoTiff <- function(file, fname, dt = "FLT4S", varNA = -9999){
  assertthat::assert_that(
    inherits(file,'Raster') || inherits(file, 'stars'),
    is.character(fname), is.character(dt),
    is.numeric(varNA)
  )
  if(!assertthat::has_extension(fname,"tif")) fname <- paste0(fname,".tif")

  # Check if layer is factor and deratify if so (causes error otherwise)
  if(any(is.factor(file))){
    file <- raster::deratify(file, complete = TRUE)
  }

  # Save output
  writeRaster(file, fname,
              format='GTiff',
              datatype = dt,
              NAflag = varNA,
              options=c("COMPRESS=DEFLATE","PREDICTOR=2","ZLEVEL=9"),
              overwrite= TRUE )
}

#' Save a raster stack to a netcdf file
#'
#' @param file A [`raster`] object to be saved
#' @param fname A [`character`] stating the output destination
#' @param varName Name for the NetCDF export variable.
#' @param varUnit Units for the NetCDF export variable.
#' @param varLong Long name for the NetCDF export variable.
#' @param dt The datatype to be written. Default is Float64
#' @param varNA The nodata value to be used. Default: -9999
#' @keywords utils, internal
#' @noRd
writeNetCDF <- function(file, fname,
                        varName, varUnit = NULL,
                        varLong = NULL, dt = "FLT4S", varNA = -9999) {
  assertthat::assert_that(
    inherits(file,'Raster'),
    is.character(fname), is.character(dt),
    is.numeric(varNA)
  )
  check_package('ncdf4')
  if(!isNamespaceLoaded("ncdf4")) { attachNamespace("ncdf4");requireNamespace('ncdf4') }
  if(!assertthat::has_extension(fname,"nc")) fname <- paste0(fname,".nc")

  # Output NetCDF file
  raster::writeRaster(x = file,
                      filename = fname,format =  "CDF", overwrite = TRUE,
                      varname = ifelse(is.null(varName),'Prediction',varName),
                      varunit = ifelse(is.null(varUnit),'',varUnit),
                      longname = ifelse(is.null(varLong),'',varLong),
                      xname = ifelse(isLonLat(ras), "Longitude","x"),
                      yname = ifelse(isLonLat(ras), "Latitude","y"),
                      zname = "Time",
                      zunit = "Years since 2000-01-01", # FIXME: Load and format date if provided
                      bylayer = FALSE, # Don't save separate layers
                      datatype = dt, NAflag = varNA
  )

  # FIXME: To be defined
  # Set the time variable
  # ncFile <- ncdf4::nc_open(outFile, write=TRUE)
  # ncdf4::ncvar_put(ncFile, "Time", dtInts)
  # ncdf4::nc_close(ncFile)
}

# ------------------------- #
#### Write summary methods ####

#' Generic function to write summary outputs from created models.
#'
#' @description
#' The [`write_summary`] function is a wrapper function to create summaries from fitted [`DistributionModel`] or
#' [`BiodiversityScenario`] objects. This function will extract parameters and statistics about the used data
#' from the input object and writes the output as either \code{'rds'} or \code{'rdata'} file. Alternative, more open file formats
#' are under consideration.
#' @note
#' No predictions or tabular data is saved through this function.
#' Use [`write_output()`] to save those.
#' @param mod Provided [`DistributionModel`] or [`BiodiversityScenario`] object.
#' @param fname A [`character`] depicting an output filename.
#' The suffix determines the file type of the output (Options: \code{'rds'}, \code{'rdata'}).
#' @param partial A [`logical`] value determining whether partial variable contributions should be calculated and added
#' to the model summary. **Note that this can be rather slow** (Default: \code{FALSE}).
#' @param verbose [`logical`] indicating whether messages should be shown. Overwrites `getOption("ibis.setupmessages")` (Default: \code{TRUE}).
#' @param ... Any other arguments passed on the individual functions.
#' @returns No R-output is created. A file is written to the target direction.
#' @examples \dontrun{
#' x <- distribution(background) %>%
#'  add_biodiversity_poipo(virtual_points, field_occurrence = 'Observed', name = 'Virtual points') %>%
#'  add_predictors(pred_current, transform = 'scale',derivates = 'none') %>%
#'  engine_xgboost(nrounds = 2000) %>%  train(varsel = FALSE, only_linear = TRUE)
#' write_summary(x, "testmodel.rds")
#' }
#' @keywords utils

#' @name write_summary
#' @aliases write_summary
#' @exportMethod write_summary
#' @export
NULL
methods::setGeneric("write_summary",
                    signature = methods::signature("mod"),
                    function(mod, fname, partial = FALSE, verbose = getOption("ibis.setupmessages"),...) standardGeneric("write_summary"))

#' @name write_summary
#' @rdname write_summary
#' @usage \S4method{write_summary}{ANY, character, logical, logical}(mod, fname, partial, verbose)
methods::setMethod(
  "write_summary",
  methods::signature(mod = "ANY"),
  function(mod, fname, partial = FALSE, verbose = getOption("ibis.setupmessages"), ...) {
    assertthat::assert_that(
      !missing(mod),
      is.character(fname),
      is.logical(partial),
      is.logical(verbose)
    )
    assertthat::assert_that(
      inherits(mod, "DistributionModel") || inherits(mod, "BiodiversityScenario"),
      msg = "Only objects created through `train()` or `project()` are supported!"
    )
    # Check writeable or not
    assertthat::assert_that(
      assertthat::is.writeable(dirname(fname)),msg = "Given input folder is not writeable!"
    )

    # Get file extension
    ext <- tolower( tools::file_ext(fname) )
    if(ext == "") ext <- "rds" # Assign rds as default
    ext <- match.arg(ext, choices = c("rds", "rdata"), several.ok = FALSE)
    fname <- paste0(tools::file_path_sans_ext(fname), ".", ext)
    if(file.exists(fname) && (verbose && getOption('ibis.setupmessages'))) myLog('[Output]','yellow','Overwriting existing file...')
    assertthat::assert_that(assertthat::is.writeable(dirname(fname)))
    # --- #
    # Gather the statistics and parameters from the provided file
    output <- list()
    if(inherits(mod, "DistributionModel")){

      # Summarize the model object
      model <- mod$model

      # Model input summary in a tibble
      output[["input"]][["extent"]] <- as.matrix( extent( model$background ) )
      output[["input"]][["predictors"]] <- model$predictors_types
      if(!is.Waiver(model$offset)) output[["input"]][["offset"]] <- names(model$offset) else output[["input"]][["offset"]] <- NA
      if(!is.Waiver(model$priors)){
        output[["input"]][["priors"]] <- model$priors$varnames()
      } else output[["input"]][["priors"]] <- NA

      # Go over biodiversity datasets
      o <- data.frame()
      for(i in 1:length(model$biodiversity)){
        o <- rbind(o,
                     data.frame(id = names(model$biodiversity)[i],
                                name = model$biodiversity[[i]]$name,
                                type = model$biodiversity[[i]]$type,
                                family = model$biodiversity[[i]]$family,
                                equation = deparse1(model$biodiversity[[i]]$equation),
                                obs_pres = sum( model$biodiversity[[i]]$observations$observed > 0 ),
                                obs_abs = sum( model$biodiversity[[i]]$observations$observed == 0 ),
                                n_predictors = length( model$biodiversity[[i]]$predictors_names )
                       )
                   )
      }
      output[["input"]][["biodiversity"]] <- o

      # Model parameters in a tibble
      output[["params"]][["id"]] <- as.character(model$id)
      output[["params"]][["runname"]] <- as.character(model$runname)
      output[["params"]][["algorithm"]] <- class(mod)[1]
      output[["params"]][["equation"]] <- mod$get_equation()
      # Collect settings and parameters if existing
      if( !is.Waiver(mod$get_data("params")) ){
        output[["params"]][["params"]] <- mod$get_data("params")
      }
      if( "settings" %in% names(mod) ){
        output[["params"]][["settings"]] <- mod$settings$data
      }

      # Model summary in a tibble and formula
      output[["output"]][["summary"]] <- mod$summary()
      if(!is.Waiver(mod$get_data("prediction") )){
        output[["output"]][["resolution"]] <- raster::res( mod$get_data("prediction") )
        output[["output"]][["prediction"]] <- names( mod$get_data("prediction") )
      } else {
        output[["output"]][["resolution"]] <- NA
        output[["output"]][["prediction"]] <- NA
      }
      # Calculate partial estimates if set
      if(partial){
        if(verbose && getOption('ibis.setupmessages')) myLog('[Export]','green',paste0('Calculating partial variable contributions...'))
        message("Not yet added") # TODO:
        output[["output"]][["partial"]] <- NA
      } else {
        output[["output"]][["partial"]] <- NA
      }

    } else if(inherits(mod, "BiodiversityScenario")){
      # Summarize the model object
      model <- mod$get_model()

      # Model input summary in a tibble
      output[["input"]][["extent"]] <- as.matrix( raster::extent( model$model$background ) )
      output[["input"]][["predictors"]] <- mod$predictors$get_names()
      output[["input"]][["timerange"]] <- mod$get_timeperiod()
      output[["input"]][["predictor_time"]] <- mod$predictors$get_time()
      # Collect settings and parameters if existing
      if( !is.Waiver(mod$get_data("constraints")) ){
        output[["input"]][["constraints"]] <- mod$get_constraints()
      }

      # Model parameters in a tibble
      output[["params"]][["id"]] <- as.character(mod$modelid)
      output[["params"]][["runname"]] <- as.character(model$model$runname)
      output[["params"]][["algorithm"]] <- class(model)[1]
      output[["params"]][["equation"]] <- mod$get_equation()
      if( "settings" %in% names(mod) ){
        output[["params"]][["settings"]] <- model$settings$data
      }

      # Model summary in a tibble and formula
      output[["output"]][["summary"]] <- mod$summary(plot = FALSE)
      if(!is.Waiver(mod$get_scenarios() )){
        sc_dim <- stars::st_dimensions(mod$get_scenarios())
        output[["output"]][["resolution"]] <- abs( c(x = sc_dim$x$delta, y = sc_dim$y$delta) )
        output[["output"]][["prediction"]] <- names(mod$get_scenarios())
      } else {
        output[["output"]][["resolution"]] <- NA
        output[["output"]][["prediction"]] <- NA
      }
    }
    assertthat::assert_that(
      is.list(output),
      length(output)>0
    )
    # --- #
    # Write the output
    if(ext == "rds"){
      saveRDS(output, fname)
    } else if(ext == "rdata") {
      save(output, file = fname)
    } else {
      message("No compatible file format found. No summary output file ignored created!")
    }
    rm(output)
    invisible()
  }
)
