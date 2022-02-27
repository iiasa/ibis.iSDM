#' Generic function to write outputs
#'
#' @description
#' The `write_output` function is a generic wrapper to writing any output files created with
#' the [`ibis.iSDM-package`]. It is possible to write outputs of fitted [`DistributionModel`],
#' [`BiodiversityScenario`] or individual [`Raster`] or [`stars`] objects. In case a [`data.frame`]
#' is supplied, the output is written as csv file.
#' @note
#' By default output files will be overwritten if already existing!
#' @param mod Provided [`DistributionModel`], [`BiodiversityScenario`], [`Raster`] or [`stars`] object.
#' @param fname A [`character`] depicting an output filename.
#' @param dt A [`character`] for the output datatype. Following the [`raster::dataType()`] options (Default: \code{'FLT4S'}).
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
                    function(mod, fname, dt = "FLT4S", ...) standardGeneric("write_output"))

#' @name write_output
#' @rdname write_output
#' @usage \S4method{write_output}{ANY, character, character}(mod, fname, dt)
methods::setMethod(
  "write_output",
  methods::signature("ANY"),
  function(mod, fname, dt = "FLT4S", ...){
    assertthat::assert_that(
      !missing(mod),
      is.character(fname),
      is.character(dt)
    )
  # This function will only capture the distribution model object and will save them separately
  if(any(class(mod) %in% getOption("ibis.engines")) ){
    # FIXME: If errors occur, check harmonization of saving among engines.
    if(getOption('ibis.setupmessages')) myLog('[Output]','green','Saving output(s)...')
    mod$save(fname = fname)
  } else if(is.Raster(mod)){
      if(getOption('ibis.setupmessages')) myLog('[Output]','green','Saving output(s)...')
      if(raster::extension(fname) %in% c('.tif', '.TIF')) {
        writeGeoTiff(file = mod, fname = fname, dt = dt)
      } else if(raster::extension(fname) %in% c('.nc', '.NC', '.ncdf', '.NCDF')){
        writeNetCDF(file = mode, fname = fname, varName = names(mod), dt = dt)
      } else {
        stop("Output type could not be determined. Currently only geoTIFF and netCDF are supported.")
      }
  } else if(is.data.frame(mod)){
    if(getOption('ibis.setupmessages')) myLog('[Output]','green','Saving output(s)...')
    write.csv(x = mod,file = fname,...)
  } else {
    # Check that a save function exists for object
    assertthat::assert_that("save" %in%names(mod),
                            msg = "No method to save the output could be found!")
    if(getOption('ibis.setupmessages')) myLog('[Output]','green','Saving output(s)...')
    # Try a generic save
    mod$save(fname = fname)
  }
  invisible()
  }
)

#' @name write_output
#' @rdname write_output
#' @usage \S4method{write_output}{BiodiversityScenario, character, character}(mod, fname, dt)
methods::setMethod(
  "write_output",
  methods::signature(mod = "BiodiversityScenario"),
  function(mod, fname, dt = "FLT4S", ...) {
    assertthat::assert_that(
      !missing(mod),
      is.character(fname),
      is.character(dt)
    )
    # Get outputs
    mod$save(fname = fname, type = tools::file_ext(fname), dt = dt)
    invisible()
  }
)

#' @name write_output
#' @rdname write_output
#' @usage \S4method{write_output}{RasterLayer, character, character}(mod, fname, dt)
methods::setMethod(
  "write_output",
  methods::signature(mod = "RasterLayer"),
  function(mod, fname, dt = "FLT4S", ...) {
    assertthat::assert_that(
      !missing(mod),
      is.Raster(mod),
      is.character(fname),
      is.character(dt)
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
#' @usage \S4method{write_output}{RasterStack, character, character}(mod, fname, dt)
methods::setMethod(
  "write_output",
  methods::signature(mod = "RasterStack"),
  function(mod, fname, dt = "FLT4S", ...) {
    assertthat::assert_that(
      !missing(mod),
      is.Raster(mod),
      is.character(fname),
      is.character(dt)
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
#' @usage \S4method{write_output}{data.frame, character, character}(mod, fname, dt)
methods::setMethod(
  "write_output",
  methods::signature(mod = "data.frame"),
  function(mod, fname, dt = "FLT4S", ...) {
    assertthat::assert_that(
      !missing(mod),
      is.data.frame(mod),
      is.character(fname)
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
