# Loading data that comes with the package
test_that('Check that data can be loaded.',{

  # Background Raster
  background <- terra::rast(system.file('extdata/europegrid_50km.tif', package='ibis.iSDM',mustWork = TRUE))
  expect_s4_class(background,'SpatRaster')

  # Get test species
  virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'points',quiet = TRUE)
  virtual_range <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'range',quiet = TRUE)
  expect_s3_class(virtual_points,'sf')
  expect_s3_class(virtual_range,'sf')
  expect_true(unique(sf::st_geometry_type(virtual_points)) == 'POINT')
  expect_true(unique(sf::st_geometry_type(virtual_range)) == 'POLYGON')

  # Get list of test predictors
  ll <- list.files(system.file('extdata/predictors/',package = 'ibis.iSDM',mustWork = TRUE),full.names = T)
  expect_gt(length(ll),0)
  expect_true(all( assertthat::has_extension(ll,'tif') ))

  # Load them as rasters
  predictors <- terra::rast(ll);names(predictors) <- tools::file_path_sans_ext(basename(ll))
  expect_s4_class(predictors,'SpatRaster')
  expect_true(terra::nlyr(predictors)>1)
  expect_true(is_comparable_raster(background, predictors))
})

test_that('Check that test scenarios can be loaded.',{
  # Load the scenario data
  skip_if_not_installed('abind')
  skip_on_os(os = "mac") # Added since stars throws errors here?

  requireNamespace("stars")
  requireNamespace("ncmeta")
  requireNamespace("abind")

  ll <- list.files(system.file('extdata/predictors_presfuture/',package = 'ibis.iSDM',mustWork = TRUE),full.names = TRUE)
  expect_vector(ll)
  expect_length(ll,9)
  expect_true(all( assertthat::has_extension(ll,'nc') ))
  expect_true(all( file.exists(ll) ))

  # Load stars files
  suppressWarnings( sc <- stars::read_stars(ll) )
  # Still having warnings for the bioclimatic files
  expect_type(sc[1], 'list')
  # Average
  o <- mean(sc[[1]],na.rm = TRUE)
  expect_true(o > 10)

  # Get time attributes
  tt <- stars::st_get_dimension_values(sc, 'Time')
  expect_length(tt, 86)
})
