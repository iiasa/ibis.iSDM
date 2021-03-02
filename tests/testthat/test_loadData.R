context('Check if testing data can be loaded.')

# Loading data that comes with the package
test_that('Loading data',{
  # Background Raster
  background <- raster::raster(system.file('extdata/europegrid_50km.tif', package='ibis'))
  expect_s4_class(background,'Raster')

  # Get test species
  virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis'),'points',quiet = TRUE)
  virtual_range <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis'),'range',quiet = TRUE)
  expect_s3_class(virtual_points,'sf')
  expect_s3_class(virtual_range,'sf')
  expect_true(unique(sf::st_geometry_type(virtual_points)) == 'POINT')
  expect_true(unique(sf::st_geometry_type(virtual_range)) == 'POLYGON')

  # Get list of test predictors
  ll <- list.files(system.file('extdata/predictors/',package = 'ibis'),full.names = T)
  expect_gt(length(ll),0)
  expect_true(all( assertthat::has_extension(ll,'tif') ))

  # Load them as rasters
  predictors <- raster::stack(ll);names(predictors) <- tools::file_path_sans_ext(basename(ll))
  expect_s4_class(predictors,'Raster')
  expect_s4_class(predictors,'RasterStack')
  expect_true(nlayers(predictors)>1)
  expect_true(is_comparable_raster(background,predictors))
})
