
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

# Setting up a distribution model
test_that('Setting up a distribution model',{
  # Background Raster
  background <- raster::raster(system.file('extdata/europegrid_50km.tif', package='ibis'))
  # Get test species
  virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis'),'points',quiet = TRUE)
  virtual_range <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis'),'range',quiet = TRUE)
  # Get list of test predictors
  ll <- list.files(system.file('extdata/predictors/',package = 'ibis'),full.names = T)
  # Load them as rasters
  predictors <- raster::stack(ll);names(predictors) <- tools::file_path_sans_ext(basename(ll))

  # Now set them one up step by step
  x <- distribution(background)
  expect_s3_class(x,'BiodiversityDistribution')
  expect_s4_class(x$background,'Raster')
  expect_error(x$biodiversity$get_data())
  expect_equal(x$biodiversity$length(),0)
  expect_type(x$show_background_info(),'list')

  # Now add one variable
  x <- x %>% add_biodiversity_poipo(virtual_points,field_occurrence = 'Observed',name = 'Virtual points')
  expect_message(x$biodiversity,NA)
  expect_equal(x$biodiversity$length(),1)
  expect_equal(x$biodiversity$get_equations()[[1]],'<Default>')
  # And another
  x <- x %>% add_biodiversity_polpo(virtual_range,field_occurrence = 'Observed',name = 'Virtual range')
  expect_equal(x$biodiversity$length(),2)
  expect_equal(sum(x$biodiversity$get_observations()),210)

  # Add Predictors
  x <- x %>% add_predictors(predictors)
  expect_equal(x$predictors$length(),14)

  x <- x %>% engine_inla()
  expect_output(x$predictor_names(),'vector')
  expect_s3_class(x$engine$data$mesh,'inla.mesh')

})
