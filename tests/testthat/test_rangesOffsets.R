# First check that INLA works
test_that('Load ranges and add them to distribution object', {
  skip_on_travis()
  skip_on_cran()
  skip_if_not_installed('INLA')

  options("ibis.setupmessages" = FALSE)

  # --- #
  # Load data
  # Background Raster
  background <- raster::raster(system.file('extdata/europegrid_50km.tif', package='ibis.iSDM'))
  # Get test species
  virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM'),'points',quiet = TRUE)
  virtual_range <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM'),'range',quiet = TRUE)
  # Get list of test predictors
  ll <- list.files(system.file('extdata/predictors/',package = 'ibis.iSDM'),full.names = T)
  # Load them as rasters
  predictors <- raster::stack(ll);names(predictors) <- tools::file_path_sans_ext(basename(ll))

  # Now set them one up step by step
  x <- distribution(background)
  # This won't work since not aligned
  expect_warning(
    expect_error(  x %>% add_predictor_range(virtual_range,method = 'distance') )
  )

  # Try and add a range as raster
  virtual_range_ras <- raster::rasterize(virtual_range, background)
  expect_type(virtual_range_ras,'S4')
  expect_s4_class(virtual_range_ras,'Raster')

  # Add the rasterized range
  y <- x %>% add_predictor_range(virtual_range_ras)
  expect_vector(y$get_predictor_names(),'precomputed_range')

  # Artificially aggregate the range
  virtual_range_ras <- raster::aggregate(virtual_range_ras, 5)
  expect_s3_class( x %>% add_predictor_range(virtual_range_ras),class = "BiodiversityDistribution" )

})
