# First check that INLA works
test_that('Load ranges and add them to distribution object', {
  skip_on_travis()
  skip_on_cran()
  skip_if_not_installed('INLA')
  skip_if_not_installed('igraph')

  require(igraph)

  options("ibis.setupmessages" = FALSE)

  # --- #
  # Load data
  # Background Raster
  background <- raster::raster(system.file('extdata/europegrid_50km.tif', package='ibis.iSDM',mustWork = TRUE))
  # Get test species
  virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'points',quiet = TRUE)
  virtual_range <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'range',quiet = TRUE)
  # Get list of test predictors
  ll <- list.files(system.file('extdata/predictors/',package = 'ibis.iSDM',mustWork = TRUE),full.names = T)
  # Load them as rasters
  predictors <- raster::stack(ll);names(predictors) <- tools::file_path_sans_ext(basename(ll))

  # Now set them one up step by step
  x <- distribution(background)
  # This will raise a warning since projection is different
  suppressMessages(
    expect_warning(
      x %>% add_predictor_range(virtual_range, method = 'distance') )
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

  # Add bias variable
  y <- x |> add_control_bias(layer = predictors$hmi_mean_50km)
  expect_type(y$bias, 'list')
  expect_length(y$get_biascontrol(), 3)

})
