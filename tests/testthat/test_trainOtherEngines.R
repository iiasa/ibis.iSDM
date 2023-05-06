# Train a full distribution model with INLA
test_that('Train a distribution model with XGboost', {

  skip_if_not_installed('xgboost')
  skip_on_travis()
  skip_on_cran()

  suppressWarnings( requireNamespace('xgboost', quietly = TRUE) )

  # Load data
  # Background Raster
  background <- terra::rast(system.file('extdata/europegrid_50km.tif', package='ibis.iSDM',mustWork = TRUE))
  # Get test species
  virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'points',quiet = TRUE)
  virtual_range <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'range',quiet = TRUE)
  # Get list of test predictors
  ll <- list.files(system.file('extdata/predictors/',package = 'ibis.iSDM',mustWork = TRUE),full.names = T)
  # Load them as rasters
  predictors <- terra::rast(ll);names(predictors) <- tools::file_path_sans_ext(basename(ll))

  # Now set them one up step by step
  x <- distribution(background) |>
    add_biodiversity_poipo(virtual_points, field_occurrence = 'Observed', name = 'Virtual points') |>
    add_predictors(predictors, transform = 'none',derivates = 'none') |>
    engine_xgboost()

  # Train the model
  suppressWarnings(
    mod <- train(x, "test", inference_only = FALSE, only_linear = TRUE,
                 varsel = "none", verbose = FALSE)
  )

  # Expect summary
  expect_s3_class(summary(mod), "data.frame")
  expect_s3_class(mod$show_duration(), "difftime")
  expect_equal(length(mod$show_rasters()), 1) # Now predictions found

  # --- #
  # Some checks
  expect_true("get_data" %in% names(mod))
  expect_true("plot" %in% names(mod))
  expect_true("summary" %in% names(mod))

  # Test some basic non-sense calculations
  tr <- threshold(mod)
  expect_type(tr$get_thresholdvalue(), "double")
  ex <- ensemble(mod, mod)
  expect_s4_class(ex, "SpatRaster")

})

# ---- #
# Train a full distribution model with INLA
test_that('Train a distribution model with Breg', {

  skip_if_not_installed('BoomSpikeSlab')
  skip_on_travis()
  skip_on_cran()

  suppressWarnings( requireNamespace('BoomSpikeSlab', quietly = TRUE) )

  # Load data
  # Background Raster
  background <- terra::rast(system.file('extdata/europegrid_50km.tif', package='ibis.iSDM',mustWork = TRUE))
  # Get test species
  virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'points',quiet = TRUE)
  virtual_range <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'range',quiet = TRUE)
  # Get list of test predictors
  ll <- list.files(system.file('extdata/predictors/',package = 'ibis.iSDM',mustWork = TRUE),full.names = T)
  # Load them as rasters
  predictors <- terra::rast(ll);names(predictors) <- tools::file_path_sans_ext(basename(ll))

  # Now set them one up step by step
  x <- distribution(background) |>
    add_biodiversity_poipo(virtual_points, field_occurrence = 'Observed', name = 'Virtual points') |>
    add_predictors(predictors, transform = 'none',derivates = 'none') |>
    engine_breg(iter = 100)

  # Train the model
  suppressWarnings(
    mod <- train(x, "test", inference_only = FALSE, only_linear = TRUE,
                 varsel = "none", verbose = FALSE)
  )

  # Expect summary
  expect_s3_class(summary(mod), "data.frame")
  expect_s3_class(mod$show_duration(), "difftime")
  expect_equal(length(mod$show_rasters()), 1) # Now predictions found

  # --- #
  # Some checks
  expect_true("get_data" %in% names(mod))
  expect_true("plot" %in% names(mod))
  expect_true("summary" %in% names(mod))

  # Test some basic non-sense calculations
  tr <- threshold(mod)
  expect_type(tr$get_thresholdvalue(), "double")
  ex <- ensemble(mod, mod)
  expect_s4_class(ex, "SpatRaster")

})
