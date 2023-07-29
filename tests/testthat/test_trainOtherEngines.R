# Train a full distribution model with INLA
test_that('Train a distribution model with XGboost', {

  skip_if_not_installed('xgboost')
  skip_if_not_installed('pdp')
  skip_on_travis()
  skip_on_cran()

  suppressWarnings( requireNamespace('xgboost', quietly = TRUE) )
  suppressWarnings( requireNamespace('pdp', quietly = TRUE) )

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

  ex_sd <- ensemble(mod, mod, uncertainty = "sd")
  ex_range <- ensemble(mod, mod, uncertainty = "range")
  ex_pca <- ensemble(mod, mod, uncertainty = "pca")

  expect_named(object = ex_sd, expected = c("ensemble_mean", "sd_mean"))
  expect_named(object = ex_range, expected = c("ensemble_mean", "range_mean"))
  expect_named(object = ex_pca, expected = c("ensemble_mean", "pca_mean"))

})

# ---- #
# Train a full distribution model with Breg
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

  # Does limiting raster work?
  suppressMessages( expect_s4_class(limiting(mod, plot = FALSE), "SpatRaster") )

})

# ---- #
# Train a full distribution model with gdb
test_that('Train a distribution model with GDB', {

  skip_if_not_installed('mboost')
  skip_on_travis()
  skip_on_cran()

  suppressWarnings( requireNamespace('mboost', quietly = TRUE) )

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
    engine_gdb(iter = 100)

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

  # Does limiting raster work?
  suppressMessages( expect_s4_class(limiting(mod, plot = FALSE), "SpatRaster") )

})

# ---- #
# Train a full distribution model with glmnet
test_that('Train a distribution model with glmnet', {

  skip_if_not_installed('glmnet')
  skip_if_not_installed('pdp')
  skip_on_travis()
  skip_on_cran()

  suppressWarnings( requireNamespace('glmnet', quietly = TRUE) )
  suppressWarnings( requireNamespace('pdp', quietly = TRUE) )

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
    engine_glmnet(alpha = 1)

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

  # Added here to tests as it is quick
  expect_no_error( partial_density(mod = mod, x.var = "elevation_mean_50km", df = FALSE))
  expect_s3_class( partial_density(mod = mod, x.var = "elevation_mean_50km", df = TRUE), "data.frame")

  # Does limiting raster work?
  suppressMessages( expect_s4_class(limiting(mod, plot = FALSE), "SpatRaster") )

  # ------- #
  # Create some mcps and collect data using some of the internal functions
  suppressMessages(
    suppressWarnings(
      m <- create_mcp(biod = mod$model,limits = list(mcp_buffer = 10))
    )
  )
  expect_true(inherits(m, "sf"))
  expect_equal(nrow(m), 1)

  p <- collect_occurrencepoints(mod$model,tosf = TRUE)
  expect_true(nrow(p)>0)
})

# ---- #
# Train a full distribution model with bart
test_that('Train a distribution model with bart', {

  skip_if_not_installed('dbarts')
  skip_on_travis()
  skip_on_cran()

  suppressWarnings( requireNamespace('dbarts', quietly = TRUE) )

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

  # Add pseudo absence points
  poa <- add_pseudoabsence(virtual_points,field_occurrence = "Observed", template =  background,
                           settings = pseudoabs_settings(method = "mcp"))

  # Now set them one up step by step
  x <- distribution(background) |>
    add_biodiversity_poipa(poa, field_occurrence = 'Observed', name = 'Virtual points') |>
    add_predictors(predictors, transform = 'none',derivates = 'none') |>
    engine_bart(iter = 100)

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
  ex <- ensemble(mod, tr)
  expect_s4_class(ex, "SpatRaster")

})

# TODO: Engine stan to be tested
