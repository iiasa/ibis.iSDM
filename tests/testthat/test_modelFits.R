# Train a full distribution model with INLA
test_that('Add further tests for model fits', {

  skip_if_not_installed('glmnet')
  skip_if_not_installed('pdp')
  skip_on_travis()
  skip_on_cran()

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

  # Add pseudo absence
  abs <- pseudoabs_settings(nrpoints = 0,min_ratio = 1,method = "mcp")
  suppressMessages(
    virtual_points <- add_pseudoabsence(virtual_points,template = background, field_occurrence = "Observed", settings = abs)
  )

  # Create testing and training data
  ind <- sample(1:nrow(virtual_points), 70)
  train_data <- virtual_points[-ind,]
  test_data <- virtual_points[ind,]

  # Now set them one up step by step
  x <- distribution(background) %>%
    add_biodiversity_poipa(train_data, field_occurrence = 'Observed', name = 'Virtual points') %>%
    add_predictors(predictors, transform = 'none',derivates = 'none') %>%
    engine_glmnet()

  # Train the model
  suppressWarnings(
    mod <- train(x, "test", inference_only = FALSE, only_linear = TRUE, varsel = "none", verbose = FALSE)
  )
  expect_s4_class(mod$get_data(), "RasterLayer")

  # Threshold with independent data
  mod <- threshold(mod,method = "perc",format = "bin")
  expect_gt(mod$get_thresholdvalue(),0)
  expect_length(mod$show_rasters(), 2)

  # Summarize model
  expect_s3_class( summary(mod), "data.frame" )
  expect_s3_class( coef(mod), "data.frame" )

  # Validate
  val <- validate(mod, method = "cont")
  expect_s3_class(val, "data.frame")
  # Validate discrete
  val <- validate(mod, method = "disc")
  expect_s3_class(val, "data.frame")

  # Validate with withold data
  val <- validate(mod, method = "disc", point = test_data,point_column = "Observed")
  expect_s3_class(val, "data.frame")

  # ----------- #
  # Partial stuff
  pp <- partial(mod,x.var = "bio19_mean_50km",plot = FALSE)
  expect_s3_class(pp, "data.frame")

  # Spartial
  pp <- spartial(mod,x.var = "bio19_mean_50km",plot = FALSE)
  expect_s4_class(pp, "RasterLayer")


  # ----------- #
  # Write model outputs
  # expect_snapshot_file(write_summary(mod, "test.rds"), "test.rds")

})
