# Ensure that future setup works
# This might be tricky in CI as different setups make use of different parallel
# configurations
test_that('Testing parallel setup', {

  # Set to verbose
  options("ibis.setupmessages" = FALSE)

  skip_if_not_installed("future")
  skip_if_not_installed("doFuture")

  # Load data
  # Background Raster
  background <- terra::rast(system.file('extdata/europegrid_50km.tif', package='ibis.iSDM',mustWork = TRUE))
  # Get test species
  virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'points',quiet = TRUE)
  # Get list of test predictors
  ll <- list.files(system.file('extdata/predictors/',package = 'ibis.iSDM',mustWork = TRUE),full.names = T)
  # Load them as rasters
  predictors <- terra::rast(ll);names(predictors) <- tools::file_path_sans_ext(basename(ll))

  # Add pseudo absence
  abs <- pseudoabs_settings(nrpoints = 0,min_ratio = 1,method = "mcp")
  suppressMessages(
    poa <- add_pseudoabsence(virtual_points,template = background, field_occurrence = "Observed", settings = abs)
  )

  # Now set them one up step by step
  x <- distribution(background) |>
    add_biodiversity_poipa(poipa = poa,field_occurrence = 'Observed',docheck = FALSE) |>
    add_predictors(predictors, transform = 'none',derivates = 'none') |>
    engine_glm()

  expect_no_error(
    fit1 <- train(x, "test", inference_only = FALSE)
  )

  # Now enable parallel
  expect_no_error(
    expect_invisible(
      ibis_enable_parallel()
    )
  )
  # Set nr of threads
  expect_no_error(
    expect_invisible(
      ibis_set_threads(2)
    )
  )

  # Set strategy
  expect_no_error(
    expect_invisible(
      ibis_set_strategy(strategy = "sequential")
    )
  )

  # --- #
  # Now define a plan
  ibis_future()

  expect_no_error(
    fit2 <- train(x, "test", inference_only = FALSE)
  )

  # Try with multi-session
  ibis_future(strategy = "multisession")

  expect_no_error(
    fit3 <- train(x, "test", inference_only = FALSE)
  )

  # Assume they are all identical
  expect_gte(
    cor(fit1$get_coefficients()[,2], fit2$get_coefficients()[,2]),
    0.99
  )
  expect_gte(
    cor(fit1$get_coefficients()[,2], fit3$get_coefficients()[,2]),
    0.99
  )

  # Set parallel to FALSE again
  options('ibis.runparallel' = FALSE)
})
