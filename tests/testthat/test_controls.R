# ---- #
# Train a full distribution model with glm base model
test_that('Test controls and limits', {

  skip_on_cran()

  # No messages
  options(ibis.setupmessages = FALSE)

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
    add_predictors(predictors, transform = 'none',derivates = 'none')

  # Expect error
  expect_error(x |> add_limits_extrapolation(method = "nomethod"))

  # Add the various controls to see that they are working
  zones <- terra::as.factor( predictors$koeppen_50km )
  y <- x |> add_limits_extrapolation(layer = zones, method = "zones")
  expect_false(y$get_limits() |> is.Waiver())
  expect_length(y$get_limits(), 4) # Those are saved in limits in the object
  expect_true(x$get_limits() |> is.Waiver()) # Check that original object remains correct
  z <- y |> rm_limits()
  expect_true(z$get_limits() |> is.Waiver()) # Check that altered object remains correct
  expect_false(y$get_limits() |> is.Waiver()) # Check that original object remains correct

  # Add mcp limits
  y <- x |> add_limits_extrapolation(method = "mcp")
  expect_false(y$get_limits() |> is.Waiver())

  # Add NT2
  y <- x |> add_limits_extrapolation(method = "nt2")
  expect_false(y$get_limits() |> is.Waiver())

  # Add MESS
  y <- x |> add_limits_extrapolation(method = "mess")
  expect_false(y$get_limits() |> is.Waiver())

  # Train the model with limits set
  x <-  x |> add_limits_extrapolation(layer = zones, method = "zones")
  suppressWarnings(
    mod <- train(x |> engine_glm(), "test", inference_only = FALSE, only_linear = TRUE,
                 varsel = "none", verbose = FALSE)
  )

  # Run a check (should work without errors at least)
  expect_no_error( suppressMessages( check(mod) ) )

  # Expect true
  expect_true( mod$has_limits() )

  # --- #
  # Also try mess
  x <- x$rm_limits()
  x <-  x |> add_limits_extrapolation(method = "mess")
  suppressWarnings(
    mod <- train(x |> engine_glm(), "test", inference_only = FALSE, only_linear = TRUE,
                 varsel = "none", verbose = FALSE)
  )

  # Expect true
  expect_true( mod$has_limits() )

  # Create a scenario object and reuse limits
  expect_no_error( scenario(mod, reuse_limits = TRUE) )

  # --- #
  # Bias control checks
  x <- x$rm_limits()
  x <- x |> add_control_bias(layer = predictors$hmi_mean_50km)
  # Train with bias
  expect_no_error(
    suppressWarnings(
      mod <- train(x |> engine_glm(), "test", inference_only = FALSE, only_linear = TRUE,
                   varsel = "none", verbose = FALSE)
    )
  )
  expect_length( x$get_control(), 4 )
  settings <- mod$settings
  expect_equal(settings$get("bias_variable"), "hmi_mean_50km")

  # Remove control
  y <- x |> rm_control()
  expect_length(x$get_control(), 4)
  expect_length(y$get_control(), 0)
})

