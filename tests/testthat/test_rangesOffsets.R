# First check that offsets work
test_that('Load ranges and add them to distribution object', {

  # Igraph should be by default installed, but check
  skip_if_not_installed('igraph')
  suppressWarnings( requireNamespace("igraph", quietly = TRUE) )
  suppressWarnings( requireNamespace("terra", quietly = TRUE) )

  # Set to verbose
  options("ibis.setupmessages" = FALSE)

  # --- #
  # Load data
  # Background Raster
  background <- terra::rast(system.file('extdata/europegrid_50km.tif', package='ibis.iSDM',mustWork = TRUE))
  # Get test species
  virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'points',quiet = TRUE)
  virtual_range <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'range',quiet = TRUE)
  # Get list of test predictors
  ll <- list.files(system.file('extdata/predictors/', package = 'ibis.iSDM', mustWork = TRUE),full.names = TRUE)
  # Load them as rasters
  predictors <- terra::rast(ll);names(predictors) <- tools::file_path_sans_ext(basename(ll))

  # Now set them one up step by step
  x <- distribution(background)
  # This will raise a warning since projection is different
  suppressMessages(
        x |> add_predictor_range(virtual_range, method = 'distance')
    )

  # Try and add a range as raster
  virtual_range_ras <- terra::rasterize(virtual_range, background)
  expect_type(virtual_range_ras, 'S4')
  expect_s4_class(virtual_range_ras, 'SpatRaster')

  # Add the rasterized range
  suppressWarnings( y <- x |> add_predictor_range(virtual_range_ras) )
  expect_vector(y$get_predictor_names(),'precomputed_range')

  # Artificially aggregate the range
  virtual_range_ras <- terra::aggregate(virtual_range_ras, 5)
  suppressWarnings( expect_s3_class( x |> add_predictor_range(virtual_range_ras), class = "BiodiversityDistribution" ) )

  # --------- #
  # Add offsets
  suppressMessages( y <- x |> add_offset_bias(layer = predictors$hmi_mean_50km) )
  expect_equal(y$get_offset(), "hmi_mean_50km")

  suppressMessages( y <- x |> add_offset(layer = virtual_range) )
  expect_s4_class(y$offset, "SpatRaster")

  suppressWarnings(
    y <- x |> add_offset_elevation(elev = predictors$elevation_mean_50km,pref = c(100,800))
  )
  expect_s4_class(y$offset, "SpatRaster")
  # --------- #

  # --------- #
  # Build a full model with various elements
  suppressWarnings(
    x <- distribution(background) |>
      add_predictors(predictors) |>
      add_biodiversity_poipo(virtual_points,field_occurrence = "Observed") |>
      add_offset_range(virtual_range,distance_function = "linear",distance_max = 300,
                       distance_clip = TRUE) |>
      add_offset_elevation(elev = predictors$elevation_mean_50km,pref = c(100,800)) |>
      add_offset_bias(layer = predictors$hmi_mean_50km) |>
      engine_glm()
  )
  expect_length(x$get_offset(), 3)

  # Train
  suppressWarnings(
    fit <- train(x,only_linear = T)
  )
  expect_s4_class(fit$get_data(), "SpatRaster")
  expect_true(fit$has_offset())
  # --------- #

})
