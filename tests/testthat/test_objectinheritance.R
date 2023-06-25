# Most of the testing code is obsolete. These tests were created to assess that
# model objects are self-contained
test_that('Check that distribution objects are properly inherited', {
  skip_if_not_installed('igraph')
  skip_if_not_installed('abind')

  skip_if_not_installed("cmdstanr")
  skip_if(condition = tryCatch(expr = cmdstanr::cmdstan_path(), error = function(e) return(TRUE)),
          message = "No cmdstan path")

  # Load packages
  suppressWarnings( requireNamespace("terra", quietly = TRUE) )
  suppressWarnings( requireNamespace("sf", quietly = TRUE) )

  options("ibis.setupmessages" = FALSE)

  # Get background
  background <- terra::rast(system.file('extdata/europegrid_50km.tif', package='ibis.iSDM',mustWork = TRUE))

  # Get test species
  virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'points',quiet = TRUE)
  virtual_range <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'range',quiet = TRUE)
  ll <- list.files(system.file('extdata/predictors/',package = 'ibis.iSDM',mustWork = TRUE),full.names = T)

  predictors <- terra::rast(ll);names(predictors) <- tools::file_path_sans_ext(basename(ll))

  # Define distribution object
  x <- distribution(background)

  # Biodiversity
  expect_equal(x$biodiversity$length(),0)
  x |> add_biodiversity_poipo(virtual_points, field_occurrence = 'Observed', name = 'Virtual points')
  expect_equal(x$biodiversity$length(),0)
  # Multiple
  x |> add_biodiversity_poipo(virtual_points, field_occurrence = 'Observed', name = 'Virtual points') |>
    add_biodiversity_polpo(virtual_range, field_occurrence = 'Observed', name = 'Virtual points')
  expect_equal(x$biodiversity$length(),0)
  x |> add_biodiversity_poipo(virtual_points, field_occurrence = 'Observed', name = 'Virtual points') |>
    add_biodiversity_polpo(virtual_range, field_occurrence = 'Observed', name = 'Virtual points')
  expect_equal(x$biodiversity$length(),0)

  # Offsets
  suppressMessages( suppressWarnings( x |> add_offset_range(virtual_range) ) )
  expect_s3_class(x$offset, "Waiver")

  # -- #
  # Call predictors
  add_predictors(x, predictors, transform = 'none',derivates = 'none', priors = NULL)
  expect_true(is.Waiver(x$predictors))

  y <- x |> add_predictors(predictors)
  expect_length(x$get_predictor_names(), 0)
  expect_length(y$get_predictor_names(), 14)

  # Add elevation
  y <- x |> add_predictor_elevationpref(predictors$elevation_mean_50km, 500, 1000)
  expect_length(y$get_predictor_names(), 2)
  y <- x |> add_predictors(predictors) |>
    add_predictor_elevationpref(predictors$elevation_mean_50km, 500, 1000)
  expect_length(y$get_predictor_names(), 16)

  # Add range
  y <- x |> add_predictors(predictors) |>
    add_predictor_range(virtual_range, method = "binary")
  expect_length(y$get_predictor_names(), 15)
  z <- y |> add_predictor_range(virtual_range, method = "distance")
  expect_length(y$get_predictor_names(), 16)

  # -- #
  # Latent effect check
  x |> add_latent_spatial(method = "spde",priors = NULL)
  expect_true(is.Waiver(x$latentfactors))

  # Engine
  x |> engine_gdb(boosting_iterations = 500)
  expect_true(is.Waiver(x$engine))
  x |> engine_stan()
  expect_true(is.Waiver(x$engine))

  # Priors
  x |> add_predictors(predictors, transform = 'none',derivates = 'none',priors = priors(INLAPrior(names(predictors)[1],'normal')))
  expect_true(is.Waiver(x$priors))
  x |> add_latent_spatial(method = "spde", priors = priors(INLAPrior('spde','prior.range')))
  expect_true(is.Waiver(x$priors))
  # Two different priors
  x |>
    add_predictors(predictors, transform = 'none',derivates = 'none',priors = priors(INLAPrior(names(predictors)[1],'normal'))) |>
    add_latent_spatial(method = "spde", priors = priors(INLAPrior('spde','prior.range')))
  expect_true(is.Waiver(x$priors))

  # Check variable removal
  xx <- x |> add_predictors(predictors)
  xx |> rm_predictors(names = "hmi_mean_50km")
  expect_length(xx$get_predictor_names(), 14)

  # --- #
  # Create a settings object
  # Define settings object for any other information
  settings <- bdproto(NULL, Settings)
  settings$set('test', 1)
  expect_equal(settings$length(), 1)
  settings2 <- settings
  settings2 <- settings2$set('test2', 1, copy =TRUE) # This returns a settings object
  expect_equal(settings$length(), 1)

  # --- #
  # Add multiple offsets
  xx <- x |> add_predictors(predictors) |>
    add_offset_range(layer = virtual_range)
  expect_length(xx$get_offset(), 1)
  xx <- xx |> add_offset_bias(layer = predictors$hmi_mean_50km)
  expect_length(xx$get_offset(), 2)
})
