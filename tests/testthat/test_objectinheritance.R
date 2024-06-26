# Most of the testing code is obsolete. These tests were created to assess that
# model objects are self-contained
test_that('Check that distribution objects are properly inherited', {
  skip_if_not_installed('igraph')
  skip_if_not_installed('abind')
  skip_if_not_installed('INLA')

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

  # For Poipa
  pa <- virtual_points |> add_pseudoabsence(field_occurrence = "Observed",template = background)
  x |> add_biodiversity_poipa(pa, field_occurrence = "Observed",docheck = FALSE)
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

  # Remove a predictor
  y |> rm_predictors("ndvi_mean_50km")
  expect_length(y$get_predictor_names(), 14)
  y <- y |> rm_predictors("ndvi_mean_50km")
  expect_length(y$get_predictor_names(), 13)
  expect_error(y |> rm_predictors("ndvi_mean_50km")) # Trying to remove it again should lead to an error

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
  expect_length(y$get_predictor_names(), 15)
  expect_length(z$get_predictor_names(), 16)

  # -- #
  # Latent effect check
  x |> add_latent_spatial(method = "spde",priors = NULL)
  expect_true(is.Waiver(x$latentfactors))

  # Engine
  x |> engine_glm()
  expect_true(is.Waiver(x$engine))

  # Priors
  x |> add_predictors(predictors, transform = 'none',derivates = 'none',
                      priors = priors(GDBPrior(names(predictors)[1],'increasing')))
  expect_true(is.Waiver(x$priors))
  x |> add_latent_spatial(method = "spde",
                          priors = priors(INLAPrior('spde','prior.range')))
  expect_true(is.Waiver(x$priors))
  # Two different priors
  x |>
    add_predictors(predictors, transform = 'none',derivates = 'none',
                   priors = priors(INLAPrior(names(predictors)[1],'normal'))) |>
    add_latent_spatial(method = "spde",
                       priors = priors(INLAPrior('spde','prior.range')))
  expect_true(is.Waiver(x$priors))

  # Check variable removal
  xx <- x |> add_predictors(predictors)
  expect_length(x$get_predictor_names(), 0)
  xx |> rm_predictors(names = "hmi_mean_50km")
  expect_length(xx$get_predictor_names(), 14)

  # --- #
  # Create a settings object
  # Define settings object for any other information
  settings <- Settings$new()
  settings$set('test', 1)
  expect_equal(settings$length(), 1)
  settings2 <- settings$clone()
  settings2 <- settings2$set('test2', 1) # This returns a settings object
  expect_equal(settings$length(), 1)

  # --- #
  # Add multiple offsets
  xx <- x |> add_predictors(predictors) |>
    add_offset_range(layer = virtual_range)
  expect_length(xx$get_offset(), 1)
  xx <- xx |> add_offset_bias(layer = predictors$hmi_mean_50km)
  expect_length(xx$get_offset(), 2)
})

# Modify predictors objects
test_that('Modify predictors objects', {

  # Load packages
  suppressWarnings( requireNamespace("terra", quietly = TRUE) )
  suppressWarnings( requireNamespace("sf", quietly = TRUE) )

  options("ibis.setupmessages" = FALSE)

  # Get background
  background <- terra::rast(system.file('extdata/europegrid_50km.tif', package='ibis.iSDM',mustWork = TRUE))

  # Get test species
  virtual_range <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'range',quiet = TRUE)
  ll <- list.files(system.file('extdata/predictors/',package = 'ibis.iSDM',mustWork = TRUE),full.names = T)

  predictors <- terra::rast(ll);names(predictors) <- tools::file_path_sans_ext(basename(ll))

  # Define distribution object
  expect_no_error(
    x <- distribution(background) |> add_predictors(predictors)
  )
  expect_s3_class(x$predictors, "PredictorDataset")

  # Clone and alter
  pp <- x$predictors$clone()
  pp$rm_data('bio01_mean_50km')
  expect_equal(pp$length(),13)
  # Check for original
  expect_equal(x$predictors$length(),14)

  # Crop
  expect_no_warning(
    pp$crop_data(virtual_range)
  )
  expect_gt(x$predictors$ncell(), pp$ncell())

  # Mask
  expect_no_warning(
    pp$mask(virtual_range)
  )
  expect_gt(x$predictors$ncell(), pp$ncell())
  pp <- x$predictors$clone()
  expect_equal(x$predictors$ncell(), pp$ncell())

  # Manual reset
  y <- x$clone(deep = TRUE)
  y$predictors$crop_data(virtual_range)
  expect_gt(x$predictors$ncell(), y$predictors$ncell())

})

# Modify scenario objects
test_that('Modify scenario objects', {

  # Load packages
  suppressWarnings( requireNamespace("terra", quietly = TRUE) )
  suppressWarnings( requireNamespace("sf", quietly = TRUE) )
  suppressWarnings( requireNamespace("stars", quietly = TRUE) )

  options("ibis.setupmessages" = FALSE)

  # Get background
  background <- terra::rast(system.file('extdata/europegrid_50km.tif', package='ibis.iSDM',mustWork = TRUE))

  # Get test species
  virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'points',quiet = TRUE)
  virtual_range <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'range',quiet = TRUE)
  ll <- list.files(system.file('extdata/predictors/',package = 'ibis.iSDM',mustWork = TRUE),full.names = T)

  predictors <- terra::rast(ll);names(predictors) <- tools::file_path_sans_ext(basename(ll))

  # Make a dummy model
  expect_no_error(
    fit <- distribution(background) |> add_biodiversity_poipo(virtual_points,field_occurrence = "Observed") |> add_predictors(predictors) |> engine_glm() |> train()
  )

  # Define scenario object
  expect_no_error(
    sc <- scenario(fit, copy_model = TRUE)
  )
  expect_equal(sc$modelid, fit$id)

  # Add some predictors
  future_dummy <- predictors
  terra::time(future_dummy) <- rep("2020-01-01", terra::nlyr(future_dummy)) |> as.Date()
  expect_no_error(
    sc |> add_predictors(future_dummy)
  )
  expect_s3_class(sc$get_predictors(), 'Waiver')
  sc <- sc |> add_predictors(future_dummy)

  # Add constraint
  sc |> add_constraint(method = 'sdd_nexp', value = 1e3)
  expect_s3_class(sc$get_constraints(), 'Waiver')
  # Try and add the same constraint. By default this should replace the previous one
  sc <- sc |> add_constraint(method = 'sdd_nexp', value = 1e3)
  expect_type(sc$get_constraints(), "list")
  expect_length(sc$get_constraints(),1)
  sc2 <- sc |> add_constraint(method = 'sdd_nexp', value = 1e3)
  expect_type(sc2$get_constraints(), "list")
  expect_length(sc2$get_constraints(),1)
  sc2 <- sc |> add_constraint_adaptability(method = "nichelimit")
  expect_type(sc2$get_constraints(), "list")
  expect_length(sc2$get_constraints(),2)

})
