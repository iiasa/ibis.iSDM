# Most of the testing code is obsolete. These tests were created to assess that
# model objects are self-contained
test_that('Check that objects are properly inherited', {
  # Load packages
  require(raster)
  require(sf)
  skip_if_not_installed('igraph')
  skip_if_not_installed('abind')

  options("ibis.setupmessages" = FALSE)

  # Get background
  background <- raster::raster(system.file('extdata/europegrid_50km.tif', package='ibis.iSDM',mustWork = TRUE))

  # Get test species
  virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'points',quiet = TRUE)
  virtual_range <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'range',quiet = TRUE)
  ll <- list.files(system.file('extdata/predictors/',package = 'ibis.iSDM',mustWork = TRUE),full.names = T)

  predictors <- raster::stack(ll);names(predictors) <- tools::file_path_sans_ext(basename(ll))

  # Define distribution object
  x <- distribution(background)

  # Biodiversity
  expect_equal(x$biodiversity$length(),0)
  x %>% add_biodiversity_poipo(virtual_points, field_occurrence = 'Observed', name = 'Virtual points')
  expect_equal(x$biodiversity$length(),0)
  # Multiple
  x %>% add_biodiversity_poipo(virtual_points, field_occurrence = 'Observed', name = 'Virtual points') %>%
    add_biodiversity_polpo(virtual_range, field_occurrence = 'Observed', name = 'Virtual points')
  expect_equal(x$biodiversity$length(),0)
  x %>% add_biodiversity_poipo(virtual_points, field_occurrence = 'Observed', name = 'Virtual points') %>%
    add_biodiversity_polpo(virtual_range, field_occurrence = 'Observed', name = 'Virtual points')
  expect_equal(x$biodiversity$length(),0)

  # Offsets
  suppressWarnings( x %>% add_offset_range(virtual_range) )
  expect_s3_class(x$offset, "Waiver")

  # Call predictors
  add_predictors(x, predictors, transform = 'none',derivates = 'none',priors = NULL)
  expect_true(is.Waiver(x$predictors))

  # Latent effect check
  x %>% add_latent_spatial(method = "spde",priors = NULL)
  expect_true(is.Waiver(x$latentfactors))

  # Engine
  x %>% engine_gdb(boosting_iterations = 500)
  expect_true(is.Waiver(x$engine))
  x %>% engine_stan()
  expect_true(is.Waiver(x$engine))

  # Priors
  x %>% add_predictors(predictors, transform = 'none',derivates = 'none',priors = priors(INLAPrior(names(predictors)[1],'normal')))
  expect_true(is.Waiver(x$priors))
  x %>% add_latent_spatial(method = "spde", priors = priors(INLAPrior('spde','prior.range')))
  expect_true(is.Waiver(x$priors))
  # Two different priors
  x %>%
    add_predictors(predictors, transform = 'none',derivates = 'none',priors = priors(INLAPrior(names(predictors)[1],'normal'))) %>%
    add_latent_spatial(method = "spde", priors = priors(INLAPrior('spde','prior.range')))
  expect_true(is.Waiver(x$priors))

  # Check variable removal
  xx <- x %>% add_predictors(predictors)
  xx %>% rm_predictors("hmi_mean_50km")
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

})
