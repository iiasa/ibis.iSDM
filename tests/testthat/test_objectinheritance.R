test_that('Check that objects are properly inherited', {
  # Load packages
  require(raster)
  require(sf)

  # Get background
  background <- raster::raster(system.file('extdata/europegrid_50km.tif', package='ibis.iSDM'))

  # Get test species
  virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM'),'points',quiet = TRUE)
  virtual_range <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM'),'range',quiet = TRUE)
  ll <- list.files('inst/extdata/predictors/',full.names = T)
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
    add_biodiversity_polpo(virtual_range, field_occurrence = 'Observed', name = 'Virtual points') %>%
    add_biodiversity_polpo(virtual_range, field_occurrence = 'Observed', name = 'Virtual points',simulate = TRUE)
  expect_equal(x$biodiversity$length(),0)

  # Call predictors
  add_predictors(x, predictors, transform = 'none',derivates = 'none',priors = NULL)
  expect_true(is.Waiver(x$predictors))

  # Latent effect check
  x %>% add_latent_spatial(priors = NULL)
  expect_true(is.Waiver(x$latentfactors))

  # Engine
  x %>% engine_gdb(boosting_iterations = 500)
  expect_true(is.Waiver(x$engine))
  x %>% engine_inla()
  expect_true(is.Waiver(x$engine))

  # Priors
  x %>% add_predictors(predictors, transform = 'none',derivates = 'none',priors = priors(INLAPrior(names(predictors)[1],'abc')))
  expect_true(is.Waiver(x$priors))
  x %>% add_latent_spatial(priors = priors(INLAPrior('spde','prior.range')))
  expect_true(is.Waiver(x$priors))
  # Two different priors
  x %>%
    add_predictors(predictors, transform = 'none',derivates = 'none',priors = priors(INLAPrior(names(predictors)[1],'abc'))) %>%
    add_latent_spatial(priors = priors(INLAPrior('spde','prior.range')))
  expect_true(is.Waiver(x$priors))

  # Check variable removal
  xx <- x %>% add_predictors(predictors)
  xx %>% rm_predictors("hmi_mean_50km")
  expect_length(xx$get_predictor_names(), 14)

})
