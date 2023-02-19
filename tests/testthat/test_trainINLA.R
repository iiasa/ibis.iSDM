# First check that INLA works
test_that('Check that INLA works', {
  skip_on_travis()
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed('INLA')

  suppressWarnings(
    suppressPackageStartupMessages( library(INLA) )
  )
  options("ibis.setupmessages" = FALSE)

  # Use test data that comes with INLA
  data(Epil)
  observed <- Epil[1:30, 'y']

  Epil <- rbind(Epil, Epil[1:30, ])
  Epil[1:30, 'y'] <- NA

  # Set up formula and train
  formula = y ~ Trt + Age + V4 + f(Ind, model="iid") + f(rand,model="iid")
  result = inla(formula, family="poisson", data = Epil, control.predictor = list(compute = TRUE, link = 1))

  expect_type(result,'list')
  expect_null(result$waic)
  expect_true(result$ok)
  expect_equal(nrow(result$summary.random$rand), 236)

})

# Train a full distribution model with INLA
test_that('Train a distribution model with INLA', {

  skip_if_not_installed('INLA')
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

  # Now set them one up step by step
  x <- distribution(background) %>%
    add_biodiversity_poipo(virtual_points, field_occurrence = 'Observed', name = 'Virtual points') %>%
    add_predictors(predictors[[c('slope_mean_50km','bio01_mean_50km','CLC3_132_mean_50km')]], transform = 'none',derivates = 'none') %>%
    engine_inla(
      max.edge = c(.5, 3),
      offset = c(0.5, 1),
      cutoff = 0.5,
      proj_stepsize = 1
    )
  # Train the model
  suppressWarnings(
    mod <- train(x, "test", inference_only = TRUE,only_linear = TRUE, varsel = "none", verbose = FALSE)
  )

  # Expect summary
  expect_s3_class(summary(mod), "data.frame")
  expect_s3_class(mod$show_duration(), "difftime")
  expect_equal(length(mod$show_rasters()), 0) # Now predictions found
  # Fit with predictions
  suppressWarnings(
    mod <- train(x, "test", inference_only = FALSE,only_linear = TRUE, varsel = "none", verbose = FALSE)
  )
  expect_equal(length(mod$show_rasters()), 1) # Now predictions found

})
