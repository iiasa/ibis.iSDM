# Other generic functions in the package
test_that('Test various thresholds calculations', {

  suppressWarnings(
    requireNamespace("terra", quietly = TRUE)
  )
  # Dummy layer
  r1 <- terra::rast(nrows = 10, ncols = 10, res = 0.05, xmin = -1.5, xmax = 1.5, ymin = -1.5, ymax = 1.5, vals = rnorm(3600,mean = .5,sd = .1))

  # Make some dummy points
  ss <- terra::spatSample(r1, size = 10, as.points = TRUE) |> sf::st_as_sf()
  expect_s3_class(ss, "sf")
  ss$observed <- rbinom(10, 1, .5)
  if(length(unique(ss$observed))<2) ss$observed <- c(0,0,0,1,0,1,1,1,0,0) # Resample

  expect_error( tr <- threshold(r1, method = "fixed") )
  tr1 <- threshold(r1, method = "fixed", value = .5)
  expect_s4_class(tr1, "SpatRaster")

  # Try different formats
  tr1a <- threshold(r1, method = "fixed", value = .5,format = "binary")
  expect_s4_class(tr1a, "SpatRaster")
  # Check that attributes are there
  expect_equal(attr(tr1a,'threshold'), 0.5)
  tr1b <- threshold(r1, method = "fixed", value = .5,format = "norm")
  expect_s4_class(tr1b, "SpatRaster")
  tr1c <- threshold(r1, method = "fixed", value = .5,format = "perc")
  expect_s4_class(tr1c, "SpatRaster")
  expect_length(unique(tr1c)[,1], 10)


  expect_error( tr2 <- threshold(r1, method = "perc", value = .1) )
  tr2 <- threshold(r1, method = "perc", value = .1, point = ss)
  expect_s4_class(tr2, "SpatRaster")

  tr3 <- threshold(r1, method = "mtp", point = ss)
  expect_s4_class(tr3, "SpatRaster")

  expect_error(tr4 <- threshold(r1, method = "min.cv", point = ss))
  tr4 <- threshold(r1, method = "min.cv",value = .2, point = ss)
  expect_s4_class(tr4, "SpatRaster")

  # Rename the point information and test again
  ss <- ss |> dplyr::rename(new = observed)
  expect_no_error(
    tr5 <- threshold(r1, method = "mtp", point = ss, field_occurrence = "new")
  )

  # --- #
  # Skip below if modeva not installed
  skip_if_not_installed("modEvA")
  suppressWarnings(
    requireNamespace("modEvA", quietly = TRUE)
  )
  expect_no_error(
    tr6 <- threshold(r1, method = "TSS", point = ss, field_occurrence = "new")
  )
  expect_s4_class(tr6, 'SpatRaster')
  # --- #
})

# --- #
# Test threshold on fitted model
test_that('Test Thresholds on fitted models', {

  suppressWarnings(
    requireNamespace("terra", quietly = TRUE)
  )

  # Don't print out as many messages
  options("ibis.setupmessages" = FALSE)

  # Load files
  background <- terra::rast(system.file('extdata/europegrid_50km.tif', package='ibis.iSDM',mustWork = TRUE))
  virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'points',quiet = TRUE)
  # Get list of test predictors
  ll <- list.files(system.file('extdata/predictors',package = 'ibis.iSDM',mustWork = TRUE),full.names = T)
  predictors <- terra::rast(ll);names(predictors) <- tools::file_path_sans_ext(basename(ll))

  # --- #
  # Fit a dummy model
  fit <- distribution(background) |>
    add_biodiversity_poipo(virtual_points,field_occurrence = "Observed") |>
    add_predictors(predictors) |>
    engine_glm() |>
    train()

  # Now do thresholds on the fitted model
  tr <- threshold(fit, method = "mtp")
  expect_gt(tr$get_thresholdvalue(), 0)
  expect_length(tr$show_rasters(), 2)
  rm(tr)
  expect_no_error(tr <- threshold(fit, method = "perc") )

  # Try with an externally supplied point dataset
  expect_error( tr <- threshold(fit, method = "mtp", point = virtual_points) )
  # Renamed
  expect_no_error( tr <- threshold(fit, method = "mtp",
                                   point = virtual_points,
                                   field_occurrence = "Observed") )

  # --- #
  # Skip below if modeva not installed
  skip_if_not_installed("modEvA")
  suppressWarnings(
    requireNamespace("modEvA", quietly = TRUE)
  )
  poipa <- add_pseudoabsence(virtual_points, field_occurrence = "Observed",template = background)

  fit <- distribution(background) |>
    add_biodiversity_poipa(poipa, field_occurrence = "Observed") |>
    add_predictors(predictors) |>
    engine_glm() |>
    train()

  expect_error(
    tr2 <- threshold(fit, method = "TSS", point = poipa, field_occurrence = "new")
  )
  expect_no_error(
    tr2 <- threshold(fit, method = "TSS", point = poipa, field_occurrence = "Observed")
  )
  expect_length(tr2$show_rasters(), 2)

})

# --- #
# Test validations
# Other generic functions in the package
test_that('Test validations', {

  suppressWarnings(
    requireNamespace("terra", quietly = TRUE)
  )
  # Dummy layers
  r_cont <- terra::rast(nrows = 10, ncols = 10, res = 0.05, xmin = -1.5, xmax = 1.5, ymin = -1.5, ymax = 1.5,
                    vals = rnorm(3600, 10, 2))
  names(r_cont) <- "cont"
  r_disc <- terra::rast(nrows = 10, ncols = 10, res = 0.05, xmin = -1.5, xmax = 1.5, ymin = -1.5, ymax = 1.5,
                        vals = rbinom(3600,1,.5))
  names(r_disc) <- "disc"

  # Make some dummy points
  s_cont <- terra::spatSample(r_cont, size = 100, as.points = TRUE) |> sf::st_as_sf()
  s_disc <- terra::spatSample(r_disc, size = 100, as.points = TRUE) |> sf::st_as_sf()
  expect_s3_class(s_cont, "sf")
  expect_true(hasName(s_cont, "cont"))
  expect_s3_class(s_disc, "sf")
  expect_true(hasName(s_disc, "disc"))

  # Continious tests
  expect_error(validate(r_cont, method = "disc"))
  suppressMessages(
    expect_error(validate(r_cont, method = "discrete", point = s_cont))
  )
  test <- validate(r_cont, method = "cont", point = s_cont, point_column = "cont")
  expect_s3_class(test, "data.frame")
  expect_true("rmse" %in% test$metric)

  # --- #
  # Discrete
  suppressMessages(
    expect_error(
      validate(r_disc, method = "disc", point = s_cont, point_column = "cont")
    )
  )

  # Skip if not installed
  skip_if_not_installed("modEvA")

  test <- validate(r_disc, method = "disc", point = s_disc, point_column = "disc")
  expect_s3_class(test, "data.frame")
  expect_true("f1" %in% test$metric)

})
