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
  tr1b <- threshold(r1, method = "fixed", value = .5,format = "norm")
  expect_s4_class(tr1b, "SpatRaster")
  tr1c <- threshold(r1, method = "fixed", value = .5,format = "perc")
  expect_s4_class(tr1c, "SpatRaster")
  expect_true(is.factor(tr1c))

  expect_error( tr2 <- threshold(r1, method = "perc", value = .1) )
  tr2 <- threshold(r1, method = "perc", value = .1, point = ss)
  expect_s4_class(tr2, "SpatRaster")

  tr3 <- threshold(r1, method = "mtp", point = ss)
  expect_s4_class(tr3, "SpatRaster")

  expect_error(tr4 <- threshold(r1, method = "min.cv", point = ss))
  tr4 <- threshold(r1, method = "min.cv",value = .2, point = ss)
  expect_s4_class(tr4, "SpatRaster")

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
