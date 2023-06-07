# Test scenario creation and constraints
test_that('Scenarios and constraints', {

  skip_if_not_installed('stars')
  skip_if_not_installed('glmnet')
  skip_on_travis()
  skip_on_cran()

  suppressWarnings( requireNamespace('glmnet', quietly = TRUE) )
  suppressWarnings( requireNamespace('igraph', quietly = TRUE) )
  suppressWarnings( requireNamespace('stars', quietly = TRUE) )

  options("ibis.setupmessages" = FALSE) # Be less chatty
  options("ibis.seed" = 1234)
  set.seed(1234)

  # Load data
  # Background Raster
  background <- terra::rast(system.file('extdata/europegrid_50km.tif', package='ibis.iSDM',mustWork = TRUE))
  # Get test species
  virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'points',quiet = TRUE)
  virtual_range <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'range',quiet = TRUE)

  # Convert points to pseudo-absence
  ab <- pseudoabs_settings(nrpoints = 1000, min_ratio = 1, method = "random")
  virtual_points <- add_pseudoabsence(virtual_points, field_occurrence = 'Observed', template = background, settings = ab)

  # Load present and future predictors
  ll <- list.files(system.file('extdata/predictors_presfuture/',package = 'ibis.iSDM',mustWork = TRUE),full.names = T)

  # Load the same files future ones
  suppressWarnings(
    pred_future <- stars::read_stars(ll) |> stars:::slice.stars('Time', seq(1, 86, by = 10))
  )
  sf::st_crs(pred_future) <- sf::st_crs(4326)

  pred_current <- stars_to_raster(pred_future, 1)[[1]]
  names(pred_current) <- names(pred_future)

  # Basic validity checks
  expect_length(pred_future, 9)
  expect_s3_class(pred_future, "stars")
  names(pred_future) <- names(pred_current) # Ensure that the names are identical, which is necessary for the projections

  # --------------- #
  # Fit a model and add a threshold to it
  fit <- distribution(background) |>
    add_biodiversity_poipa(virtual_points, field_occurrence = 'Observed', name = 'Virtual points') |>
    add_predictors(pred_current) |>
    engine_glmnet(alpha = 0) |>
    train("test", inference_only = FALSE, verbose = FALSE) |>
    threshold(method = "perc", value = .3)

  # Expect summary
  expect_s3_class(summary(fit), "data.frame")
  expect_s3_class(fit$show_duration(), "difftime")
  expect_equal(length(fit$show_rasters()), 2) # Should have prediction and threshold
  expect_s3_class(fit$get_equation(), "formula")

  # -- #
  expect_gt(fit$get_thresholdvalue(), 0)
  expect_s3_class(fit$settings, "Settings")
  expect_type(fit$model, "list")
  expect_length(fit$get_resolution(), 2)
  # --------------- #

  ## Now Create scenario objects ##
  sc <- scenario(fit)
  expect_s3_class(sc$get_data(), "Waiver")
  expect_type(sc$modelobject, "character")
  sc <- scenario(fit, copy_model = TRUE)
  expect_s3_class(sc$get_model(),"DistributionModel")# Model correctly inherited?
  expect_equal(sc$modelid, fit$model$id)

  # Add covariates in various transformations
  x <- sc |> add_predictors(pred_future, transform = "none")
  expect_length(x$get_predictor_names(), 9)
  x <- sc |> add_predictors(pred_future, transform = "scale")
  expect_length(x$get_predictor_names(), 9)
  x <- sc |> add_predictors(pred_future, transform = "norm")
  expect_length(x$get_predictor_names(), 9)
  # x <- sc |> add_predictors(pred_future, transform = "pca")
  # expect_length(x$get_predictor_names(), 9)

  expect_equal(x$get_predictor_names(), names(pred_current))
  expect_length(x$get_timeperiod(), 2)
  expect_gt(x$get_timeperiod()[2],2050) # This might fail if I try to reformat the date
  #  Check that predictors are right
  expect_s3_class(x$get_predictors()$get_data(), "stars")
  expect_s3_class(x$get_predictors()$get_data(df = TRUE), "data.frame")

  invisible( x$rm_predictors() )
  expect_length(x$get_predictor_names(), 9) # Properly inherited?
  x <- x$rm_predictors()
  expect_length(x$get_predictor_names(), 0) # Properly inherited?

  # Try and add current raster Layers for the projection
  obj <- pred_current
  # Set some Z values and correct projection
  terra::time(obj) <- rep(as.Date("2015-01-01"), terra::nlyr(obj))
  terra::set.crs(obj, terra::crs( "+proj=longlat +datum=WGS84") )
  expect_false(is.na( terra::crs(obj) ))

  x <- sc |> add_predictors(obj, transform = "none")
  expect_length(x$get_predictor_names(), 9)
  expect_equal(x$get_predictor_names(), names(obj))
  expect_lte(as.numeric( diff(x$get_timeperiod()) ), 1)
  # Test train
  mod <- x |> project()
  expect_s3_class(mod$get_data(), "stars")

  # Also check that it works with single SpatRaster layers
  x <- sc |> add_predictors(obj[[5]], transform = "none")
  expect_length(x$get_predictor_names(), 1)
  # Test train (should be an error as predictors are missing)
  expect_error( mod <- x |> project() )

  # Apply some transformations
  x <- sc |> add_predictors(obj, transform = "norm")
  expect_length(x$get_predictor_names(), 9)

  # Predict
  mod <- x |> project()
  expect_s3_class(mod$get_data(), "stars")

  # Make a first projection
  mod <- sc |> add_predictors(pred_future) |> project()
  suppressWarnings( expect_s3_class(summary(mod), "data.frame") )
  invisible(
    suppressWarnings( expect_s3_class(mod$calc_scenarios_slope(plot = FALSE), "stars") )
  )
  expect_length(mod$get_predictors()$get_time(), 9)

  # These will throw errors as we haven't added thresholds
  expect_error(mod$plot_relative_change())
  expect_error(mod$summary_beforeafter())
  # Now add threshold
  mod <- sc |> add_predictors(pred_future) |> threshold() |> project()
  expect_s3_class(mod$summary_beforeafter(), "data.frame")
  expect_true(inherits(mod$plot_relative_change(plot=FALSE), "SpatRaster"))

  # identical
  expect_equal(as.numeric(mod$get_threshold()), fit$get_thresholdvalue())
  expect_invisible(mod$verify())

  # --- #
  # Finally add the various constraints
  mod <- sc |> add_predictors(pred_future) |> threshold()
  expect_invisible(mod$verify())

  # Boundary
  mod1 <- mod |> add_constraint_boundary(virtual_range) |> project()
  expect_type(mod1$get_constraints(), "list")
  mod1b <- mod |> add_constraint(method = "boundary", layer = virtual_range)   # Generic constraint
  expect_type(mod1b$get_constraints(), "list")
  expect_equal(names(mod1$get_constraints()), names(mod1b$get_constraints()))

  # Dispersal simple
  expect_error(mod |> add_constraint_dispersal(method = "sdd_nexpkernel"))
  mod2 <- mod |> add_constraint_dispersal(method = "sdd_nexpkernel", value = 2e4) |> project()
  expect_type(mod2$get_constraints(), "list")
  # Add two constraints
  mod2b <- mod |> add_constraint_dispersal(method = "sdd_nexpkernel", value = 2e4) |>
    add_constraint_boundary(virtual_range) |> project()
  expect_type(mod2b$get_constraints(), "list")
  expect_length(mod2b$get_constraints(), 2)

  # Connectivity stuff
  res <- pred_current$urban
  mod2 <- mod |> add_constraint_connectivity(method = "resistance", resistance = res)
  expect_equal(names(mod2$get_constraints()), "connectivity")
  expect_true(is.Raster(mod2$get_constraints()$connectivity$params$resistance))

  # --- #
  # Check that stabilization works
  mods <- mod |> project(stabilize = TRUE)
  expect_invisible(mods$verify())
})
