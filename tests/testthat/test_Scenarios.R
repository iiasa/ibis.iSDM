# Test scenario creation and constraints
test_that('Testing data prep functions for spatial-temporal data in stars', {

  skip_if_not_installed('geosphere')
  skip_if_not_installed('cubelyr')
  skip_if_not_installed('lwgeom')

  suppressWarnings( requireNamespace('stars', quietly = TRUE) )
  suppressWarnings( requireNamespace('cubelyr', quietly = TRUE) )

  options("ibis.setupmessages" = FALSE) # Be less chatty

  # Load some stars rasters
  ll <- list.files(system.file('extdata/predictors_presfuture/',
                               package = 'ibis.iSDM',
                               mustWork = TRUE),full.names = TRUE)

  # Load the same files future ones
  suppressWarnings(
    pred_future <- stars::read_stars(ll) |> dplyr::slice('Time', seq(1, 86, by = 10))
  )
  sf::st_crs(pred_future) <- sf::st_crs(4326)

  expect_true(is.stars(pred_future))

  # Do some aggregations
  new <- st_reduce(pred_future, c("primf","primn"), newname = "combined",fun = "sum")
  expect_true("combined" %in% names(new))
  expect_false("primn" %in% names(new))

  # Other aggregation methods
  expect_no_error(st_reduce(pred_future, c("primf","primn"), newname = "combined",fun = "mean") )
  expect_no_error(st_reduce(pred_future, c("primf","primn"), newname = "combined",fun = "multiply") )
  expect_no_error(st_reduce(pred_future, c("primf","primn"), newname = "combined",fun = "divide") )
  expect_no_error(st_reduce(pred_future, c("primf","primn"), newname = "combined",fun = "subtract") )
  expect_no_error(st_reduce(pred_future, c("primf","primn"), newname = "combined",
                            fun = "weighted.mean", weights = c("bio01","bio01")) ) # Weight by different layer

  # Stars to raster check
  expect_type(stars_to_raster(pred_future,which = 1), "list")
  # Also try the inverse
  o <- stars_to_raster(pred_future,which = 1)[[1]] |>
    raster_to_stars()
  expect_s3_class(o, "stars")
  expect_length(o, 9)

  # Add a new raster to stars object
  o <- stars_to_raster(pred_future,which = 1)[[1]]['bio01']
  suppressMessages(
    expect_length(st_add_raster(pred_future, o),10)
  )

  # Summarize a single variable
  expect_s3_class(
    summarise_projection(pred_future[1],fun = "mean"),
    "data.frame"
  )

  # Make a simple interpolation
  expect_no_error(
    new <- interpolate_gaps(pred_future, date_interpolation = "annual")
  )
  expect_length(new, 9)
  expect_length(stars::st_get_dimension_values(new, "time"), 81)

  # --- #
  # Make transformations
  test <- pred_future[1]
  expect_no_error(tt <- predictor_transform(test, option = 'scale'))
  expect_length(tt, 1)
  expect_no_error(tt <- predictor_transform(test, option = 'norm'))
  expect_length(tt, 1)
  expect_lte(max(tt[[1]], na.rm = TRUE), 1)
  expect_no_error(tt <- predictor_transform(test, option = 'windsor'))
  expect_length(tt, 1)
  expect_no_error(tt <- predictor_transform(test, option = 'percentile'))
  expect_length(tt, 1)

  # --- #
  # Create derivates of stars data for testing
  test <- pred_future[1]
  expect_length(test, 1)
  # Nothing
  expect_no_error(new <- predictor_derivate(test, option = "none"))
  expect_length(new, 1); expect_true(utils::hasName(new,"bio01"))
  # Quadratic transform
  expect_no_error(new <- predictor_derivate(test, option = "quad"))
  expect_length(new, 2); expect_true(utils::hasName(new,"bio01"))
  # Hinge transform
  expect_no_error(new <- predictor_derivate(test, option = "hinge", nknots = 4))
  expect_length(new, 5); expect_true(utils::hasName(new,"bio01"))
  # Threshold transform
  expect_no_error(new <- predictor_derivate(test, option = "thresh",nknots = 4))
  expect_length(new, 4);expect_true(utils::hasName(new,"bio01"))
  # Binning
  expect_no_error(new <- predictor_derivate(test, option = "bin",nknots = 4))
  expect_length(new, 2);expect_true(utils::hasName(new,"bio01"))
  # Interaction
  test <- pred_future[1:2]
  expect_length(test, 2)
  expect_error(new <- predictor_derivate(test, option = "interaction"))
  expect_no_error(new <- predictor_derivate(test, option = "interaction", int_variables = names(test)))
  expect_length(new, 3)
  expect_true(utils::hasName(new,"bio01"));expect_true(utils::hasName(new,"bio12"))

  # --- #
  # Create threshold
  expect_no_error(
    tr <- stars_to_raster(pred_future,1)[[1]] |>
      terra::subset("crops") |>
      threshold(method = "fixed",value = 100)
  )
  names(tr) <- "threshold"
  # Apply minimum size filter
  expect_no_error( st_minsize(tr, 2, unit = "pixel") )
  expect_no_error( st_minsize(tr, 2000, unit = "km2") )

  # Replicate
  # st_rep(pred_future[1], stars::st_dimensions(pred_future[1]), "Time")

})

# Test scenario creation and constraints
test_that('Scenarios and constraints', {

  skip_if_not_installed('geosphere')
  skip_if_not_installed('cubelyr')
  skip_if_not_installed('lwgeom')

  skip_on_cran()

  suppressWarnings( requireNamespace('igraph', quietly = TRUE) )
  suppressWarnings( requireNamespace('stars', quietly = TRUE) )
  suppressWarnings( requireNamespace('geosphere', quietly = TRUE) )
  suppressWarnings( requireNamespace('cubelyr', quietly = TRUE) )

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

  expect_true(is.stars(pred_future))

  pred_current <- stars_to_raster(pred_future, 1)[[1]]
  names(pred_current) <- names(pred_future)

  # Also test the reverse
  test <- raster_to_stars(pred_current)
  expect_length(test, 9)
  expect_equal(names(test), names(pred_current))
  expect_length(stars::st_get_dimension_values(test,"time"), 1)

  # Basic validity checks
  expect_length(pred_future, 9)
  expect_s3_class(pred_future, "stars")
  names(pred_future) <- names(pred_current) # Ensure that the names are identical, which is necessary for the projections

  # --------------- #
  # Fit a model and add a threshold to it
  fit <- distribution(background) |>
    add_biodiversity_poipa(virtual_points, field_occurrence = 'Observed', name = 'Virtual points') |>
    add_predictors(pred_current,transform = 'scale') |>
    engine_glm() |>
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
  suppressMessages( sc <- scenario(fit, copy_model = TRUE) )
  expect_s3_class(sc$get_model(),"DistributionModel")# Model correctly inherited?
  expect_equal(sc$modelid, fit$model$id)

  # Add covariates in various transformations
  x <- sc |> add_predictors(pred_future, transform = "none")
  expect_length(x$get_predictor_names(), 9)
  x <- sc |> add_predictors(pred_future, transform = "scale")
  expect_length(x$get_predictor_names(), 9)
  # Error as we used a different transform earlier
  expect_error( x <- sc |> add_predictors(pred_future, transform = "norm") )

  expect_equal(x$get_predictor_names(), names(pred_current))
  expect_length(x$get_timeperiod(), 2)
  expect_gt(x$get_timeperiod()[2],2050) # This might fail if I try to reformat the date
  #  Check that predictors are right
  expect_s3_class(x$get_predictors()$get_data(), "stars")
  expect_s3_class(x$get_predictors()$get_data(df = TRUE), "data.frame")

  expect_length(x$get_predictor_names(), 9)
  invisible( x$rm_predictors() )
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

  # Apply a mask
  expect_no_error( mod$mask(virtual_range) )
  mod <- x |> project() # Project anew

  # Calculate centroids
  expect_s3_class(mod$get_centroid(), "sf")

  # Apply a manual threshold and check that works
  mod <- threshold(mod, tr = .5)
  expect_true("threshold" %in% names(get_data(mod)))

  # Calculate centroid on threshold
  expect_s3_class(mod$get_centroid(), "sf")
  expect_true(nrow(mod$get_centroid(patch=TRUE))>1)

  # Also check that it works with single SpatRaster layers
  x <- sc |> add_predictors(obj[[5]], transform = "none")
  expect_length(x$get_predictor_names(), 1)
  # Test train (should be an error as predictors are missing)
  expect_error( mod <- x |> project() )

  # Apply some transformations
  expect_error( x <- sc |> add_predictors(obj, transform = "norm") )
  # This uses the correct transformation
  expect_no_error( x <- sc |> add_predictors(obj, transform = "scale") )
  expect_length(x$get_predictor_names(), 9)

  # Predict
  mod <- x |> project()
  # Get layer
  expect_s3_class(mod$get_data(), "stars")
  expect_s3_class(mod |> get_data(), "stars")

  # Make a first projection
  expect_no_error(
    mod <- sc |> add_predictors(pred_future) |> project()
  )
  suppressWarnings( expect_s3_class(summary(mod), "data.frame") )
  invisible(
    suppressWarnings( expect_s3_class(mod$calc_scenarios_slope(plot = FALSE), "stars") )
  )
  expect_length(mod$get_predictors()$get_time(), 9)

  # These will throw errors as we haven't added thresholds
  expect_error(mod$plot_relative_change())
  expect_error(mod$summary_beforeafter())

  # Now add threshold
  expect_no_error(
    mod <- sc |> add_predictors(pred_future) |> threshold() |> project()
  )
  expect_s3_class(mod$summary_beforeafter(), "data.frame")
  expect_true(inherits(mod$plot_relative_change(plot=FALSE), "SpatRaster"))

  # identical
  expect_equal(as.numeric(mod$get_threshold()), fit$get_thresholdvalue())
  expect_invisible(mod$verify())

  # Reapply a different threshold and check that works
  mod <- threshold(mod, tr = .25)
  expect_equal(as.numeric(mod$get_threshold()), .25)

  # ----------- #
  # Finally add the various constraints
  mod <- sc |> add_predictors(pred_future) |> threshold()
  expect_invisible(mod$verify())

  # Check summary
  mod0 <- mod |> project()
  expect_s3_class(mod0$summary_beforeafter(), 'data.frame')

  # Boundary constraint
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

  # Minimum size constraint
  mod3 <- mod |> threshold(value = .2) |>
    add_constraint_minsize(value = 4, unit = "pixel") |> project()
  expect_equal(names(mod3$get_constraints()), "min_size")
  expect_s3_class( mod3$get_data()['threshold'], "stars")

  # Threshold constraint
  expect_no_error(
    mod4 <- mod |> threshold(value = .2) |> add_constraint_threshold() |> project()
  )
  expect_equal(names(mod4$get_constraints()), "threshold")

  # --- #
  # Check that stabilization works
  mods <- mod |> project(stabilize = TRUE)
  expect_invisible(mods$verify())
})
