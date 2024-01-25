# Setting up a distribution model
test_that('Setting up a distribution model',{

  skip_if_not_installed('igraph')

  skip_on_cran()

  suppressWarnings( requireNamespace("terra", quietly = TRUE) )
  suppressWarnings( requireNamespace("sf", quietly = TRUE) )
  suppressWarnings( requireNamespace("igraph", quietly = TRUE) )

  options("ibis.setupmessages" = FALSE)
  # Background Raster
  background <- terra::rast(system.file('extdata/europegrid_50km.tif', package='ibis.iSDM',mustWork = TRUE))
  # Get test species
  virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'points',quiet = TRUE)
  virtual_range <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'range',quiet = TRUE)
  # Get list of test predictors
  ll <- list.files(system.file('extdata/predictors',package = 'ibis.iSDM',mustWork = TRUE),full.names = T)
  # Load them as rasters
  predictors <- terra::rast(ll);names(predictors) <- tools::file_path_sans_ext(basename(ll))

  expect_equal(nrow(virtual_points), 208)
  expect_equal(terra::ncell(predictors), 7313)

  # Expect errors if there are non-unique data
  background2 <- background
  background2 <- background2 * terra::cellSize(background2)
  expect_error(x <- distribution(background2))

  # Now set them one up step by step
  x <- distribution(background)
  expect_s3_class(x,'BiodiversityDistribution')
  expect_s3_class(x$background,'sf')
  expect_error(x$biodiversity$get_data())
  expect_equal(x$biodiversity$length(),0)
  expect_type(x$show_background_info(),'list')
  expect_null(x$get_engine())
  expect_vector(names(x))

  # Add small error check by setting crs to NULL
  background2 <- background
  terra::crs(background2) <- NULL
  expect_error(x <- distribution(background2))

  # Now add one variable
  x <- x |> add_biodiversity_poipo(virtual_points,field_occurrence = 'Observed',name = 'Virtual points')
  expect_message(x$biodiversity,NA)
  expect_equal(x$biodiversity$length(),1)
  expect_equal(x$biodiversity$get_equations()[[1]],'<Default>')
  expect_true(is.Waiver(x$engine))
  expect_error(train(x)) # Try to solve without solver
  expect_s3_class(x$biodiversity, "BiodiversityDatasetCollection")

  # Apply a mask
  expect_no_error( x$biodiversity$mask(virtual_range) )
  x <- distribution(background) |> add_biodiversity_poipo(virtual_points,
                                                         field_occurrence = 'Observed',
                                                         name = 'Virtual points',docheck = FALSE)

  # And a range off
  invisible( suppressWarnings( suppressMessages(y <- x |> add_offset_range(virtual_range))) )
  expect_equal(y$get_offset(),'range_distance')
  expect_s4_class(y$offset,'SpatRaster')

  # remove again
  x <- y$rm_offset()
  expect_true(is.Waiver( x$get_offset() ) )

  # Add biodiversity data and remove again
  y <- x$clone(deep = TRUE)
  y <- rm_biodiversity(y, id = y$get_biodiversity_ids()[[1]])
  expect_equal(y$show_biodiversity_length(), 0)

  # Try also different bias controls
  expect_no_error(x |> add_control_bias(predictors$hmi_mean_50km,method = "partial"))
  expect_no_error(x |> add_control_bias(predictors$hmi_mean_50km,method = "partial",bias_value = 0))
  expect_no_error(y <- x |> add_control_bias(predictors$hmi_mean_50km,method = "offset",bias_value = 0))
  expect_equal(y$get_offset(),"hmi_mean_50km")
  expect_no_error(y <- x |> add_control_bias(method = "proximity",bias_value = 0))
  expect_length(y$control$bias_value, 2)

  # Add Predictors
  x <- x |> add_predictors(predictors)
  expect_s3_class(x$predictors, "PredictorDataset")
  expect_equal(x$predictors$length(),14)
  expect_true(is.vector(x$get_predictor_names()))
  # Try removing one
  x <- x |> rm_predictors('bio01_mean_50km')
  expect_equal(x$predictors$length(),13)
  expect_error( rm_predictors(x,'bio20_mean_50km') )
  # Finally select all predictors with CLC3
  n <- grep('CLC',x$get_predictor_names(),value = TRUE)
  x <- x |> sel_predictors(n)
  expect_equal(x$predictors$length(),5)
  expect_equal(x$get_predictor_names(), n)

  # Add brick object make derivatives
  x <- distribution(background) |> add_predictors(predictors$aspect_mean_50km, derivates = 'quadratic')
  testthat::expect_equal(x$predictors$length(),2)
  x <- distribution(background) |> add_predictors(predictors$aspect_mean_50km, derivates = 'hinge')
  testthat::expect_equal(x$predictors$length(),5)
  x <- distribution(background) |> add_predictors(predictors$aspect_mean_50km, derivates = 'bin')
  testthat::expect_equal(x$predictors$length(),5)
  x <- distribution(background) |> add_predictors(predictors$aspect_mean_50km, derivates = 'none')
  testthat::expect_equal(x$predictors$length(),1)
  x <- distribution(background) |> add_predictors(predictors$aspect_mean_50km, derivates = 'thresh')
  testthat::expect_equal(x$predictors$length(),4)
  # For interactions
  expect_error( x <- distribution(background) |> add_predictors(predictors$aspect_mean_50km, derivates = 'interaction') )
  expect_error( x <- distribution(background) |> add_predictors(predictors$aspect_mean_50km, derivates = 'interaction',int_variables = "aspect_mean_50km"))
  x <- distribution(background) |> add_predictors(predictors, derivates = 'interaction',int_variables = c("aspect_mean_50km", "bio01_mean_50km"))
  testthat::expect_equal(x$predictors$length(), terra::nlyr(predictors) + 1)

  # Try multiple on all
  x <- distribution(background) |> add_predictors(predictors, derivates = c('quadratic','hinge'))
  testthat::expect_equal(x$predictors$length(),84)

  # Interactions
  suppressMessages( expect_error(x |> add_predictors(predictors, derivates = "interaction")) )
  suppressMessages(y <- (x |> add_predictors(predictors, derivates = "interaction",int_variables = c(1,2)) ) )
  testthat::expect_s3_class(y, "BiodiversityDistribution")
  rm(y)

  suppressWarnings( x <- x |> engine_glm() )

  # Do a check
  suppressMessages( expect_no_error(check(x)) )

  # Mesh is not created yet
  expect_s3_class(x$engine$get_data("mesh"),'Waiver')
  expect_equal(x$engine$name,'<GLM>')
  expect_error(x$engine$calc_stack_poipo()) # Nothing to train on

  expect_s3_class(x$get_priors(),'Waiver')
  expect_s3_class(x$get_limits(), 'Waiver')

  # Add latent effect and see whether the attributes is changed
  y <- x |> add_latent_spatial(method = "spde")
  expect_vector( attr(y$get_latent(),'method'),'spde')

  # ---- #
  # Check that all engines can be added with default options
  # Also add a zonal layer for limits
  zones <- terra::as.factor( predictors$koeppen_50km )
  x <- distribution(background, limits = zones)
  expect_type(x$get_limits(), "list")
  expect_s3_class(x$get_limits()$layer, "sf")

  # Alternatively with MCP
  x <- distribution(background, limits_method = "mcp", mcp_buffer = 1000)
  expect_type(x$get_limits(), "list")
  expect_equal(x$get_limits()$limits_method, "mcp")
  expect_null( x$get_engine() )

  y <- x |> engine_bart()
  expect_equal( y$get_engine(), "<BART>")
  y <- x |> engine_breg()
  expect_equal( y$get_engine(), "<BREG>")
  y <- x |> engine_gdb()
  expect_equal( y$get_engine(), "<GDB>")
  y <- x |> engine_xgboost()
  expect_equal( y$get_engine(), "<XGBOOST>")

  # Normal x should still be none
  expect_null( x$get_engine() )

  # MJ: INLA call last to avoid errors upfront.
  skip_if_not_installed('INLA')
  y <- x |> engine_inla()
  expect_equal( y$get_engine(), "<INLA>")
  y <- x |> engine_inlabru()
  expect_equal( y$get_engine(), "<INLABRU>")

  # MH: skip if no cmd stan path can be found, only quick-and-dirty fix for now
  skip_if_not_installed("cmdstanr")
  skip_if(condition = tryCatch(expr = cmdstanr::cmdstan_path(), error = function(e) return(TRUE)),
          message = "No cmdstan path")

  y <- x |> engine_stan()
  expect_equal( y$get_engine(), "<STAN>")

})
