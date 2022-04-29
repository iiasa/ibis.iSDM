# Setting up a distribution model
test_that('Setting up a distribution model',{
  testthat::skip_on_cran()

  options("ibis.setupmessages" = FALSE)
  # Background Raster
  background <- raster::raster(system.file('extdata/europegrid_50km.tif', package='ibis.iSDM'))
  # Get test species
  virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM'),'points',quiet = TRUE)
  virtual_range <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM'),'range',quiet = TRUE)
  # Get list of test predictors
  ll <- list.files(system.file('extdata/predictors',package = 'ibis.iSDM'),full.names = T)
  # Load them as rasters
  predictors <- raster::stack(ll);names(predictors) <- tools::file_path_sans_ext(basename(ll))

  # Now set them one up step by step
  x <- distribution(background)
  expect_s3_class(x,'BiodiversityDistribution')
  expect_s3_class(x$background,'sf')
  expect_error(x$biodiversity$get_data())
  expect_equal(x$biodiversity$length(),0)
  expect_type(x$show_background_info(),'list')
  expect_equal(x$get_engine(),'None')
  expect_vector(names(x))

  # Now add one variable
  x <- x %>% add_biodiversity_poipo(virtual_points,field_occurrence = 'Observed',name = 'Virtual points')
  expect_message(x$biodiversity,NA)
  expect_equal(x$biodiversity$length(),1)
  expect_equal(x$biodiversity$get_equations()[[1]],'<Default>')
  expect_true(is.Waiver(x$engine))
  expect_error(train(x)) # Try to solve without solver

  # And a range off
  x <- x %>% add_offset_range(virtual_range)
  expect_equal(x$get_offset(),'range_distance')
  expect_s4_class(x$offset,'Raster')
  # remove again
  x <- x$rm_offset()
  expect_true(is.Waiver( x$get_offset() ) )

  # Add Predictors
  x <- x %>% add_predictors(predictors)
  expect_equal(x$predictors$length(),14)
  expect_true(is.vector(x$get_predictor_names()))
  # Try removing one
  x <- x %>% rm_predictors('bio01_mean_50km')
  expect_equal(x$predictors$length(),13)
  expect_error( rm_predictors(x,'bio20_mean_50km') )
  # Finally select all predictors with CLC3
  n <- grep('CLC',x$get_predictor_names(),value = TRUE)
  x <- x %>% sel_predictors(n)
  expect_equal(x$predictors$length(),5)
  expect_equal(x$get_predictor_names(), n)

  # Add brick object make derivatives
  pb <- raster::brick(predictors)
  x <- distribution(background) %>% add_predictors(pb$aspect_mean_50km,derivates = 'quadratic')
  testthat::expect_equal(x$predictors$length(),2)
  x <- distribution(background) %>% add_predictors(pb, derivates = c('quadratic','hinge'))
  testthat::expect_equal(x$predictors$length(),84)

  x <- x %>% engine_inla()
  expect_s3_class(x$engine$data$mesh,'inla.mesh')
  expect_equal(x$engine$name,'<INLA>')
  expect_error(x$engine$calc_stack_poipo()) # Nothing to train on

  expect_s3_class(x$get_priors(),'Waiver')
  expect_null(x$get_limits())

  expect_type(x$engine$get_data('mesh.area'),'double')
  expect_gt(sum(x$engine$get_data('mesh.area')),800)

  # Add latent effect and see whether the attributes is changed
  y <- x %>% add_latent_spatial(method = "spde")
  expect_vector( attr(y$get_latent(),'method'),'spde')

  # ---- #
  # Check that all engines can be added with default options
  # Also add a zonal layer for limits
  zones <- raster::ratify( predictors$koeppen_50km )
  x <- distribution(background,limits = zones)
  expect_s3_class(x$get_limits(), "sf")

  y <- x %>% engine_bart()
  expect_equal( y$engine$name, "<BART>")
  y <- x %>% engine_breg()
  expect_equal( y$engine$name, "<BREG>")
  y <- x %>% engine_gdb()
  expect_equal( y$engine$name, "<GDB>")
  y <- x %>% engine_inla()
  expect_equal( y$engine$name, "<INLA>")
  y <- x %>% engine_inlabru()
  expect_equal( y$engine$name, "<INLABRU>")
  y <- x %>% engine_stan()
  expect_equal( y$engine$name, "<STAN>")
  y <- x %>% engine_xgboost()
  expect_equal( y$engine$name, "<XGBOOST>")

})
