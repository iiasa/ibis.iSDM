# Test any other functions in ibis that do not fit into the other tests
# These will be particular important as we switch to terra in the future
test_that('Custom functions - Test gridded transformations and ensembles', {
  # skip_on_travis()
  # skip_on_cran()

  suppressWarnings(
    requireNamespace("terra", quietly = TRUE)
  )

  # --- #
  # Manipulating raster files #
  # create dummy predictions

  r1 <- terra::rast(nrows = 10, ncols = 10, res = 0.05, xmin = -1.5, xmax = 1.5, ymin = -1.5, ymax = 1.5, vals = rnorm(3600,mean = .5,sd = .1))
  r2 <- terra::rast(nrows = 10, ncols = 10, res = 0.05, xmin = -1.5, xmax = 1.5, ymin = -1.5, ymax = 1.5, vals = rnorm(3600,mean = .5,sd = .5))
  r3 <- terra::rast(nrows = 10, ncols = 10, res = 0.05, xmin = -1.5, xmax = 1.5, ymin = -1.5, ymax = 1.5, vals = rnorm(3600,mean = .2,sd = .7))

  expect_s4_class(r1, "SpatRaster")
  expect_true(is.Raster(r1))
  expect_true( is_comparable_raster(r1,r2) )

  # Predictor transformations
  t0 <- predictor_transform(r1, option = "none")
  expect_true(is_comparable_raster(r1,t0))
  expect_true(identical(r1,t0))

  # Create transformations for each and check that those are raster values with valid outputs
  t1 <- predictor_transform(r1, option = "norm")
  expect_s4_class(t1, "SpatRaster")
  expect_equal(terra::minmax(t1)[1,1], 0)
  expect_equal(terra::minmax(t1)[2,1], 1)

  # Scaling
  t2 <- predictor_transform(r1, option = "scale")
  expect_s4_class(t2, "SpatRaster")

  # PCA
  expect_error( predictor_transform(r1, option = "pca") )
  t3 <- predictor_transform(c(r1, r2), option = "pca")
  expect_s4_class(t3, "SpatRaster")
  t3b <- predictor_transform(c(r1, r2, r3), option = "pca",pca.var = 1)
  expect_s4_class(t3b, "SpatRaster")
  expect_equal(terra::nlyr(t3b), 3)

  # windsorization
  t4 <- predictor_transform(r1, option = "windsor")
  expect_s4_class(t4, "SpatRaster")
  expect_gt(terra::global(r1,"max"), terra::global(t4,"max"))

  # Percentile
  t5 <- predictor_transform(r1, option = "percentile")
  expect_length( terra::levels(t5)[[1]][,1] , 10)

  # Reverse jackknife
  t6 <- predictor_transform(r1, option = "revjack")
  expect_s4_class(t6, "SpatRaster") # Likely does not do much

  # ---- #
  # Check that homogenization works
  rr <- r1
  rr[sample(1:3600,100)] <- NA # introduce artificial NAs
  # Homogenize NA values among all the raster layers
  tt <- predictor_homogenize_na(c(rr, r2, r3))
  expect_true(anyNA(terra::values(tt[[2]])))
  expect_equal(which(terra::values(is.na(tt[[2]]))),
               which(terra::values(is.na(tt[[3]]))))

  # Align missing values
  tt <- predictor_homogenize_na(c(rr, r2, r3),fill = FALSE)
  expect_equal(object = terra::global(is.na(rr), "sum")[,1] ==100,
               expected = terra::global( is.na(tt[[3]]), "sum")[,1] ==100)

  # ---- #
  ## Create derivates
  # quad
  t1 <- predictor_derivate(r2, option = "quadratic")
  expect_equal(terra::global(t1, "max")[,1], (terra::global(r2, "max")^2)[,1] )

  # hinge
  suppressWarnings( t2 <- predictor_derivate(r2, option = "hinge") )
  expect_equal(terra::nlyr(t2), 4)
  # With different number of knots
  suppressWarnings( t2b <- predictor_derivate(r2, option = "hinge",nknots = 6) )
  expect_equal(terra::nlyr(t2b), 4*2)

  # Threshold
  suppressWarnings( t3 <- predictor_derivate(r2, option = "thresh") )
  expect_true(all(unique(t3)>=0),
              all(unique(t3)<=1))

  # Bins
  suppressWarnings( t4 <- predictor_derivate(r2, option = "bin") )
  expect_s4_class(t4, "SpatRaster")

  # Interactions
  expect_error( predictor_derivate(r2, option = "int") )
  s <- c(r1,r2,r3)
  names(s) <- c("lyra", "lyrb", "lyrc")
  suppressWarnings( t5 <- predictor_derivate(s, option = "int",int_variables = c("lyra", "lyrc")) )
  expect_s4_class(t5, "SpatRaster")

  # --- #
  # Finally do some ensemble calculations
  ex <- ensemble(r1, r2, r3, layer = "lyr.1")
  expect_equal(terra::nlyr(ex), 2)
  expect_lte( terra::global(ex, "max", na.rm = TRUE)[1,1], max( terra::global( c(r1, r2, r3), "max", na.rm = TRUE) ))
  expect_named(object = ex, expected = c("ensemble_lyr.1", "cv_lyr.1"))

  ex <- ensemble(r1, r2, r3, layer = "lyr.1", normalize = TRUE)
  expect_lte( terra::global(ex, "max")[1,1], 1)

  ex_sd <- ensemble(r1, r2, r3, layer = "lyr.1", uncertainty = "sd")
  ex_range <- ensemble(r1, r2, r3, layer = "lyr.1", uncertainty = "range")

  expect_named(object = ex_sd, expected = c("ensemble_lyr.1", "sd_lyr.1"))
  expect_named(object = ex_range, expected = c("ensemble_lyr.1", "range_lyr.1"))

  expect_error(ensemble(r1, r2, r3, layer = "lyr.1", uncertainty = "pca"),
               regexp = "Currently, uncertainty = 'pca' is not implemented for SpatRaster input.")

  # Check centroid calculation
  expect_s3_class(raster_centroid(r1), "sf")
  expect_s3_class(raster_centroid(r1,patch = TRUE), "sf")
  # Make binary
  b <- r1; b[b<.5] <- 0; b[b>0] <- 1
  expect_s3_class(raster_centroid(b,patch = FALSE), "sf")
  expect_true(nrow( raster_centroid(b,patch = TRUE) )>1)

  # --- #
  # Calculate threshold on new rasters
  o <- ex[[1]]
  pp <- terra::spatSample(o,20,as.points = TRUE) |> sf::st_as_sf()
  pp[[1]] <- sample(c(0,1),size = nrow(pp),replace = TRUE)
  names(pp)[1] <- "observed"
  # Dummy points

  expect_error(tr <- threshold(o))
  expect_no_error(tr <- threshold(o,value = .2))
  expect_s4_class(tr, "SpatRaster")
  # Return only threshold
  expect_no_error(tr <- threshold(o,value = .2,return_threshold = TRUE))
  expect_type(tr, "double")

  expect_no_error(tr <- threshold(o,method = "perc",point = pp,return_threshold = TRUE))
  expect_type(tr, "double")

  # --- #
})



# ---- #
# Other generic functions in the package
test_that('Test other generic functions', {

  # Options working?
  expect_type(ibis_options(), "list")

  # Colours working
  expect_type(ibis.iSDM:::ibis_colours, "list")

  # New id
  expect_no_error(id <- new_id() )
  expect_true(is.Id(id))
  expect_type(as.character(id), "character")
  expect_true( as.Id(as.character(id)) |> is.Id() )

  # New waivers
  expect_no_error(new <- new_waiver())
  expect_true(is.Waiver(new))

  # New form
  expect_no_error(form <- y ~ x)
  expect_true(is.formula(form))

  # Form functions
  expect_type(logistic(5), "double")
  expect_warning(logit(5))
  suppressWarnings( expect_true( is.nan(logit(5)) ) )
  expect_length(logisticRichard(rpois(10,.5)), 10)

})

# ---- #
# Other generic functions in the package
test_that('Test pseudo-absence options', {

  abs <- pseudoabs_settings()

  expect_s3_class(abs, "Settings")
  expect_equal(abs, ibis_options()$ibis.pseudoabsence) # Should be identical to default

  # Wrong method
  expect_error(pseudoabs_settings(method = "cool"))

  # Add custom options
  abs <- pseudoabs_settings(nrpoints = 1000, min_ratio = 1, method = "buffer",inside = FALSE, buffer_distance = 100)
  expect_s3_class(abs, "Settings")
  expect_equal(abs$get("nrpoints"), 1000)

})

# ---- #
# Test data
test_that('Test data preparation convenience functions', {

  # Load files
  background <- terra::rast(system.file('extdata/europegrid_50km.tif', package='ibis.iSDM',mustWork = TRUE))
  virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'points',quiet = TRUE)
  virtual_range <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'range',quiet = TRUE)
  # Get list of test predictors
  ll <- list.files(system.file('extdata/predictors',package = 'ibis.iSDM',mustWork = TRUE),full.names = T)
  predictors <- terra::rast(ll);names(predictors) <- tools::file_path_sans_ext(basename(ll))

  # --- #
  # Apply thinning methods
  pp1 <- thin_observations(df = virtual_points, background = background,
                                      method = "random", minpoints = 100,verbose = FALSE)
  expect_gt(nrow(pp1),0)
  # - #
  expect_error(pp2 <- thin_observations(df = virtual_points, background = background,
                           method = "bias", verbose = FALSE))
  pp2 <- thin_observations(df = virtual_points, background = background,
                           env = predictors$hmi_mean_50km,
                           method = "bias", verbose = FALSE)
  expect_gt(nrow(pp2),0)
  # - #
  pp3 <- thin_observations(df = virtual_points, background = background,
                           zones = as.factor(predictors$koeppen_50km),
                           method = "zones", verbose = FALSE)
  expect_gt(nrow(pp3),0)
  # - #
  pp4 <- thin_observations(df = virtual_points, background = background,
                           env = predictors,
                           method = "environmental", verbose = FALSE)
  expect_gt(nrow(pp4),0)
  # pp5 <- thin_observations(df = virtual_points, background = background,
  #                          env = predictors,
  #                          method = "spatial", verbose = FALSE)
  # expect_gt(nrow(pp5),0)

  # --- #
  # Aggregation functionalities
  pp1 <- aggregate_observations2grid(virtual_points, background,
                                                field_occurrence = "Observed")
  expect_s3_class(pp1, "sf")
  expect_gt(nrow(pp1), 1)
  expect_true(any(pp1$observed>1))


  # -- #
  # Zonal masking functions
  suppressMessages(
    lim <- create_zonaloccurrence_mask(df = virtual_points, buffer_width = 3,
                                     template = background)
  )
  expect_s4_class(lim, "SpatRaster")

  # --- #
  o <- st_kde(points = virtual_points,background = background)
  expect_s4_class(o, "SpatRaster")

  # --- #
  # Correct variable names
  vars <- c("Climate-temperature2015", "Elevation__sealevel", "Landuse.forest..meanshare")
  newvars <- sanitize_names(vars)
  expect_length(newvars, 3)
  expect_type(newvars, "character")

  # --- #
  # Similarity function test
  x <- distribution(background) |> add_predictors(predictors)
  expect_error( m1 <- similarity(x, method = "mess") )
  x <- x |> add_biodiversity_poipo(virtual_points)
  m1 <- similarity(x, method = "mess", plot = FALSE)
  expect_true(is.Raster(m1))
  expect_true(any(is.factor(m1)))
  m2 <- similarity(x, method = "nt", plot = FALSE)
  expect_true(is.Raster(m2))
  expect_true(any(is.factor(m2)))

})

