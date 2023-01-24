# Test any other functions in ibis that do not fit into the other tests
# These will be particular important as we switch to terra in the future
test_that('Custom functions - Test gridded transformations and ensembles', {
  skip_on_travis()
  skip_on_cran()
  library(raster)

  # --- #
  # Manipulating raster files #
  # create dummy predictions

  r1 <- raster(nrows = 10, ncols = 10, res = 0.05, xmn = -1.5, xmx = 1.5, ymn = -1.5, ymx = 1.5, vals = rnorm(3600,mean = .5,sd = .1))
  r2 <- raster(nrows = 10, ncols = 10, res = 0.05, xmn = -1.5, xmx = 1.5, ymn = -1.5, ymx = 1.5, vals = rnorm(3600,mean = .5,sd = .5))
  r3 <- raster(nrows = 10, ncols = 10, res = 0.05, xmn = -1.5, xmx = 1.5, ymn = -1.5, ymx = 1.5, vals = rnorm(3600,mean = .2,sd = .7))

  expect_s4_class(r1, "RasterLayer")

  # Predictor transformations
  t0 <- predictor_transform(r1, option = "none")
  expect_true(raster::compareRaster(r1,t0))
  expect_true(identical(r1,t0))

  # Create transformations for each and cechk that those are raster values with valid outputs
  t1 <- predictor_transform(r1, option = "norm")
  expect_s4_class(t1, "RasterLayer")
  expect_equal(min(t1[]), 0)
  expect_equal(max(t1[]), 1)

  # Scaling
  t2 <- predictor_transform(r1, option = "scale")
  expect_s4_class(t2, "RasterLayer")

  # pca
  expect_error( predictor_transform(r1, option = "pca") )
  t3 <- predictor_transform(raster::stack(r1, r2), option = "pca")
  expect_s4_class(t3, "RasterLayer")
  t3b <- predictor_transform(raster::stack(r1, r2, r3), option = "pca",pca.var = 1)
  expect_s4_class(t3b, "RasterBrick")
  expect_equal(raster::nlayers(t3b), 3)

  # windsorization
  t4 <- predictor_transform(r1, option = "windsor")
  expect_s4_class(t4, "RasterLayer")
  expect_gt(raster::cellStats(r1,"max"), raster::cellStats(t4,"max"))

  # Percentile
  t5 <- predictor_transform(r1, option = "percentile")
  expect_length(unique(t5), 10)

  # Reverse jackknife
  t6 <- predictor_transform(r1, option = "revjack")
  expect_s4_class(t6, "RasterStack") # Likely does not do much

  # ---- #
  # Check that homogenization works
  rr <- r1
  rr[sample(1:3600,100)] <- NA # introduce artificial NAs
  # Homogenize NA values among all the raster layers
  tt <- predictor_homogenize_na(raster::stack(rr, r2, r3))
  expect_true(anyNA(raster::values(tt[[2]])))
  expect_equal(which(raster::values(is.na(tt[[2]]))),
               which(raster::values(is.na(tt[[3]]))))

  # Fill missing values
  # TODO: Implement
  # tt <- predictor_homogenize_na(raster::stack(rr, r2, r3),fill = TRUE)

  # ---- #
  ## Create derivates
  # quad
  t1 <- predictor_derivate(r2, option = "quadratic")
  expect_equal(raster::cellStats(t1, "max"), raster::cellStats(r2, "max")^2 )

  # hinge
  t2 <- predictor_derivate(r2, option = "hinge")
  expect_equal(raster::nlayers(t2), 4)
  # With different number of knots
  t2b <- predictor_derivate(r2, option = "hinge",nknots = 6)
  expect_equal(raster::nlayers(t2b), 4*2)

  # Threshold
  t3 <- predictor_derivate(r2, option = "thresh")
  expect_true(all(unique(t3)>=0),
              all(unique(t3)<=1))

  # Bins
  t4 <- predictor_derivate(r2, option = "bin")
  expect_s4_class(t4, "RasterStack")

  # Interactions
  expect_error( predictor_derivate(r2, option = "int") )
  s <- raster::stack(r1,r2,r3)
  names(s) <- c("lyra", "lyrb", "lyrc")
  t5 <- predictor_derivate(s, option = "int",int_variables = c("lyra", "lyrc"))
  expect_s4_class(t5, "RasterStack")

  # --- #
  # Finally do some ensemble calculations
  ex <- ensemble(r1, r2, r3,layer = "layer")
  expect_equal(raster::nlayers(ex), 2)
  expect_lte(raster::cellStats(ex, "max")[[1]], max(cellStats(raster::stack(r1, r2, r3), "max")))

  ex <- ensemble(r1, r2, r3, layer = "layer", normalize = TRUE)
  expect_lte(raster::cellStats(ex[[1]], "max"), 1)

})
