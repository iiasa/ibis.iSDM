# Further tests for model fits
test_that('Add further tests for model fits', {

  # Set to verbose
  options("ibis.setupmessages" = FALSE)

  # Load data
  # Background Raster
  background <- terra::rast(system.file('extdata/europegrid_50km.tif', package='ibis.iSDM',mustWork = TRUE))
  # Get test species
  virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'points',quiet = TRUE)
  virtual_range <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM',mustWork = TRUE),'range',quiet = TRUE)
  # Get list of test predictors
  ll <- list.files(system.file('extdata/predictors/',package = 'ibis.iSDM',mustWork = TRUE),full.names = T)
  # Load them as rasters
  predictors <- terra::rast(ll);names(predictors) <- tools::file_path_sans_ext(basename(ll))

  # Add pseudo absence
  abs <- pseudoabs_settings(nrpoints = 0,min_ratio = 1,method = "mcp")
  suppressMessages(
    virtual_points2 <- add_pseudoabsence(virtual_points,template = background, field_occurrence = "Observed", settings = abs)
  )

  # Create testing and training data
  ind <- sample(1:nrow(virtual_points2), 70)
  train_data <- virtual_points2[-ind,]
  test_data <- virtual_points2[ind,]

  # Now set them one up step by step
  x <- distribution(background) |>
    add_predictors(predictors, transform = 'none',derivates = 'none') |>
    engine_glm()

  # Train 2 model
  suppressWarnings(
    mod <- train(x |> add_biodiversity_poipa(train_data, field_occurrence = 'Observed',
                                             name = 'Virtual points',docheck = F),
                 "test", inference_only = FALSE, only_linear = TRUE, varsel = "none", verbose = FALSE)
  )
  suppressWarnings(
    mod_poipo <- train(x |> add_biodiversity_poipo(virtual_points, field_occurrence = 'Observed',
                                                   name = 'Virtual points',docheck = F),
                 "test", inference_only = FALSE, only_linear = TRUE, varsel = "none", verbose = FALSE)
  )
  expect_s4_class(mod$get_data(), "SpatRaster")

  # Threshold with independent data
  suppressMessages(
    mod <- threshold(mod, method = "perc", format = "bin")
  )
  expect_no_error(mod_poipo <- threshold(mod_poipo, method = "perc", value = .2))
  expect_gt(mod$get_thresholdvalue(),0)
  expect_length(mod$show_rasters(), 2)

  # Make an ensemble and check that thresholds are also present
  expect_no_error( ens <- ensemble(mod, mod_poipo, method = "mean"))
  expect_length(names(ens), 3)

  # Summarize model
  expect_s3_class( summary(mod), "data.frame" )
  expect_s3_class( coef(mod), "data.frame" )

  # Validate
  val <- validate(mod, method = "cont")
  expect_s3_class(val, "data.frame")
  # Validate discrete
  val <- validate(mod, method = "disc")
  expect_s3_class(val, "data.frame")

  # Validate with withold data
  val <- validate(mod, method = "disc", point = test_data,point_column = "Observed")
  expect_s3_class(val, "data.frame")

  # ----------- #
  # Project with separate data
  pp <- mod$project(predictors |> as.data.frame(xy = TRUE, na.rm = FALSE))
  expect_s4_class(pp, "SpatRaster")

  # ----------- #
  # Create a suitability index
  o <- mod$calc_suitabilityindex()
  expect_s4_class(o, "SpatRaster")

  # ----------- #
  # Clip the projected data with an external mask
  expect_no_error( mod$mask(virtual_range) )

  # ----------- #
  # Other functions
  pp <- mod$get_centroid()
  expect_s3_class(pp, "sf")
  expect_length(mod$show_rasters(), 2)

  # ----------- #
  # Adding a range of latent constraints
  y <- x |> add_biodiversity_poipa(train_data, field_occurrence = 'Observed',
                                   name = 'Virtual points',docheck = F)
  suppressWarnings(
    mod1 <- train(y |> add_latent_spatial(method = "poly"), "test", inference_only = FALSE,
                  only_linear = TRUE, varsel = "none", verbose = FALSE)
  )
  expect_true(length( grep("spatialtrend", mod1$get_coefficients()[,1] ) )>0)
  suppressWarnings(
    mod2 <- train(y |> add_latent_spatial(method = "nnd"), "test", inference_only = FALSE,
                  only_linear = TRUE, varsel = "none", verbose = FALSE)
  )
  expect_true(length( grep("nearestpoint", mod2$get_coefficients()[,1] ) )>0)
  suppressWarnings(
    mod3 <- train(y |> add_latent_spatial(method = "kde"), "test", inference_only = FALSE,
                  only_linear = TRUE, varsel = "none", verbose = FALSE)
  )
  expect_true(length( grep("kde", mod3$get_coefficients()[,1] ) )>0)

  # Make an ensemble
  expect_no_error(
    o <- ensemble(mod1, mod2, mod3, method = "median",uncertainty = "range")
  )
  expect_s4_class(o, "SpatRaster")
  expect_length(names(o), 2) # Should be at maximum 2 layers

  # ----------- #
  # Partial stuff
  skip_if_not_installed("pdp")
  pp <- partial(mod, x.var = "bio19_mean_50km",plot = FALSE)
  expect_s3_class(pp, "data.frame")

  # Spartial
  pp <- spartial(mod,x.var = "bio19_mean_50km",plot = FALSE)
  expect_s4_class(pp, "SpatRaster")

  # ----------- #
  # Write model outputs
  tf <- base::tempfile()
  expect_no_error(write_summary(mod, paste0(tf, ".rds")))
  expect_no_error(write_model(mod, paste0(tf, ".rds")))

})
