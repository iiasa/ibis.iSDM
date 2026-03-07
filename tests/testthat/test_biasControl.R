# Test bias control with GLM and GLMNET engines
# Specifically test the case where the bias layer is NOT already in the
# predictor set, which triggers the code path that adds it during train().
test_that('Bias control works with GLM and GLMNET', {

  skip_on_cran()

  options(ibis.setupmessages = FALSE)

  # Load data
  background <- terra::rast(system.file('extdata/europegrid_50km.tif', package='ibis.iSDM', mustWork = TRUE))
  virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg', package='ibis.iSDM', mustWork = TRUE), 'points', quiet = TRUE)
  ll <- list.files(system.file('extdata/predictors/', package = 'ibis.iSDM', mustWork = TRUE), full.names = TRUE)
  predictors <- terra::rast(ll)
  names(predictors) <- tools::file_path_sans_ext(basename(ll))

  # Use hmi as bias layer and remove it from the predictor set
  bias_layer <- predictors[['hmi_mean_50km']]
  predictors <- terra::subset(predictors, setdiff(names(predictors), 'hmi_mean_50km'))

  # Base distribution with scaled predictors and separate bias layer
  x <- distribution(background) |>
    add_biodiversity_poipo(virtual_points, field_occurrence = 'Observed', name = 'Virtual points') |>
    add_predictors(predictors, transform = 'scale', derivates = 'none') |>
    add_control_bias(layer = bias_layer, bias_value = 0)

  # --- GLM ---
  expect_no_error(
    suppressWarnings(
      mod_glm <- train(x |> engine_glm(), "test_bias_glm",
                       inference_only = FALSE, only_linear = TRUE,
                       varsel = "none", verbose = FALSE)
    )
  )
  expect_true(mod_glm$settings$get("bias_variable") == "hmi_mean_50km")

  # --- GLMNET ---
  skip_if_not_installed('glmnet')
  expect_no_error(
    suppressWarnings(
      mod_glmnet <- train(x |> engine_glmnet(), "test_bias_glmnet",
                          inference_only = FALSE, only_linear = TRUE,
                          varsel = "none", verbose = FALSE)
    )
  )
  expect_true(mod_glmnet$settings$get("bias_variable") == "hmi_mean_50km")
})
