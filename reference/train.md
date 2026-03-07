# Train the model from a given engine

This function trains a
[`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
model with the specified engine and furthermore has some generic options
that apply to all engines (regardless of type). See Details with regards
to such options.

Users are advised to check the help files for individual engines for
advice on how the estimation is being done.

## Usage

``` r
train(
  x,
  runname,
  filter_predictors = "none",
  optim_hyperparam = FALSE,
  inference_only = FALSE,
  only_linear = TRUE,
  method_integration = "predictor",
  keep_models = TRUE,
  aggregate_observations = TRUE,
  clamp = FALSE,
  verbose = getOption("ibis.setupmessages", default = TRUE),
  ...
)

# S4 method for class 'BiodiversityDistribution'
train(
  x,
  runname,
  filter_predictors = "none",
  optim_hyperparam = FALSE,
  inference_only = FALSE,
  only_linear = TRUE,
  method_integration = "predictor",
  keep_models = TRUE,
  aggregate_observations = TRUE,
  clamp = TRUE,
  verbose = getOption("ibis.setupmessages", default = TRUE),
  ...
)
```

## Arguments

- x:

  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  (i.e.
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md))
  object).

- runname:

  A [`character`](https://rdrr.io/r/base/character.html) name of the
  trained run.

- filter_predictors:

  A [`character`](https://rdrr.io/r/base/character.html) defining if and
  how highly correlated predictors are to be removed prior to any model
  estimation. Available options are:

  - `"none"` No prior variable removal is performed (Default).

  - `"pearson"`, `"spearman"` or `"kendall"` Makes use of pairwise
    comparisons to identify and remove highly collinear predictors
    (Pearson's `r >= 0.7`).

  - `"abess"` A-priori adaptive best subset selection of covariates via
    the `"abess"` package (see References). Note that this effectively
    fits a separate generalized linear model to reduce the number of
    covariates.

  - `"boruta"` Uses the `"Boruta"` package to identify non-informative
    features.

- optim_hyperparam:

  Parameter to tune the model by iterating over input parameters or
  selection of predictors included in each iteration. Can be set to
  `TRUE` if extra precision is needed (Default: `FALSE`).

- inference_only:

  By default the engine is used to create a spatial prediction of the
  suitability surface, which can take time. If only inferences of the
  strength of relationship between covariates and observations are
  required, this parameter can be set to `TRUE` to ignore any spatial
  projection (Default: `FALSE`).

- only_linear:

  Fit model only on linear baselearners and functions. Depending on the
  engine setting this option to `FALSE` will result in non-linear
  relationships between observations and covariates, often increasing
  processing time (Default: `TRUE`). How non-linearity is captured
  depends on the used engine.

- method_integration:

  A [`character`](https://rdrr.io/r/base/character.html) with the type
  of integration that should be applied if more than one
  [`BiodiversityDataset`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDataset-class.md)
  object is provided in `x`. Particular relevant for engines that do not
  support the integration of more than one dataset. Integration methods
  are generally sensitive to the order in which they have been added to
  the
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md)
  object. Available options are:

  - `"predictor"` The predicted output of the first (or previously
    fitted) models are added to the predictor stack and thus are
    predictors for subsequent models (Default).

  - `"offset"` The predicted output of the first (or previously fitted)
    models are added as spatial offsets to subsequent models. Offsets
    are back-transformed depending on the model family. This option
    might not be supported for every
    [`Engine`](https://iiasa.github.io/ibis.iSDM/reference/Engine-class.md).

  - `"interaction"` Instead of fitting several separate models, the
    observations from each dataset are combined and incorporated in the
    prediction as a factor interaction with the "weaker" data source
    being partialed out during prediction. Here the first dataset added
    determines the reference level (see Leung et al. 2019 for a
    description).

  - `"prior"` In this option we only make use of the coefficients from a
    previous model to define priors to be used in the next model. Might
    not work with any engine!

  - `"weight"` This option only works for multiple biodiversity datasets
    with the same type (e.g. `"poipo"`). Individual weight multipliers
    can be determined while setting up the model (**Note**: Default is
    1). Datasets are then combined for estimation and weighted
    respectively, thus giving for example presence-only records less
    weight than survey records. **Note** that this parameter is ignored
    for engines that support joint likelihood estimation.

- keep_models:

  [`logical`](https://rdrr.io/r/base/logical.html) if true and
  `method_integration = "predictor"`, all models are stored in the
  `.internal` list of the model object.

- aggregate_observations:

  [`logical`](https://rdrr.io/r/base/logical.html) on whether
  observations covering the same grid cell should be aggregated
  (Default: `TRUE`).

- clamp:

  [`logical`](https://rdrr.io/r/base/logical.html) whether predictions
  should be clamped to the range of predictor values observed during
  model fitting (Default: `FALSE`).

- verbose:

  Setting this [`logical`](https://rdrr.io/r/base/logical.html) value to
  `TRUE` prints out further information during the model fitting
  (Default: `FALSE`).

- ...:

  further arguments passed on.

## Value

A
[DistributionModel](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
object.

## Details

This function acts as a generic training function that - based on the
provided
[`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md)
object creates a new distribution model. The resulting object contains
both a `"fit_best"` object of the estimated model and, if
`inference_only` is `FALSE` a
[terra::SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
object named `"prediction"` that contains the spatial prediction of the
model. These objects can be requested via `object$get_data("fit_best")`.

Other parameters in this function:

- `"filter_predictors"` The parameter can be set to various options to
  remove highly correlated variables or those with little additional
  information gain from the model prior to any estimation. Available
  options are `"none"` (Default) `"pearson"` for applying a `0.7`
  correlation cutoff, `"abess"` for the regularization framework by Zhu
  et al. (2020), or `"RF"` or `"randomforest"` for removing the least
  important variables according to a randomForest model. **Note**: This
  function is only applied on predictors for which no prior has been
  provided (e.g. potentially non-informative ones).

- `"optim_hyperparam"` This option allows to make use of hyper-parameter
  search for several models, which can improve prediction accuracy
  although through the a substantial increase in computational cost.

- `"method_integration"` Only relevant if more than one
  [`BiodiversityDataset`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDataset-class.md)
  is supplied and when the engine does not support joint integration of
  likelihoods. See also Miller et al. (2019) in the references for more
  details on different types of integration. Of course, if users want
  more control about this aspect, another option is to fit separate
  models and make use of the
  [add_offset](https://iiasa.github.io/ibis.iSDM/reference/add_offset.md),
  [add_offset_range](https://iiasa.github.io/ibis.iSDM/reference/add_offset_range.md)
  and
  [ensemble](https://iiasa.github.io/ibis.iSDM/reference/ensemble.md)
  functionalities.

- `"clamp"` Boolean parameter to support a clamping of the projection
  predictors to the range of values observed during model training.

## Note

There are no silver bullets in (correlative) species distribution
modelling and for each model the analyst has to understand the
objective, workflow and parameters than can be used to modify the
outcomes. Different predictions can be obtained from the same data and
parameters and not all necessarily make sense or are useful.

## References

- Miller, D.A.W., Pacifici, K., Sanderlin, J.S., Reich, B.J., 2019. The
  recent past and promising future for data integration methods to
  estimate species’ distributions. Methods Ecol. Evol. 10, 22–37.
  https://doi.org/10.1111/2041-210X.13110

- Zhu, J., Wen, C., Zhu, J., Zhang, H., & Wang, X. (2020). A polynomial
  algorithm for best-subset selection problem. Proceedings of the
  National Academy of Sciences, 117(52), 33117-33123.

- Leung, B., Hudgins, E. J., Potapova, A. & Ruiz‐Jaen, M. C. A new
  baseline for countrywide α‐diversity and species distributions:
  illustration using \>6,000 plant species in Panama. Ecol. Appl. 29,
  1–13 (2019).

## See also

[engine_gdb](https://iiasa.github.io/ibis.iSDM/reference/engine_gdb.md),
[engine_xgboost](https://iiasa.github.io/ibis.iSDM/reference/engine_xgboost.md),
[engine_bart](https://iiasa.github.io/ibis.iSDM/reference/engine_bart.md),
[engine_inla](https://iiasa.github.io/ibis.iSDM/reference/engine_inla.md),
[engine_inlabru](https://iiasa.github.io/ibis.iSDM/reference/engine_inlabru.md),
[engine_breg](https://iiasa.github.io/ibis.iSDM/reference/engine_breg.md),
[engine_stan](https://iiasa.github.io/ibis.iSDM/reference/engine_stan.md),
[engine_glm](https://iiasa.github.io/ibis.iSDM/reference/engine_glm.md)

## Examples

``` r
 # Load example data
 background <- terra::rast(system.file('extdata/europegrid_50km.tif',
 package='ibis.iSDM',mustWork = TRUE))
 # Get test species
 virtual_points <- sf::st_read(system.file('extdata/input_data.gpkg',
 package='ibis.iSDM',mustWork = TRUE),'points',quiet = TRUE)

 # Get list of test predictors
 ll <- list.files(system.file('extdata/predictors/', package = 'ibis.iSDM',
 mustWork = TRUE),full.names = TRUE)
 # Load them as rasters
 predictors <- terra::rast(ll);names(predictors) <- tools::file_path_sans_ext(basename(ll))

 # Use a basic GLM to fit a SDM
 x <- distribution(background) |>
        # Presence-only data
        add_biodiversity_poipo(virtual_points, field_occurrence = "Observed") |>
        # Add predictors and scale them
        add_predictors(env = predictors, transform = "scale", derivates = "none") |>
        # Use GLM as engine
        engine_glm()
#> [Setup] 2026-03-07 11:17:15.400146 | Creating distribution object...
#> [Setup] 2026-03-07 11:17:15.401107 | Adding poipo dataset...
#> [Setup] 2026-03-07 11:17:15.492463 | Adding predictors...
#> [Setup] 2026-03-07 11:17:15.498464 | Transforming predictors...

 # Train the model, Also filter out co-linear predictors using a pearson threshold
 mod <- train(x, only_linear = TRUE, filter_predictors = 'pearson')
#> [Estimation] 2026-03-07 11:17:15.550872 | Collecting input parameters.
#> [Estimation] 2026-03-07 11:17:15.591254 | Filtering predictors via pearson...
#> [Estimation] 2026-03-07 11:17:15.597068 | Adding engine-specific parameters.
#> [Estimation] 2026-03-07 11:17:15.602683 | Engine setup.
#> [Estimation] 2026-03-07 11:17:15.719097 | Starting fitting: 805867b1
#> [Estimation] 2026-03-07 11:17:15.765491 | Starting prediction...
#> [Done] 2026-03-07 11:17:15.825371 | Completed after 0.27 secs
 mod
#> Trained GLM-Model (Unnamed run)
#>   Strongest summary effects:
#>      Positive: CLC3_112_mean_50km, CLC3_132_mean_50km, CLC3_211_mean_50km, ... (7)
#>      Negative: aspect_mean_50km, bio03_mean_50km, slope_mean_50km (3)
#>   Prediction fitted: yes
```
