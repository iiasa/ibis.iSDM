# Changelog

## ibis.iSDM 0.1.8 (current dev branch)

##### New features

##### Minor improvements and bug fixes

- Replacement of the `st_kde` function with an internal implementation
  of kernel density estimation in c++

## ibis.iSDM 0.1.7

##### New features

- Implementation of small ensemble of models for linear engines
  ([`add_control_esm()`](https://iiasa.github.io/ibis.iSDM/reference/add_control_esm.md))
  [\#141](https://github.com/iiasa/ibis.iSDM/issues/141)

##### Minor improvements and bug fixes

- Raise informative error if bias layer has negative values.
  [\#148](https://github.com/iiasa/ibis.iSDM/issues/148)
- Update
  [`engine_xgboost()`](https://iiasa.github.io/ibis.iSDM/reference/engine_xgboost.md)
  to work with latest development build and remove \[pdp\] as dependency
  [\#149](https://github.com/iiasa/ibis.iSDM/issues/149)
- Improvement of document grammar and readability. Few more vignette
  examples based on latest development.
- 🐛 Fixing of
  [`engine_inlabru()`](https://iiasa.github.io/ibis.iSDM/reference/engine_inlabru.md)
  to work with latest version
  [\#150](https://github.com/iiasa/ibis.iSDM/issues/150)
  [\#145](https://github.com/iiasa/ibis.iSDM/issues/145)
- 🐛 Fixes to
  [`add_control_bias()`](https://iiasa.github.io/ibis.iSDM/reference/add_control_bias.md)
  to have it working again.
- 🐛 Refactor code to avoid importing namespaces
  [\#95](https://github.com/iiasa/ibis.iSDM/issues/95)

## ibis.iSDM 0.1.6

##### New features

- Support for \[`data.frame`\] as predictors in
  [`add_predictors()`](https://iiasa.github.io/ibis.iSDM/reference/add_predictors.md)
  [\#136](https://github.com/iiasa/ibis.iSDM/issues/136)
- Convenience function to allow \[`data.frame`\] and \[`SpatRaster`\] to
  be supplied directly via
  [`project()`](https://iiasa.github.io/ibis.iSDM/reference/project.md)
  [\#136](https://github.com/iiasa/ibis.iSDM/issues/136)
- Small helper function to summarize values in a zone
  \[[`predictor_summarize_zones()`](https://iiasa.github.io/ibis.iSDM/reference/predictor_summarize_zones.md)\].
- New scenario projection constraint option in
  \[[`add_constraint_adaptability()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_adaptability.md)\]
  for (thermal) limitation of a variable
  [\#137](https://github.com/iiasa/ibis.iSDM/issues/137)

##### Minor improvements and bug fixes

- Fix `kissmig` to work with version 2.0
  [\#140](https://github.com/iiasa/ibis.iSDM/issues/140)
- IIASA internal functionalities such as preparation of GLOBIOM data
  have been transferred to [BNRTools](https://github.com/iiasa/BNRTools)
- Small bug fixed related to manual provision of scenario thresholds.
- [`predictor_filter()`](https://iiasa.github.io/ibis.iSDM/reference/predictor_filter.md)
  now also accepts \[`SpatRaster`\] objects as inputs.
- Small fix so that bounding box extent is correctly printed
  [\#138](https://github.com/iiasa/ibis.iSDM/issues/138)
- Small 🐛 fixes to `engine_gdb` in dev branch.

## ibis.iSDM 0.1.5

##### New features

- New visualization function
  [`nicheplot()`](https://iiasa.github.io/ibis.iSDM/reference/nicheplot.md)
  to visualize suitability across 2 axes
  [\#87](https://github.com/iiasa/ibis.iSDM/issues/87).
- Support for ‘modal’ value calculations in
  [`ensemble()`](https://iiasa.github.io/ibis.iSDM/reference/ensemble.md).
- Support for ‘superlearner’ in
  [`ensemble()`](https://iiasa.github.io/ibis.iSDM/reference/ensemble.md).
- Support for ‘kmeans’ derived threshold calculation in
  [`threshold()`](https://iiasa.github.io/ibis.iSDM/reference/threshold.md)
  and
  [`predictor_derivate()`](https://iiasa.github.io/ibis.iSDM/reference/predictor_derivate.md).
- Support for future processing streamlined. See FAQ section for
  instructions [\#18](https://github.com/iiasa/ibis.iSDM/issues/18).

##### Minor improvements and bug fixes

- Now overwriting temporary data by default in
  [`predictor_transform()`](https://iiasa.github.io/ibis.iSDM/reference/predictor_transform.md)
  and similar functions.
- Minor 🐛 fix related to misaligned thresholds and negative exponential
  kernels.
- 🔥 🐛 fix for scenario projections that use different grain sizes than
  for inference.

## ibis.iSDM 0.1.4

##### New features

- Support for carnying over latent spatial effects
  ([`add_latent_spatial()`](https://iiasa.github.io/ibis.iSDM/reference/add_latent_spatial.md))
  to
  [`scenario()`](https://iiasa.github.io/ibis.iSDM/reference/scenario.md)
  projections.
- Convenience functions to remove limits and controls
  [`rm_limits()`](https://iiasa.github.io/ibis.iSDM/reference/rm_limits.md)/[`rm_control()`](https://iiasa.github.io/ibis.iSDM/reference/rm_control.md)
  [\#121](https://github.com/iiasa/ibis.iSDM/issues/121)
- 🔥 Enable stars and multi-temporal SpatRaster zones for
  [`scenario()`](https://iiasa.github.io/ibis.iSDM/reference/scenario.md)
  and
  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  [\#121](https://github.com/iiasa/ibis.iSDM/issues/121)

##### Minor improvements and bug fixes

- 🐛 fix for support of factor x continuous variable interaction
  [\#131](https://github.com/iiasa/ibis.iSDM/issues/131)
- Renamed `add_control_extrapolation` to
  [`add_limits_extrapolation()`](https://iiasa.github.io/ibis.iSDM/reference/add_limits_extrapolation.md).
- 🐛 fix to `engine_gdb` also to support non-linear smooth functions
  (again).
- Small fix to support deprecated `field_occurrence` field in `validate`
  for convenience.
- 🐛 fix that prevented `BART` models to be saved/loaded from disk
  [\#127](https://github.com/iiasa/ibis.iSDM/issues/127).
- 🐛 fixes related to `factor` handling for all engines.
- 🐛 fixes related to `is_comparable_raster` and
  `add_predictors`/`add_predictors_range`
  [\#130](https://github.com/iiasa/ibis.iSDM/issues/130)
- 🐛 fix related to `partial` for `engine_gdb` and priors

## ibis.iSDM 0.1.3

##### New features

- Add functions that creates HTML file base on `DistributionModel`.
- Added new engine
  [`engine_scampr()`](https://iiasa.github.io/ibis.iSDM/reference/engine_scampr.md)
  for model-based integration.
- Allow projection of models using `method_integration = "predictor"`

##### Minor improvements and bug fixes

- Small fixes to ensure `boruta` filtering works (again)?
- Small fix to parameter in
  [`train()`](https://iiasa.github.io/ibis.iSDM/reference/train.md)
  [\#102](https://github.com/iiasa/ibis.iSDM/issues/102)
  [@jeffreyhanson](https://github.com/jeffreyhanson)
- Small helper function for combining 2 different formula objects
  [`combine_formulas()`](https://iiasa.github.io/ibis.iSDM/reference/combine_formulas.md)
- Small bug fixes dealing with
  [`scenario()`](https://iiasa.github.io/ibis.iSDM/reference/scenario.md)
  projections and limits, plus unit tests
  [\#104](https://github.com/iiasa/ibis.iSDM/issues/104)
- Bug fixes with adding
  [`predictor_derivate()`](https://iiasa.github.io/ibis.iSDM/reference/predictor_derivate.md)
  to scenario predictors and added unit tests
  [\#106](https://github.com/iiasa/ibis.iSDM/issues/106)
- Several fixes related to different engines and priors.
- Changed default output for netcdf files to multidimensional arrays
  [\#109](https://github.com/iiasa/ibis.iSDM/issues/109)
- 🔥 hot fixes for scenario scaling and normalization issue
  [\#113](https://github.com/iiasa/ibis.iSDM/issues/113)
- 🐛 fix so that projection works with different extents than used for
  inference.

## ibis.iSDM 0.1.2

##### New features

- Switched object structure to `R6` throughout for improved data and
  memory handling [\#44](https://github.com/iiasa/ibis.iSDM/issues/44)
- Implemented a convenience function ro remove biodiversity datasets
  ([`rm_biodiversity()`](https://iiasa.github.io/ibis.iSDM/reference/rm_biodiversity.md)).

##### Minor improvements and bug fixes

- Added a logical parameter to
  [`ensemble()`](https://iiasa.github.io/ibis.iSDM/reference/ensemble.md)
  enabling compositing of thresholds if set
  [\#84](https://github.com/iiasa/ibis.iSDM/issues/84)
- Support of multi-band rasters in
  [`ensemble()`](https://iiasa.github.io/ibis.iSDM/reference/ensemble.md)
  for convenience.
- Fix of bug in
  [`threshold()`](https://iiasa.github.io/ibis.iSDM/reference/threshold.md)
  for supplied point data and improved error messages.
- Cleaner docs and structure
- Adding `wrap_model`/`unwrap_model` functions
- Added default parameters for all ibis specific options
  [\#90](https://github.com/iiasa/ibis.iSDM/issues/90)
- Changing behaviour of weights in
  [`engine_inlabru()`](https://iiasa.github.io/ibis.iSDM/reference/engine_inlabru.md)
  [\#93](https://github.com/iiasa/ibis.iSDM/issues/93)

## ibis.iSDM 0.1.1

##### New features

- Added default
  [`engine_glm()`](https://iiasa.github.io/ibis.iSDM/reference/engine_glm.md)
  for dependency-free inference and projection.
- Harmonized controls settings and added option to contrain
  extrapolation `add_control_extrapolation()`
- Adding a function for temporal interpolation of predictors
  [\#52](https://github.com/iiasa/ibis.iSDM/issues/52)

##### Minor improvements and bug fixes

- Minor corrective fixes and additions to
  [`add_offset()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset.md).
- Switch to
  [`engine_glm()`](https://iiasa.github.io/ibis.iSDM/reference/engine_glm.md)
  in many of the unittests for better coverage.
- Several bug fixes and improvements in
  [`thin_observations()`](https://iiasa.github.io/ibis.iSDM/reference/thin_observations.md)
- `global`, `probs`, and `centers` argument for better control of
  [`thin_observations()`](https://iiasa.github.io/ibis.iSDM/reference/thin_observations.md)
- Harmonization of parameters for
  [`spartial()`](https://iiasa.github.io/ibis.iSDM/reference/spartial.md)
  and addressing [\#80](https://github.com/iiasa/ibis.iSDM/issues/80)

## ibis.iSDM 0.1.0

##### New features

- Added a small convenience wrapper to add model outputs to another
  model
  [`add_predictors_model()`](https://iiasa.github.io/ibis.iSDM/reference/add_predictors_model.md)
- Started adding mechanistic SDM vignette
  [\#67](https://github.com/iiasa/ibis.iSDM/issues/67)
- Wrapper for *steps* implemented via
  [`simulate_population_steps()`](https://iiasa.github.io/ibis.iSDM/reference/simulate_population_steps.md)
  [\#68](https://github.com/iiasa/ibis.iSDM/issues/68)

##### Minor improvements and bug fixes

- Added R-universe installation option as alternative to github
  [\#38](https://github.com/iiasa/ibis.iSDM/issues/38)
- Minor bug fixes in
  [`scenario()`](https://iiasa.github.io/ibis.iSDM/reference/scenario.md)
  object, and MigClim and Kissmig wrappers.
- Bug fix related to CRS classes of sp and sf
- Bug fix related to blas.num.threads
- Bug fix that crashed
  [`write_summary()`](https://iiasa.github.io/ibis.iSDM/reference/write_summary.md)
  outputs when no prediction was made.
- Bug fix related to CRS in
  [`engine_inla()`](https://iiasa.github.io/ibis.iSDM/reference/engine_inla.md)
- Bug fix in
  [`engine_stan()`](https://iiasa.github.io/ibis.iSDM/reference/engine_stan.md)
  related to background layer
- Class of biodiversity data is identical for PO and PA
- Bug fix in `built_formula_glmnet()` and response
- Bug fix in `built_formula_gdb()` and response
- Each model\$biodiversity stores only predictors of current ID
- Bug fix in `built_formula_inla()` for INLABRU

## ibis.iSDM 0.0.9

##### New features

- Added new vignette on available functions for data preparation
  [\#67](https://github.com/iiasa/ibis.iSDM/issues/67)
- Addition of small
  [`mask()`](https://iiasa.github.io/ibis.iSDM/reference/mask.md)
  function that emulates the for `terra`.

##### Minor improvements and bug fixes

- Small fix to
  [`ensemble()`](https://iiasa.github.io/ibis.iSDM/reference/ensemble.md)
  so that ensembles of future scenarios use correct standardization.
- Small fix to
  [`threshold()`](https://iiasa.github.io/ibis.iSDM/reference/threshold.md)
  now returning threshold values correctly.
- Bug fix and error catching to
  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  and
  [`ensemble_partial()`](https://iiasa.github.io/ibis.iSDM/reference/ensemble_partial.md),[`ensemble_spartial()`](https://iiasa.github.io/ibis.iSDM/reference/ensemble_spartial.md)
- Further checks added to
  [`check()`](https://iiasa.github.io/ibis.iSDM/reference/check.md)
  [\#45](https://github.com/iiasa/ibis.iSDM/issues/45)
- Small fix to
  [`alignRasters()`](https://iiasa.github.io/ibis.iSDM/reference/alignRasters.md).
- Small fix to harmonize field_column throughout.
- Improved error messages and handling of formula’s.

## ibis.iSDM 0.0.8

##### New features

- Implemented min size constraint
  ([`add_constraint_minsize()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_minsize.md))
  [\#56](https://github.com/iiasa/ibis.iSDM/issues/56)
- Added a function for estimating partial effects of ensembles
  [`ensemble_spartial()`](https://iiasa.github.io/ibis.iSDM/reference/ensemble_spartial.md).

##### Minor improvements and bug fixes

- Added warnings and checks for missing crs in supplied layers
  [\#65](https://github.com/iiasa/ibis.iSDM/issues/65)
- Smaller bug and code harmonizations to
  [`ensemble_partial()`](https://iiasa.github.io/ibis.iSDM/reference/ensemble_partial.md),
  [`partial()`](https://iiasa.github.io/ibis.iSDM/reference/partial.md)
  and
  [`spartial()`](https://iiasa.github.io/ibis.iSDM/reference/spartial.md).
- Smaller bug fixes to
  [`threshold()`](https://iiasa.github.io/ibis.iSDM/reference/threshold.md)
  in
  [`scenario()`](https://iiasa.github.io/ibis.iSDM/reference/scenario.md)
  projections.
- Improved error messages in several functions.
- Further documentation fixes towards CRAN submission
  [\#38](https://github.com/iiasa/ibis.iSDM/issues/38)
- Allow to specify location of biodiversity point records in
  [`threshold()`](https://iiasa.github.io/ibis.iSDM/reference/threshold.md).

## ibis.iSDM 0.0.7

##### New features

- Added method proximity to
  [`add_control_bias()`](https://iiasa.github.io/ibis.iSDM/reference/add_control_bias.md)
  to place lower weights on points closer to another.
- Added helper functions
  [`get_data()`](https://iiasa.github.io/ibis.iSDM/reference/get_data.md)
  and the option to apply
  [`threshold()`](https://iiasa.github.io/ibis.iSDM/reference/threshold.md)
  directly on BiodiversityScenarios.
- Added centroid function to BiodiversityScenarios and
  DistributionModels
  [\#29](https://github.com/iiasa/ibis.iSDM/issues/29)

##### Minor improvements and bug fixes

- Add Error message for background data of different units easier to
  understand.
- Added warning message to the threshold creation to use independent
  data where possible.
- Fixed min.cv bug in
  [`threshold()`](https://iiasa.github.io/ibis.iSDM/reference/threshold.md)
  introduced by [\#17](https://github.com/iiasa/ibis.iSDM/issues/17)
- Fixed
  [`add_offset()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset.md)
  function now also allowing sf objects as input.
- Fixed bug with writing outputs in
  [`write_output()`](https://iiasa.github.io/ibis.iSDM/reference/write_output.md)
- Fixed a bug so that prediction limits work correctly again
  (`distribution(...,lim = x)`)

## ibis.iSDM 0.0.6

##### New features

- [`partial_density()`](https://iiasa.github.io/ibis.iSDM/reference/partial_density.md)
  function implemented
  [\#57](https://github.com/iiasa/ibis.iSDM/issues/57)
- Re-specification of limits with implementation of minimum convex
  polygon limits to
  [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md).
- Added
  [`check()`](https://iiasa.github.io/ibis.iSDM/reference/check.md)
  function for assessing assumptions and fits for various objects
  [\#45](https://github.com/iiasa/ibis.iSDM/issues/45)
- Added minor internal helper functions to duplicate `stars` objects via
  `st_rep`.
- Implemented local limiting factor function
  ([`limiting()`](https://iiasa.github.io/ibis.iSDM/reference/limiting.md))
  [\#37](https://github.com/iiasa/ibis.iSDM/issues/37)

##### Minor improvements and bug fixes

- Further smaller documentation fixes towards a CRAN submission
  [\#38](https://github.com/iiasa/ibis.iSDM/issues/38)
- Bug fix to method `buffer` in pseudo-absence settings.
- Minor bug fixes to
  [`ensemble()`](https://iiasa.github.io/ibis.iSDM/reference/ensemble.md)
  uncertainty calculations.

## ibis.iSDM 0.0.5

##### New features

- Addition of 5 parameter logistic curve offsets with parameter search
  to
  [`add_offset()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset.md).

##### Minor improvements and bug fixes

- Further smaller documentation fixes towards a CRAN submission
  [\#38](https://github.com/iiasa/ibis.iSDM/issues/38)
- Bug with with
  [`write_model()`](https://iiasa.github.io/ibis.iSDM/reference/write_model.md),
  now converting `terra` objects to `data.frame` between import/export.
- Smaller bug fixes, for example in
  [`similarity()`](https://iiasa.github.io/ibis.iSDM/reference/similarity.md),
  addition of variable name sanitization to predictors by default.

## ibis.iSDM 0.0.4

##### Minor improvements and bug fixes

- Smaller bug fixes with regards to writing outputs and adding
  pseudo-absences.
- Added short convenience function to convert prediction outputs
  [\#48](https://github.com/iiasa/ibis.iSDM/issues/48)
- Converted from `raster` to `terra`
  [\#17](https://github.com/iiasa/ibis.iSDM/issues/17)
- Updated and added further unit checks and tests

## ibis.iSDM 0.0.3

##### New features

- Aded Boruta for iterative feature selection of predictor variables.

##### Minor improvements and bug fixes

- Removed Magittr dependency
  [\#41](https://github.com/iiasa/ibis.iSDM/issues/41)
- Smaller improvements to documentation and removing of CRAN preventing
  function calls.
- Made the separation from hyperparameter search functions clearer and
  added new option to filter highly correlated covariates via
  [`train()`](https://iiasa.github.io/ibis.iSDM/reference/train.md).

## ibis.iSDM 0.0.2

##### Minor improvements and bug fixes

- Smaller documentation fixes, including to make sure examples and
  returns are in all exported function documentations.
- Preparation for cran release
  [\#38](https://github.com/iiasa/ibis.iSDM/issues/38), including fixing
  some common issues and checks.
- Some smaller bug fixes to
  [`validate()`](https://iiasa.github.io/ibis.iSDM/reference/validate.md)
  to make Boyce more robust.
- Change of the logo. Thanks to [@elliwoto](https://github.com/elliwoto)
- Added warning to validate call for users to be aware of
  non-independent validation.
- Further fixes on github actions and tests by
  [@mhesselbarth](https://github.com/mhesselbarth)

## ibis.iSDM 0.0.1

- Initial public release version! Finding and fixing further bugs…
