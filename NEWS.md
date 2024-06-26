# ibis.iSDM 0.1.4 (current dev branch)

#### New features

#### Minor improvements and bug fixes
* :bug: fix to `engine_gdb` also to support non-linear smooth functions (again).
* Small fix to support deprecated `field_occurrence` field in `validate` for convenience.
* :bug: fix that prevented `BART` models to be saved/loaded from disk #127.
* :bug: fixes related to `factor` handling for all engines.
* :bug: fixes related to `is_comparable_raster` and `add_predictors`/`add_predictors_range` #130
* :bug: fix related to `partial` for `engine_gdb` and priors

# ibis.iSDM 0.1.3

#### New features
* Add functions that creates HTML file base on `DistributionModel`.
* Added new engine `engine_scampr()` for model-based integration.
* Allow projection of models using `method_integration = "predictor"`

#### Minor improvements and bug fixes
* Small fixes to ensure `boruta` filtering works (again)?
* Small fix to parameter in `train()` #102 @jeffreyhanson
* Small helper function for combining 2 different formula objects `combine_formulas()`
* Small bug fixes dealing with `scenario()` projections and limits, plus unit tests #104
* Bug fixes with adding `predictor_derivate()` to scenario predictors and added unit tests #106
* Several fixes related to different engines and priors.
* Changed default output for netcdf files to multidimensional arrays #109
* :fire: hot fixes for scenario scaling and normalization issue #113
* :bug: fix so that projection works with different extents than used for inference.

# ibis.iSDM 0.1.2

#### New features
* Switched object structure to `R6` throughout for improved data and memory handling #44
* Implemented a convenience function ro remove biodiversity datasets (`rm_biodiversity()`).

#### Minor improvements and bug fixes
* Added a logical parameter to `ensemble()` enabling compositing of thresholds if set #84
* Support of multi-band rasters in `ensemble()` for convenience.
* Fix of bug in `threshold()` for supplied point data and improved error messages.
* Cleaner docs and structure
* Adding `wrap_model`/`unwrap_model` functions
* Added default parameters for all ibis specific options #90
* Changing behaviour of weights in `engine_inlabru()` #93

# ibis.iSDM 0.1.1

#### New features
* Added default `engine_glm()` for dependency-free inference and projection.
* Harmonized controls settings and added option to contrain extrapolation `add_control_extrapolation()`
* Adding a function for temporal interpolation of predictors #52

#### Minor improvements and bug fixes
* Minor corrective fixes and additions to `add_offset()`.
* Switch to `engine_glm()` in many of the unittests for better coverage.
* Several bug fixes and improvements in `thin_observations()`
* `global`, `probs`, and `centers` argument for better control of `thin_observations()`
* Harmonization of parameters for `spartial()` and addressing #80

# ibis.iSDM 0.1.0

#### New features
* Added a small convenience wrapper to add model outputs to another model `add_predictors_model()`
* Started adding mechanistic SDM vignette #67
* Wrapper for *steps* implemented via `simulate_population_steps()` #68

#### Minor improvements and bug fixes
* Added R-universe installation option as alternative to github #38
* Minor bug fixes in `scenario()` object, and MigClim and Kissmig wrappers.
* Bug fix related to CRS classes of sp and sf
* Bug fix related to blas.num.threads
* Bug fix that crashed `write_summary()` outputs when no prediction was made.
* Bug fix related to CRS in `engine_inla()`
* Bug fix in `engine_stan()` related to background layer
* Class of biodiversity data is identical for PO and PA
* Bug fix in `built_formula_glmnet()` and response
* Bug fix in `built_formula_gdb()` and response
* Each model$biodiversity stores only predictors of current ID
* Bug fix in `built_formula_inla()` for INLABRU

# ibis.iSDM 0.0.9

#### New features
* Added new vignette on available functions for data preparation #67
* Addition of small `mask()` function that emulates the for `terra`.

#### Minor improvements and bug fixes
* Small fix to `ensemble()` so that ensembles of future scenarios use correct standardization.
* Small fix to `threshold()` now returning threshold values correctly. 
* Bug fix and error catching to `distribution()` and `ensemble_partial()`,`ensemble_spartial()`
* Further checks added to `check()` #45
* Small fix to `alignRasters()`.
* Small fix to harmonize field_column throughout.
* Improved error messages and handling of formula's.

# ibis.iSDM 0.0.8

#### New features
* Implemented min size constraint (`add_constraint_minsize()`) #56
* Added a function for estimating partial effects of ensembles `ensemble_spartial()`.

#### Minor improvements and bug fixes
* Added warnings and checks for missing crs in supplied layers #65
* Smaller bug and code harmonizations to `ensemble_partial()`, `partial()` and `spartial()`. 
* Smaller bug fixes to `threshold()` in `scenario()` projections.
* Improved error messages in several functions.
* Further documentation fixes towards CRAN submission #38
* Allow to specify location of biodiversity point records in `threshold()`.

# ibis.iSDM 0.0.7

#### New features
* Added method proximity to `add_control_bias()` to place lower weights on points closer to another.
* Added helper functions `get_data()` and the option to apply `threshold()` directly on BiodiversityScenarios.
* Added centroid function to BiodiversityScenarios and DistributionModels #29

#### Minor improvements and bug fixes
* Add Error message for background data of different units easier to understand.
* Added warning message to the threshold creation to use independent data where possible.
* Fixed min.cv bug in `threshold()` introduced by #17
* Fixed `add_offset()` function now also allowing sf objects as input.
* Fixed bug with writing outputs in `write_output()`
* Fixed a bug so that prediction limits work correctly again (`distribution(...,lim = x)`)

# ibis.iSDM 0.0.6

#### New features
* `partial_density()` function implemented  #57
* Re-specification of limits with implementation of minimum convex polygon limits to `distribution()`.
* Added `check()` function for assessing assumptions and fits for various objects #45
* Added minor internal helper functions to duplicate `stars` objects via `st_rep`.
* Implemented local limiting factor function (`limiting()`) #37

#### Minor improvements and bug fixes
* Further smaller documentation fixes towards a CRAN submission #38
* Bug fix to method `buffer` in pseudo-absence settings.
* Minor bug fixes to `ensemble()` uncertainty calculations.

# ibis.iSDM 0.0.5

#### New features
* Addition of 5 parameter logistic curve offsets with parameter search to `add_offset()`.

#### Minor improvements and bug fixes
* Further smaller documentation fixes towards a CRAN submission #38
* Bug with with `write_model()`, now converting `terra` objects to `data.frame` between import/export.
* Smaller bug fixes, for example in `similarity()`, addition of variable name sanitization to predictors by default.

# ibis.iSDM 0.0.4

#### Minor improvements and bug fixes
* Smaller bug fixes with regards to writing outputs and adding pseudo-absences.
* Added short convenience function to convert prediction outputs #48
* Converted from `raster` to `terra` #17
* Updated and added further unit checks and tests

# ibis.iSDM 0.0.3

#### New features
* Aded Boruta for iterative feature selection of predictor variables.

#### Minor improvements and bug fixes
* Removed Magittr dependency #41
* Smaller improvements to documentation and removing of CRAN preventing function calls.
* Made the separation from hyperparameter search functions clearer and added new option to filter highly correlated covariates via `train()`.

# ibis.iSDM 0.0.2

#### Minor improvements and bug fixes
* Smaller documentation fixes, including to make sure examples and returns are in all exported function documentations.
* Preparation for cran release #38, including fixing some common issues and checks.
* Some smaller bug fixes to `validate()` to make Boyce more robust.
* Change of the logo. Thanks to @elliwoto 
* Added warning to validate call for users to be aware of non-independent validation.
* Further fixes on github actions and tests by @mhesselbarth

# ibis.iSDM 0.0.1

* Initial public release version! Finding and fixing further bugs... 
