# ibis.iSDM 0.0.6 (current dev branch)

### New features
* `partial_density` function implemented  #57
* Re-specification of limits with implementation of minimum convex polygon limits to `distribution`.
* Added `check` function for assessing assumptions and fits for various objects #45
* Added minor internal helper functions to duplicate `stars` objects via `st_rep`.
* Implemented local limiting factor function (`limiting`) #37

### Minor improvements and bug fixes
* Further smaller documentation fixes towards a CRAN submission #38
* Bug fix to method `buffer` in pseudo-absence settings.
* Minor bug fixes to `ensemble` uncertainty calculations.

# ibis.iSDM 0.0.5

### New features
* Addition of 5 parameter logistic curve offsets with parameter search to `add_offset`.

### Minor improvements and bug fixes
* Further smaller documentation fixes towards a CRAN submission #38
* Bug with with `write_model`, now converting `terra` objects to `data.frames` between import/export.
* Smaller bug fixes, for example in `similarity`, addition of variable name sanitization to predictors by default.

# ibis.iSDM 0.0.4

### Minor improvements and bug fixes
* Smaller bug fixes with regards to writing outputs and adding pseudo-absences.
* Added short convenience function to convert prediction outputs #48
* Converted from `raster` to `terra` #17
* Updated and added further unit checks and tests

# ibis.iSDM 0.0.3

### New features
* Aded Boruta for iterative feature selection of predictor variables.

### Minor improvements and bug fixes
* Removed Magittr dependency #41
* Smaller improvements to documentation and removing of CRAN preventing function calls.
* Made the separation from hyperparameter search functions clearer and added new option to filter highly correlated covariates via `train`.

# ibis.iSDM 0.0.2

### Minor improvements and bug fixes
* Smaller documentation fixes, including to make sure examples and returns are in all exported function documentations.
* Preparation for cran release #38, including fixing some common issues and checks.
* Some smaller bug fixes to `validate` to make Boyce more robust.
* Change of the logo. Thanks to @elliwoto 
* Added warning to validate call for users to be aware of non-independent validation.
* Further fixes on github actions and tests by @mhesselbarth

# ibis.iSDM 0.0.1

* Initial public release version! Finding and fixing further bugs... 
