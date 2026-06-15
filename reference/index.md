# Package index

## Setting up and train models

Key functions for setting up species distribution models and adding
information to them. Start with distribution() and see articles examples
on how to build a model from there. Also includes functions to specify
priors for a model.

- [`distribution()`](https://iiasa.github.io/ibis.iSDM/reference/distribution.md)
  : Create distribution modelling procedure
- [`BARTPrior()`](https://iiasa.github.io/ibis.iSDM/reference/BARTPrior.md)
  : Create a tree-based split probability prior for BART
- [`BARTPriors()`](https://iiasa.github.io/ibis.iSDM/reference/BARTPriors.md)
  : Helper function when multiple variables are supplied for BART priors
- [`BREGPrior()`](https://iiasa.github.io/ibis.iSDM/reference/BREGPrior.md)
  : Create a new spike and slab prior for Bayesian generalized linear
  models
- [`BREGPriors()`](https://iiasa.github.io/ibis.iSDM/reference/BREGPriors.md)
  : Helper function when multiple variables are supplied for BREG priors
- [`GDBPrior()`](https://iiasa.github.io/ibis.iSDM/reference/GDBPrior.md)
  : Monotonic constrained priors for boosted regressions
- [`GDBPriors()`](https://iiasa.github.io/ibis.iSDM/reference/GDBPriors.md)
  : Helper function when multiple variables are supplied for GDB priors
- [`GLMNETPrior()`](https://iiasa.github.io/ibis.iSDM/reference/GLMNETPrior.md)
  : Regression penalty priors for GLMNET
- [`GLMNETPriors()`](https://iiasa.github.io/ibis.iSDM/reference/GLMNETPriors.md)
  : Helper function when multiple variables are supplied for GLMNET
  priors
- [`INLAPrior()`](https://iiasa.github.io/ibis.iSDM/reference/INLAPrior.md)
  : Create a new INLA prior
- [`INLAPriors()`](https://iiasa.github.io/ibis.iSDM/reference/INLAPriors.md)
  : Helper function when multiple variables and types are supplied for
  INLA priors
- [`STANPrior()`](https://iiasa.github.io/ibis.iSDM/reference/STANPrior.md)
  : Create a new STAN prior
- [`STANPriors()`](https://iiasa.github.io/ibis.iSDM/reference/STANPriors.md)
  : Helper function when multiple variables and types are supplied for
  Stan priors
- [`XGBInteractionPrior()`](https://iiasa.github.io/ibis.iSDM/reference/XGBInteractionPrior.md)
  : Create a new interaction prior for XGBoost
- [`XGBInteractionPriors()`](https://iiasa.github.io/ibis.iSDM/reference/XGBInteractionPriors.md)
  : Helper function when multiple interaction groups are supplied for
  XGBoost
- [`XGBPrior()`](https://iiasa.github.io/ibis.iSDM/reference/XGBPrior.md)
  : Create a new monotonic prior for boosted regressions
- [`XGBPriors()`](https://iiasa.github.io/ibis.iSDM/reference/XGBPriors.md)
  : Helper function when multiple variables are supplied for XGBoost
  priors
- [`Prior-class`](https://iiasa.github.io/ibis.iSDM/reference/Prior-class.md)
  [`Prior`](https://iiasa.github.io/ibis.iSDM/reference/Prior-class.md)
  : Base Prior class
- [`PriorList-class`](https://iiasa.github.io/ibis.iSDM/reference/PriorList-class.md)
  [`PriorList`](https://iiasa.github.io/ibis.iSDM/reference/PriorList-class.md)
  : List of Priors supplied to an class
- [`print(`*`<distribution>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<BiodiversityDistribution>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<BiodiversityDatasetCollection>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<BiodiversityDataset>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<PredictorDataset>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<DistributionModel>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<BiodiversityScenario>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<Prior>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<PriorList>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<Engine>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<Settings>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<Log>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<Id>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<tbl_df>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  : Print
- [`summary(`*`<distribution>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/summary.md)
  [`summary(`*`<DistributionModel>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/summary.md)
  [`summary(`*`<PredictorDataset>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/summary.md)
  [`summary(`*`<BiodiversityScenario>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/summary.md)
  [`summary(`*`<PriorList>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/summary.md)
  [`summary(`*`<Settings>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/summary.md)
  : Summarises a trained model or predictor object
- [`priors()`](https://iiasa.github.io/ibis.iSDM/reference/priors.md) :
  Creates a new PriorList object
- [`pseudoabs_settings()`](https://iiasa.github.io/ibis.iSDM/reference/pseudoabs_settings.md)
  : Settings for specifying pseudo-absence points within the model
  background
- [`train()`](https://iiasa.github.io/ibis.iSDM/reference/train.md) :
  Train the model from a given engine

## Add or modify data and parameters

Functions to add or modify data and parameters in a distribution object.
These can be used to add or remove biodiversity, covariates and priors
in various forms.

- [`add_biodiversity_poipa()`](https://iiasa.github.io/ibis.iSDM/reference/add_biodiversity_poipa.md)
  : Add biodiversity point dataset to a distribution object
  (presence-absence).

- [`add_biodiversity_poipo()`](https://iiasa.github.io/ibis.iSDM/reference/add_biodiversity_poipo.md)
  : Add biodiversity point dataset to a distribution object
  (presence-only)

- [`add_biodiversity_polpa()`](https://iiasa.github.io/ibis.iSDM/reference/add_biodiversity_polpa.md)
  : Add biodiversity polygon dataset to a distribution object
  (presence-absence)

- [`add_biodiversity_polpo()`](https://iiasa.github.io/ibis.iSDM/reference/add_biodiversity_polpo.md)
  : Add biodiversity polygon dataset to a distribution object
  (presence-only)

- [`add_constraint()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint.md)
  :

  Add a constraint to an existing `scenario`

- [`add_constraint_MigClim()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_MigClim.md)
  : Add constraints to the modelled distribution projection using the
  MigClim approach

- [`add_constraint_adaptability()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_adaptability.md)
  : Adds an adaptability constraint to a scenario object

- [`add_constraint_boundary()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_boundary.md)
  : Adds a boundary or zone constraint to a scenario object

- [`add_constraint_connectivity()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_connectivity.md)
  : Adds a connectivity constraint to a scenario object.

- [`add_constraint_dispersal()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_dispersal.md)
  :

  Add dispersal constraint to an existing `scenario`

- [`add_constraint_minsize()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_minsize.md)
  : Adds a size constraint on a scenario

- [`add_constraint_threshold()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_threshold.md)
  : Adds a threshold constraint to a scenario object

- [`add_control_bias()`](https://iiasa.github.io/ibis.iSDM/reference/add_control_bias.md)
  : Add a control to a BiodiversityModel object to control biases

- [`add_control_esm()`](https://iiasa.github.io/ibis.iSDM/reference/add_control_esm.md)
  : Add a control to a BiodiversityModel object to train ensembles of
  small models

- [`add_latent_spatial()`](https://iiasa.github.io/ibis.iSDM/reference/add_latent_spatial.md)
  : Add latent spatial effect to the model equation

- [`add_limits_extrapolation()`](https://iiasa.github.io/ibis.iSDM/reference/add_limits_extrapolation.md)
  : Add a control to a BiodiversityModel object to limit extrapolation

- [`add_log()`](https://iiasa.github.io/ibis.iSDM/reference/add_log.md)
  : Adds a log file to distribution or scenario object

- [`add_offset()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset.md)
  : Specify a spatial explicit offset

- [`add_offset_bias()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset_bias.md)
  : Specify a spatial explicit offset as bias

- [`add_offset_elevation()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset_elevation.md)
  : Specify elevational preferences as offset

- [`add_offset_range()`](https://iiasa.github.io/ibis.iSDM/reference/add_offset_range.md)
  : Specify a expert-based species range as offset

- [`add_predictor_elevationpref()`](https://iiasa.github.io/ibis.iSDM/reference/add_predictor_elevationpref.md)
  : Create lower and upper limits for an elevational range and add them
  as separate predictors

- [`add_predictor_range()`](https://iiasa.github.io/ibis.iSDM/reference/add_predictor_range.md)
  : Add a range of a species as predictor to a distribution object

- [`add_predictors()`](https://iiasa.github.io/ibis.iSDM/reference/add_predictors.md)
  : Add predictors to a Biodiversity distribution object

- [`add_predictors_globiom()`](https://iiasa.github.io/ibis.iSDM/reference/add_predictors_globiom.md)
  : Function to add GLOBIOM-DownScalr derived predictors to a
  Biodiversity distribution object

- [`add_predictors_model()`](https://iiasa.github.io/ibis.iSDM/reference/add_predictors_model.md)
  : Add predictions from a fitted model to a Biodiversity distribution
  object

- [`add_priors()`](https://iiasa.github.io/ibis.iSDM/reference/add_priors.md)
  : Add priors to an existing distribution object

- [`add_pseudoabsence()`](https://iiasa.github.io/ibis.iSDM/reference/add_pseudoabsence.md)
  : Add pseudo-absence points to a point data set

- [`set_priors(`*`<BiodiversityDistribution>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/set_priors-BiodiversityDistribution-method.md)
  : Add priors to an existing distribution object

- [`set_priors()`](https://iiasa.github.io/ibis.iSDM/reference/set_priors.md)
  : Add priors to an existing distribution object

- [`sel_predictors()`](https://iiasa.github.io/ibis.iSDM/reference/sel_predictors.md)
  :

  Select specific predictors from a
  [distribution](https://rdrr.io/r/stats/Distributions.html) object

- [`rm_biodiversity()`](https://iiasa.github.io/ibis.iSDM/reference/rm_biodiversity.md)
  :

  Remove specific BiodiversityDataset from a
  [distribution](https://rdrr.io/r/stats/Distributions.html) object

- [`rm_control()`](https://iiasa.github.io/ibis.iSDM/reference/rm_control.md)
  : Remove control from an existing distribution object

- [`rm_latent()`](https://iiasa.github.io/ibis.iSDM/reference/rm_latent.md)
  : Function to remove a latent effect

- [`rm_limits()`](https://iiasa.github.io/ibis.iSDM/reference/rm_limits.md)
  : Remove limits from an existing distribution object

- [`rm_offset()`](https://iiasa.github.io/ibis.iSDM/reference/rm_offset.md)
  : Function to remove an offset

- [`rm_predictors()`](https://iiasa.github.io/ibis.iSDM/reference/rm_predictors.md)
  :

  Remove specific predictors from a
  [distribution](https://rdrr.io/r/stats/Distributions.html) object

- [`rm_priors()`](https://iiasa.github.io/ibis.iSDM/reference/rm_priors.md)
  : Remove existing priors from an existing distribution object

- [`get_data()`](https://iiasa.github.io/ibis.iSDM/reference/get_data.md)
  : Small helper function to obtain predictions from an object

- [`get_ngbvalue()`](https://iiasa.github.io/ibis.iSDM/reference/get_ngbvalue.md)
  : Function to extract nearest neighbour predictor values of provided
  points

- [`get_priors()`](https://iiasa.github.io/ibis.iSDM/reference/get_priors.md)
  : Create priors from an existing distribution model

- [`get_rastervalue()`](https://iiasa.github.io/ibis.iSDM/reference/get_rastervalue.md)
  : Function to extract point values directly from a SpatRaster

## Engines

Statistical models used for estimation of species distributions.

- [`engine_bart()`](https://iiasa.github.io/ibis.iSDM/reference/engine_bart.md)
  : Engine for use of Bayesian Additive Regression Trees (BART)
- [`engine_breg()`](https://iiasa.github.io/ibis.iSDM/reference/engine_breg.md)
  : Engine for Bayesian regularized regression models
- [`engine_gdb()`](https://iiasa.github.io/ibis.iSDM/reference/engine_gdb.md)
  : Use of Gradient Descent Boosting for model estimation
- [`engine_glm()`](https://iiasa.github.io/ibis.iSDM/reference/engine_glm.md)
  : Engine for Generalized linear models (GLM)
- [`engine_glmnet()`](https://iiasa.github.io/ibis.iSDM/reference/engine_glmnet.md)
  : Engine for regularized regression models
- [`engine_inla()`](https://iiasa.github.io/ibis.iSDM/reference/engine_inla.md)
  : Use INLA as engine
- [`engine_inlabru()`](https://iiasa.github.io/ibis.iSDM/reference/engine_inlabru.md)
  : Use inlabru as engine
- [`engine_scampr()`](https://iiasa.github.io/ibis.iSDM/reference/engine_scampr.md)
  : Engine for process models using scampr
- [`engine_stan()`](https://iiasa.github.io/ibis.iSDM/reference/engine_stan.md)
  : Use Stan as engine
- [`engine_xgboost()`](https://iiasa.github.io/ibis.iSDM/reference/engine_xgboost.md)
  : Engine for extreme gradient boosting (XGBoost)

## Create spatial-temporal projections

After a model has been trained, the functions in here can be used to
create projections with scenario() objects. Constraints can be on such
scenarios to limit extrapolations.

- [`scenario()`](https://iiasa.github.io/ibis.iSDM/reference/scenario.md)
  : Create a new scenario based on trained model parameters

- [`project.BiodiversityScenario()`](https://iiasa.github.io/ibis.iSDM/reference/project.md)
  [`project(`*`<BiodiversityScenario>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/project.md)
  [`project.DistributionModel()`](https://iiasa.github.io/ibis.iSDM/reference/project.md)
  [`project(`*`<DistributionModel>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/project.md)
  : Project a fitted model to a new environment and covariates

- [`simulate_population_steps()`](https://iiasa.github.io/ibis.iSDM/reference/simulate_population_steps.md)
  : Simulate population dynamics following the steps approach

- [`add_constraint()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint.md)
  :

  Add a constraint to an existing `scenario`

- [`add_constraint_MigClim()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_MigClim.md)
  : Add constraints to the modelled distribution projection using the
  MigClim approach

- [`add_constraint_adaptability()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_adaptability.md)
  : Adds an adaptability constraint to a scenario object

- [`add_constraint_boundary()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_boundary.md)
  : Adds a boundary or zone constraint to a scenario object

- [`add_constraint_connectivity()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_connectivity.md)
  : Adds a connectivity constraint to a scenario object.

- [`add_constraint_dispersal()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_dispersal.md)
  :

  Add dispersal constraint to an existing `scenario`

- [`add_constraint_minsize()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_minsize.md)
  : Adds a size constraint on a scenario

- [`add_constraint_threshold()`](https://iiasa.github.io/ibis.iSDM/reference/add_constraint_threshold.md)
  : Adds a threshold constraint to a scenario object

## Model summary and validation

Key functions to summarize, validate or extract information from trained
models.

- [`plot(`*`<DistributionModel>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/plot.md)
  [`plot(`*`<BiodiversityDatasetCollection>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/plot.md)
  [`plot(`*`<PredictorDataset>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/plot.md)
  [`plot(`*`<Engine>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/plot.md)
  [`plot(`*`<BiodiversityScenario>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/plot.md)
  : Plot wrappers
- [`bivplot()`](https://iiasa.github.io/ibis.iSDM/reference/bivplot.md)
  : Bivariate prediction plot for distribution objects
- [`nicheplot()`](https://iiasa.github.io/ibis.iSDM/reference/nicheplot.md)
  : Niche plot for distribution objects
- [`print(`*`<distribution>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<BiodiversityDistribution>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<BiodiversityDatasetCollection>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<BiodiversityDataset>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<PredictorDataset>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<DistributionModel>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<BiodiversityScenario>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<Prior>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<PriorList>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<Engine>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<Settings>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<Log>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<Id>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<tbl_df>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  : Print
- [`summary(`*`<distribution>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/summary.md)
  [`summary(`*`<DistributionModel>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/summary.md)
  [`summary(`*`<PredictorDataset>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/summary.md)
  [`summary(`*`<BiodiversityScenario>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/summary.md)
  [`summary(`*`<PriorList>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/summary.md)
  [`summary(`*`<Settings>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/summary.md)
  : Summarises a trained model or predictor object
- [`coef(`*`<DistributionModel>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/coef.md)
  : Obtains the coefficients of a trained model
- [`validate()`](https://iiasa.github.io/ibis.iSDM/reference/validate.md)
  : Validation of a fitted distribution object
- [`similarity()`](https://iiasa.github.io/ibis.iSDM/reference/similarity.md)
  : Calculate environmental similarity of reference datasets to
  predictors.
- [`effects(`*`<DistributionModel>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/effects.md)
  : Plot effects of trained model
- [`partial()`](https://iiasa.github.io/ibis.iSDM/reference/partial.md)
  [`partial.DistributionModel()`](https://iiasa.github.io/ibis.iSDM/reference/partial.md)
  : Obtain partial effects of trained model
- [`spartial()`](https://iiasa.github.io/ibis.iSDM/reference/spartial.md)
  [`spartial.DistributionModel()`](https://iiasa.github.io/ibis.iSDM/reference/spartial.md)
  : Obtain spatial partial effects of trained model
- [`partial_density()`](https://iiasa.github.io/ibis.iSDM/reference/partial_density.md)
  : Visualize the density of the data over the environmental data
- [`limiting()`](https://iiasa.github.io/ibis.iSDM/reference/limiting.md)
  : Identify local limiting factor
- [`threshold()`](https://iiasa.github.io/ibis.iSDM/reference/threshold.md)
  : Threshold a continuous prediction to a categorical layer
- [`ensemble()`](https://iiasa.github.io/ibis.iSDM/reference/ensemble.md)
  : Function to create an ensemble of multiple fitted models
- [`ensemble_partial()`](https://iiasa.github.io/ibis.iSDM/reference/ensemble_partial.md)
  : Function to create an ensemble of partial effects from multiple
  models
- [`ensemble_spartial()`](https://iiasa.github.io/ibis.iSDM/reference/ensemble_spartial.md)
  : Function to create an ensemble of spartial effects from multiple
  models

## Utility functions

These functions are used by engines or spatial processing in the
package. Most of them are for internal use, but can be of use if input
needs to be reformatted.

- [`posterior_predict_stanfit()`](https://iiasa.github.io/ibis.iSDM/reference/posterior_predict_stanfit.md)
  : Create a posterior prediction from a Stan fit object

- [`alignRasters()`](https://iiasa.github.io/ibis.iSDM/reference/alignRasters.md)
  :

  Align a
  [`SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  object to another by harmonizing geometry and extend.

- [`emptyraster()`](https://iiasa.github.io/ibis.iSDM/reference/emptyraster.md)
  :

  Create an empty `SpatRaster` based on a template

- [`get_ngbvalue()`](https://iiasa.github.io/ibis.iSDM/reference/get_ngbvalue.md)
  : Function to extract nearest neighbour predictor values of provided
  points

- [`get_rastervalue()`](https://iiasa.github.io/ibis.iSDM/reference/get_rastervalue.md)
  : Function to extract point values directly from a SpatRaster

- [`predictor_transform()`](https://iiasa.github.io/ibis.iSDM/reference/predictor_transform.md)
  : Spatial adjustment of environmental predictors and raster stacks

- [`predictor_derivate()`](https://iiasa.github.io/ibis.iSDM/reference/predictor_derivate.md)
  : Create spatial derivative of raster stacks

- [`predictor_filter()`](https://iiasa.github.io/ibis.iSDM/reference/predictor_filter.md)
  : Filter a set of correlated predictors to fewer ones

- [`predictor_summarize_zones()`](https://iiasa.github.io/ibis.iSDM/reference/predictor_summarize_zones.md)
  : Summarize and replace values in predictors with an aggregation

- [`interpolate_gaps()`](https://iiasa.github.io/ibis.iSDM/reference/interpolate_gaps.md)
  : Approximate missing time steps between dates

- [`run_stan()`](https://iiasa.github.io/ibis.iSDM/reference/run_stan.md)
  : Fit a cmdstanr model

- [`sanitize_names()`](https://iiasa.github.io/ibis.iSDM/reference/sanitize_names.md)
  : Sanitize variable names

- [`get_data()`](https://iiasa.github.io/ibis.iSDM/reference/get_data.md)
  : Small helper function to obtain predictions from an object

- [`combine_formulas()`](https://iiasa.github.io/ibis.iSDM/reference/combine_formulas.md)
  : Combine or concatenate multiple formula objects

- [`stancode()`](https://iiasa.github.io/ibis.iSDM/reference/stancode.md)
  [`stancode.DistributionModel()`](https://iiasa.github.io/ibis.iSDM/reference/stancode.md)
  : Show the stan code from a trained model

- [`write_model()`](https://iiasa.github.io/ibis.iSDM/reference/write_model.md)
  : Save a model for later use

- [`write_output()`](https://iiasa.github.io/ibis.iSDM/reference/write_output.md)
  : Generic function to write spatial outputs

- [`write_summary()`](https://iiasa.github.io/ibis.iSDM/reference/write_summary.md)
  : Generic function to write summary outputs from created models.

- [`load_model()`](https://iiasa.github.io/ibis.iSDM/reference/load_model.md)
  : Load a pre-computed model

- [`mask.DistributionModel()`](https://iiasa.github.io/ibis.iSDM/reference/mask.md)
  [`mask.BiodiversityDatasetCollection()`](https://iiasa.github.io/ibis.iSDM/reference/mask.md)
  [`mask.PredictorDataset()`](https://iiasa.github.io/ibis.iSDM/reference/mask.md)
  [`mask.BiodiversityScenario()`](https://iiasa.github.io/ibis.iSDM/reference/mask.md)
  : Mask data with an external layer

- [`myLog()`](https://iiasa.github.io/ibis.iSDM/reference/myLog.md) :
  Custom messaging function for scripts

- [`predictor_homogenize_na()`](https://iiasa.github.io/ibis.iSDM/reference/predictor_homogenize_na.md)
  : Homogenize NA values across a set of predictors.

- [`run_parallel()`](https://iiasa.github.io/ibis.iSDM/reference/run_parallel.md)
  : Parallel computation of function

- [`thin_observations()`](https://iiasa.github.io/ibis.iSDM/reference/thin_observations.md)
  : Functionality for geographic and environmental thinning

- [`unwrap_model()`](https://iiasa.github.io/ibis.iSDM/reference/unwrap_model.md)
  : Unwrap a model for later use

- [`wrap_model()`](https://iiasa.github.io/ibis.iSDM/reference/wrap_model.md)
  : Wrap a model for later use

## Class definitions and methods

These pages document the package’s internal data structures and
functions for manipulating them—they contain information that is really
only useful when adding new functionality to the package.

- [`BiodiversityDataset-class`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDataset-class.md)
  [`BiodiversityDataset`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDataset-class.md)
  : BiodiversityDataset prototype description
- [`BiodiversityDatasetCollection-class`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDatasetCollection-class.md)
  [`BiodiversityDatasetCollection`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDatasetCollection-class.md)
  : BiodiversityDatasetCollection super class description
- [`BiodiversityDistribution-class`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md)
  [`BiodiversityDistribution`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityDistribution-class.md)
  : Biodiversity Distribution master class
- [`BiodiversityScenario-class`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
  [`BiodiversityScenario`](https://iiasa.github.io/ibis.iSDM/reference/BiodiversityScenario-class.md)
  : Class for a biodiversity scenario from a trained model
- [`DistributionModel-class`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
  [`DistributionModel`](https://iiasa.github.io/ibis.iSDM/reference/DistributionModel-class.md)
  : Class for the trained Model object
- [`Engine-class`](https://iiasa.github.io/ibis.iSDM/reference/Engine-class.md)
  [`Engine`](https://iiasa.github.io/ibis.iSDM/reference/Engine-class.md)
  : Engine class description
- [`Prior-class`](https://iiasa.github.io/ibis.iSDM/reference/Prior-class.md)
  [`Prior`](https://iiasa.github.io/ibis.iSDM/reference/Prior-class.md)
  : Base Prior class
- [`PriorList-class`](https://iiasa.github.io/ibis.iSDM/reference/PriorList-class.md)
  [`PriorList`](https://iiasa.github.io/ibis.iSDM/reference/PriorList-class.md)
  : List of Priors supplied to an class
- [`Settings-class`](https://iiasa.github.io/ibis.iSDM/reference/Settings-class.md)
  [`Settings`](https://iiasa.github.io/ibis.iSDM/reference/Settings-class.md)
  : Prototype for model settings object
- [`PredictorDataset-class`](https://iiasa.github.io/ibis.iSDM/reference/PredictorDataset-class.md)
  [`PredictorDataset`](https://iiasa.github.io/ibis.iSDM/reference/PredictorDataset-class.md)
  : PredictorDataset class description
- [`Log-class`](https://iiasa.github.io/ibis.iSDM/reference/Log-class.md)
  [`Log`](https://iiasa.github.io/ibis.iSDM/reference/Log-class.md) :
  Log prototype.

## Parallel processing

Functions to enable or set up parallel processing

- [`ibis_future()`](https://iiasa.github.io/ibis.iSDM/reference/ibis_future.md)
  : Internal function to enable (a)synchronous parallel processing
- [`ibis_enable_parallel()`](https://iiasa.github.io/ibis.iSDM/reference/ibis_enable_parallel.md)
  : Set the parallel processing flag to TRUE
- [`ibis_set_strategy()`](https://iiasa.github.io/ibis.iSDM/reference/ibis_set_strategy.md)
  : Set the number of threads for parallel processing.
- [`ibis_set_threads()`](https://iiasa.github.io/ibis.iSDM/reference/ibis_set_threads.md)
  : Set the threads for parallel processing.
- [`run_parallel()`](https://iiasa.github.io/ibis.iSDM/reference/run_parallel.md)
  : Parallel computation of function

## Miscellaneous functions

Other functions only relevant for development

- [`as.Id()`](https://iiasa.github.io/ibis.iSDM/reference/as.Id.md) : As
  Id
- [`is.Id()`](https://iiasa.github.io/ibis.iSDM/reference/is.Id.md) :
  Check whether a provided object is truly of a specific type
- [`check()`](https://iiasa.github.io/ibis.iSDM/reference/check.md) :
  Check objects in the package for common errors or issues
- [`bivplot()`](https://iiasa.github.io/ibis.iSDM/reference/bivplot.md)
  : Bivariate prediction plot for distribution objects
- [`ibis_dependencies()`](https://iiasa.github.io/ibis.iSDM/reference/ibis_dependencies.md)
  : Install ibis dependencies
- [`ibis_enable_parallel()`](https://iiasa.github.io/ibis.iSDM/reference/ibis_enable_parallel.md)
  : Set the parallel processing flag to TRUE
- [`ibis_future()`](https://iiasa.github.io/ibis.iSDM/reference/ibis_future.md)
  : Internal function to enable (a)synchronous parallel processing
- [`ibis_options()`](https://iiasa.github.io/ibis.iSDM/reference/ibis_options.md)
  : Print ibis options
- [`ibis_set_strategy()`](https://iiasa.github.io/ibis.iSDM/reference/ibis_set_strategy.md)
  : Set the number of threads for parallel processing.
- [`ibis_set_threads()`](https://iiasa.github.io/ibis.iSDM/reference/ibis_set_threads.md)
  : Set the threads for parallel processing.
- [`is.Raster()`](https://iiasa.github.io/ibis.iSDM/reference/is.Raster.md)
  : Tests if an input is a SpatRaster object.
- [`is.Waiver()`](https://iiasa.github.io/ibis.iSDM/reference/is.Waiver.md)
  : Is the provided object of type waiver?
- [`is.formula()`](https://iiasa.github.io/ibis.iSDM/reference/is.formula.md)
  : Check whether a formula is valid
- [`is.stars()`](https://iiasa.github.io/ibis.iSDM/reference/is.stars.md)
  : Tests if an input is a stars object.
- [`modal()`](https://iiasa.github.io/ibis.iSDM/reference/modal.md) :
  Calculate the mode of a provided vector
- [`new_id()`](https://iiasa.github.io/ibis.iSDM/reference/new_id.md) :
  Identifier
- [`new_waiver()`](https://iiasa.github.io/ibis.iSDM/reference/new_waiver.md)
  : Waiver
- [`nicheplot()`](https://iiasa.github.io/ibis.iSDM/reference/nicheplot.md)
  : Niche plot for distribution objects
- [`plot(`*`<DistributionModel>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/plot.md)
  [`plot(`*`<BiodiversityDatasetCollection>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/plot.md)
  [`plot(`*`<PredictorDataset>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/plot.md)
  [`plot(`*`<Engine>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/plot.md)
  [`plot(`*`<BiodiversityScenario>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/plot.md)
  : Plot wrappers
- [`print(`*`<distribution>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<BiodiversityDistribution>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<BiodiversityDatasetCollection>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<BiodiversityDataset>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<PredictorDataset>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<DistributionModel>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<BiodiversityScenario>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<Prior>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<PriorList>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<Engine>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<Settings>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<Log>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<Id>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  [`print(`*`<tbl_df>`*`)`](https://iiasa.github.io/ibis.iSDM/reference/print.md)
  : Print
- [`render_html()`](https://iiasa.github.io/ibis.iSDM/reference/render_html.md)
  : render_html
- [`run_stan()`](https://iiasa.github.io/ibis.iSDM/reference/run_stan.md)
  : Fit a cmdstanr model
- [`myLog()`](https://iiasa.github.io/ibis.iSDM/reference/myLog.md) :
  Custom messaging function for scripts
