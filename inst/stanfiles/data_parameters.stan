// Common data declarations for Stan species distribution models.
//
// The same data list is used for single-dataset and integrated models so the R
// engine can switch templates without changing the data-building path.
int<lower=1> N;                              // Total number of observations.
int<lower=1> K;                              // Number of shared covariate effects.
matrix[N, K] X;                              // Shared covariate design matrix.
array[N] int<lower=0> observed;              // Count or binary response.
vector[N] offsets;                           // Additive link-scale offsets.
vector<lower=0>[N] weights;                  // Likelihood weights.

int<lower=1> J;                              // Number of biodiversity datasets.
int<lower=1> J_intercept;                    // Number of intercept parameters.
array[N] int<lower=1, upper=J> dataset;      // Dataset id per observation.
array[J] int<lower=1, upper=2> likelihood;   // 1 = Poisson, 2 = Bernoulli.
array[J] int<lower=1, upper=3> link;         // 1 = log, 2 = logit, 3 = cloglog.
array[N] int<lower=1, upper=J_intercept> intercept_id;
