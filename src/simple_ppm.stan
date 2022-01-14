// ----------------------------- //
// Stan code for a simple Point Process Model based on a single data source
functions {
  //#include src/spatial_functions.stan // Stan functions such as CAR and GP models
  #include src/prior_functions.stan
}
data {
  #include src/data_parameters.stan
  // data for the horseshoe prior
  real<lower=0> hs_df;  // local degrees of freedom
  real<lower=0> hs_df_global;  // global degrees of freedom
  real<lower=0> hs_df_slab;  // slab degrees of freedom
  real<lower=0> hs_scale_global;  // global prior scale
  real<lower=0> hs_scale_slab;  // slab prior scale
}
transformed data {
  // int Kc = K - 1;
  // matrix[N, Kc] Xc;    // centered version of X without an intercept
  // vector[Kc] means_X;  // column means of X before centering
  // for (i in 2:K) {
  //   means_X[i - 1] = mean(X[, i]);
  //   Xc[, i - 1] = X[, i] - means_X[i - 1];
  // }
}
parameters {
  // Intercept
  // if (has_intercept) {
    // real Intercept;  // temporary intercept for centered predictors
    // local parameters for horseshoe prior
    // vector[Kc] zb;
    // vector<lower=0>[Kc] hs_local;
  // } else {
    vector[K] zb;
    vector<lower=0>[K] hs_local;
  // }
  // horseshoe shrinkage parameters
  real<lower=0> hs_global;  // global shrinkage parameters
  real<lower=0> hs_slab;  // slab regularization parameter
}
transformed parameters {
  // if (has_intercept) {
    // vector[Kc] beta;  // population-level effects
  // } else {
    vector[K] beta;  // population-level effects
  // }
  // compute actual regression coefficients
  beta = horseshoe(zb, hs_local, hs_global, hs_scale_slab^2 * hs_slab);
}

model {
  // ---------------------- //
  // likelihood
  #include src/poipo_ll_poisson.stan

  // Add priors for all covariates
  // for (j in 1:K){
  //     target += normal_lpdf(b[j] | 0, 5);
  // }

  // priors including constants
  target += std_normal_lpdf(zb);
  target += student_t_lpdf(hs_local | hs_df, 0, 1)
    - rows(hs_local) * log(0.5);

  // if (has_intercept) {
  //   target += student_t_lpdf(Intercept | 3, -1.25795845831865, 2.5);
  // }
  target += student_t_lpdf(hs_global | hs_df_global, 0, hs_scale_global)
    - 1 * log(0.5);
  target += inv_gamma_lpdf(hs_slab | 0.5 * hs_df_slab, 0.5 * hs_df_slab);
}
generated quantities {
  // // Simulated prediction. This works only if parameters are not local, transformed and thus part of outcome
  // vector[N] observedSim;             // Vector to contain the observed Preditions
  // vector[N] log_lik; // Vector to contain the log-likelihood
  // for (i in 1:N) {
  //     if (lambda[i] > 20) {
  //     // To avoid errors during warumup
  //       observedSim[i] = poisson_log_rng(20);
  //     } else {
  //       observedSim[i] = poisson_log_rng(lambda[i]);
  //     }
  //   // Save log-likelihood
  //   log_lik[i] = poisson_log_lpmf(observed[i] | lambda[i]);
  // }
  // if (has_intercept){
  //   // actual population-level intercept
  //   real beta_Intercept = Intercept - dot_product(means_X, beta);
  // }
}
