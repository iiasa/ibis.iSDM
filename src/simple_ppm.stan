// ----------------------------- //
// Stan code for a simple Point Process Model based on a single data source
functions {
  //#include src/spatial_functions.stan // Stan functions such as CAR and GP models
}
data {
  // Data Flags
  int<lower=1> N;                    // total number of observations
  int observed[N];                   // response variable
  int<lower=1> K;                    // number of population-level effects
  matrix[N, K] X;                    // population-level design matrix
  vector[N] offset_exposure;         // Vector with exposure offset
  // Generic parameters for modelling options
  int<lower=0,upper=1> has_intercept;// 1 = yes
  int<lower=0,upper=2> has_spatial;  // 0 = none, 1 = CAR, 2 = GRMF
  // // Prediction flags
  // int<lower=0,upper=1> do_prediction;// 1 = yes
  // int<lower=1> Npred;                // Number of prediction observations
  // matrix[Npred, K] Xpred;            // Prediction matrix
  // vector[Npred] offset_pred_exposure;// Prediction offset
}
transformed data {
}
parameters {
  vector[K] b;                       // population-level effects for slopes
  // real<lower = 0, upper = 5> sigma_noise;  // Sigma noise effect
  // vector[N] noise;                   //  Generic error noise
}
transformed parameters {
}
model {
  vector[N] lambda;
  vector[N] mu;
  // Initialize linear predictor term with or without intercept
  if(has_intercept == 1){
    mu = offset_exposure;
  } else {
    mu = rep_vector(0.0, N) + offset_exposure;
  }
  // Build lambda
  // lambda = mu + X * b ;
  for (n in 1:N){
    lambda[n] = (mu[n] + X[n] * b);
  }

  // Add priors for all covariates
  for (j in 1:K){
      target += normal_lpdf(b[j] | 0, 5);
  }
  // ---------------------- //

  // likelihood
  // observed_i ~ poisson( exp(lambda_i) );
  for (i in 1:N) {
    target += poisson_log_lpmf(observed[i] | lambda[i]);
    // observed[i] ~ poisson_log(lambda[i]);
  }

  // log poisson probability mass of y given log-rate alpha+x*beta
  // target += poisson_log_glm_lpmf(observed | X, mu, b);
}
generated quantities {
  // Simulated prediction. This works only if parameters are not local, transformed and thus part of the output
  // vector[N] observedSim;
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
}
