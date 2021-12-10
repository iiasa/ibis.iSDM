// ----------------------------- //
// Stan code for a simple Point Process Model based on a single data source

// ----------------------------- //
functions {
  #include src/stan_functions.stan
}
// ----------------------------- //
data{
  // dimensions
  int<lower=0> N;  // number of observations
  int<lower=0> K;  // number of predictors
  int<lower=0> y[N]; // The observed outcome variable
  // Predictors
  vector[N] x;
  // matrix[N, K] X;   // Matrix with predictors, each with a vector[N] x;
}
parameters{
  vector[N] intercept;
  real<lower=-10, upper = 10> beta1;
  real<lower = 0, upper = 5> sigma_noise;
  vector[N] noise;
}
transformed parameters {

}
model{
  //priors
  target += normal_lpdf(intercept | 0,5);
  target += normal_lpdf(beta1 | 0,10);
  target += uniform_lpdf(sigma_noise | 0,1);

  // Prior for the noise
  target += normal_lpdf(noise | 0, sigma_noise);

  // likelihood
  target += poisson_log_lpmf(y | intercept + beta1 * x + noise);
}
generated quantities{
  vector[N] lambda_rep;
  lambda_rep = exp(intercept + beta1 * x + noise);
}
