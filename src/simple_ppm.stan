#include stan_functions.stan
// ----------------------------- //
// Stan code for a simple Point Process Model based on a single data source

// ----------------------------- //
functions {

}
// ----------------------------- //
data{
  int<lower = 1> n; // Number of data points
  vector[n] x;
  int<lower = 0> observed[n]; // The observed variable
}
parameters{
  real beta0;
  real beta1;
  real<lower = 0, upper = 5> sigma_noise;
}
transformed parameters{

}
model{
  //priors
  target += normal_lpdf(beta0 | 0,5);
  target += normal_lpdf(beta1 | 0,10);
  target += uniform_lpdf(sigma_noise | 0,1);

  // likelihood
  target += poisson_log_lpmf(observed | beta0 + beta1 * x );
}
generated quantities{
  vector[n] lambda_rep;
  lambda_rep = exp(beta0 + beta1 * x );
}
