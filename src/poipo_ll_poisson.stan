// vecto<r[N] lambda;
vector[N] mu;

// Initialize linear predictor term with or without intercept
if(has_intercept == 1){
  mu = offset_exposure;
} else {
  mu = rep_vector(0.0, N) + offset_exposure;
}
// Build lambda as lambda = mu + X * b ;
// for (n in 1:N){
//   lambda[n] = (mu[n] + X[n] * b);
// }
// ---------------------- //

// likelihood
// observed_i ~ poisson( exp(lambda_i) );
// for (i in 1:N) {
//   target += poisson_log_lpmf(observed[i] | lambda[i]);
//   // observed[i] ~ poisson_log(lambda[i]);
// }

// Alternative more efficient formulation:
# log poisson probability mass of y given log-rate alpha+x*beta
target += poisson_log_glm_lpmf(observed | X, mu, beta);
