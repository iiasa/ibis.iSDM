// vecto<r[N] lambda;
vector[N] mu;

// Initialize linear predictor term with intercept
mu = Intercept + rep_vector(0.0, N) + offsets;
// log poisson probability mass of y given log-rate alpha+x*beta
target += poisson_log_glm_lpmf(observed | Xc, mu, beta);

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
