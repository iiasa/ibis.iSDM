// Fit a Cox process in Stan

data{
  int<lower = 1> n;
  vector[n] x;
  int<lower = 0> y[n];
}
parameters{
  real beta0;
  real beta1;
  real<lower = 0, upper = 5> sigma_noise;
  vector[n] noise;
}
transformed parameters{
}
model{
  //priors
  target += normal_lpdf(beta0 | 0,5);
  target += normal_lpdf(beta1 | 0,10);
  target += uniform_lpdf(sigma_noise | 0,1);

  // Prior for the noise
  target += normal_lpdf(noise | 0, sigma_noise);

  // likelihood
  target += poisson_log_lpmf(y | beta0 + beta1 * x + noise);
}
generated quantities{
 vector[n] lambda_rep;
 lambda_rep = exp(beta0 + beta1 * x + noise);
}
