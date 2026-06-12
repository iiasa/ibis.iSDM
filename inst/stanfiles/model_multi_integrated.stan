// Integrated multi-dataset species distribution model.
//
// All datasets share the same covariate effects (`beta`). Dataset/source effects
// are represented through the intercept selected by `intercept_id`.
functions {
  {{other_functions}}
  {{functions_extra}}
}
data {
  {{data_parameters}}
  {{data_extra}}
}
transformed data {
}
parameters {
  {{parameters}}
}
transformed parameters {
  {{transformed_parameters}}
}
model {
  vector[N] eta;

  Intercept ~ student_t(3, 0, 2.5);
  {{coefficient_priors}}

  eta = X * beta + offsets;
  for (n in 1:N) {
    eta[n] += Intercept[intercept_id[n]];
  }

  for (n in 1:N) {
    if (likelihood[dataset[n]] == 1) {
      target += weights[n] * poisson_log_lpmf(observed[n] | eta[n]);
    } else {
      if (link[dataset[n]] == 2) {
        target += weights[n] * bernoulli_logit_lpmf(observed[n] | eta[n]);
      } else {
        target += weights[n] * bernoulli_lpmf(observed[n] | inv_cloglog_response(eta[n]));
      }
    }
  }
}
generated quantities {
}
