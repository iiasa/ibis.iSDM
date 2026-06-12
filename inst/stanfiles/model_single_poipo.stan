// Single presence-only point process model.
//
// Uses a Poisson likelihood with a log-intensity linear predictor. Exposure
// and spatial offsets are supplied on the log/link scale via `offsets`.
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

  eta = X * beta + offsets + Intercept[1];
  for (n in 1:N) {
    target += weights[n] * poisson_log_lpmf(observed[n] | eta[n]);
  }
}
generated quantities {
}
