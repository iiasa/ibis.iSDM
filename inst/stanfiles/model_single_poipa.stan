// Single presence-absence model.
//
// Uses a Bernoulli likelihood. The default link is logit, while cloglog is
// available through the dataset link metadata.
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
    if (link[1] == 2) {
      target += weights[n] * bernoulli_logit_lpmf(observed[n] | eta[n]);
    } else {
      target += weights[n] * bernoulli_lpmf(observed[n] | inv_cloglog_response(eta[n]));
    }
  }
}
generated quantities {
}
