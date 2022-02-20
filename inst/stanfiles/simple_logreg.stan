functions {
}
data {
  #include src/data_parameters.stan
}
transformed data {
}
parameters {
  vector[K] b;  // population-level effects
}
transformed parameters {
}
model {
  // likelihood including constants
  #include src/poipa_ll_bernoulli.stan
}
generated quantities {
}
