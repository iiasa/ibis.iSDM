data{
  int nsite; // Number of sites
  int npixel;
  int npar_det;
  int npar_state;
  int npar_thin;
  int m; // total number of data points
  int CONSTANT;//
  int W;
  matrix[npixel, npar_det] v; // detection covariates for the entire region B
  matrix[npixel, npar_state] x_s; // covariates for latent state
  matrix[npixel, npar_thin] h_s; //covariates for thinning probability
  real cell_area;//log area of grid cell
  int pa_pixel[nsite]; //The pixels the presence only data occurred
  int po_pixel[m]; //The pixels the presence only data occurred
  int ones[m];
  int y[nsite];
}

transformed data{
  int obs_pa_po[npixel];
  for(prr in 1:npixel){
    int obs = 0;
    for(po in 1:m){
      if(po_pixel[po] == prr ){
          obs =1;
          }
      }
      if(obs == 0){
        for(pa in 1:nsite){
          if(pa_pixel[pa]==prr){
          if (y[pa] > 0) {
            obs =1;
            }
          }
        }
        }
      obs_pa_po[prr] = obs;
  }
}


parameters{
  vector[npar_det] a;  //det/non-det data model regression coefficients
  vector[npar_state] beta; //latent state model regression coefficients
  vector[npar_thin] cc; //presence-only data model regression coefficients
}

transformed parameters{
  vector<lower=0>[npixel] lambda;
  vector<lower=0,upper=1>[npixel] psi; // probability of Species presence in a gridcell
  vector[npixel] b; //  presence only thinning prob linear predictor
  real po_denominator; //The presence_only data model denominator
  vector[nsite] rho; // detection/non-detection data model linear predictor
  vector[m] po_p;
  vector[nsite] log_lik_occ;
  lambda = exp(x_s * beta + cell_area);
  psi = 1 - exp(-lambda);
  b = inv_logit(h_s * cc);
  po_denominator = dot_product(lambda , b) / m;
  po_p = exp( log(lambda[po_pixel] .* b[po_pixel] ) -
          log(po_denominator)
          ) / CONSTANT ;
  rho = inv_logit(v[pa_pixel, ] * a);
  for(site in 1:nsite){
      // https://mbjoseph.github.io/posts/2020-04-28-a-step-by-step-guide-to-marginalizing-over-discrete-parameters-for-ecologists-using-stan/
  if (y[site] > 0) {
      log_lik_occ[site] = log(psi[pa_pixel[site]]) + binomial_lpmf(y[site] | W, rho[site]);
    } else {
      log_lik_occ[site] = log_sum_exp(
        log(psi[pa_pixel[site]]) + binomial_lpmf(y[site] | W, rho[site]),
        log1m(psi[pa_pixel[site]])
      );
    }
}
}


model{
  // # Priors for latent state model
  to_vector(beta) ~ normal(0, 1);
  // # Priors for presence-only data model
  to_vector(cc) ~ logistic(0, 1);
  // # Priors for det/non-det data model
  to_vector(a) ~ logistic(0, 1);

  // # Loop through each presence-only datapoint
  // #  using Bernoulli one's trick. The numerator
  // #  is just the thinned poisson process for
  // #  the ith data point.
  target += bernoulli_lpmf(ones | po_p);
  target += log_lik_occ;
}


generated quantities{
  // # Derived parameter, the number of cells occupied
  vector[npixel] lp_present; // [z=1][y=0 | z=1] / [y=0] on a log scale
//   One key tweak here: setting z=1 when known, and drawing from the posterior of z if unknown:
// [z | data, params] = [data | z, params] * [z | params] / \sum_{z=0}^1  [data | z, params] * [z | params]
  // https://github.com/mbjoseph/scr-stan/blob/master/ch07/individual-heterogeneity-ranefs.stan#L84-L91
  real zsum;
  int z[npixel];
  for (p in 1:npixel) {
      if(obs_pa_po[p]){
        z[p] = 1;
        } else{
          lp_present[p] = log(psi[p]) -
          log_sum_exp(log(psi[p]) ,
                  bernoulli_lpmf(0| psi[p] ) );
                // log_sum_exp(psi[p],
                // bernoulli_lpmf(0|  ))
        z[p] = bernoulli_rng(exp(lp_present[p]));
      }
  }

  zsum = sum(z);

}


