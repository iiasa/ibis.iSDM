// Code from spatialfusion

data {
           int n_point;
           int n_area;
           int n_grid;
           int n_neighbor;
           int n_sample;
           int n_w;
           int n_point_var;
           int n_area_var;
           int n_pp_var;
           int n_ban;
           int n_norm; // number of normally distributed responses
           int idx_norm;vector[n_point] Y_point1;int p_point; // number of coefficient for point
              matrix[n_point, p_point] X_point; // design matrix for point
           int nearind[n_point + n_grid - 1, n_neighbor];
           vector[n_neighbor] sC_site_nei[n_point + n_grid - 1];
           matrix[n_neighbor, n_neighbor] sC_nei[n_point + n_grid - 1];
           int Y_area1[n_area];
              matrix[n_area, n_sample] A1; // aggregation matrix for areal
              int p_area; // number of coefficient for area
              matrix[n_area, p_area] X_area; // design matrix for area
             matrix[n_neighbor, n_neighbor] C_nei[n_sample];
             vector[n_neighbor] C_site_nei[n_sample];
             int nearind_sample[n_sample, n_neighbor];int Y_pp1[n_grid];
              real offset[n_grid, n_pp_var];
              real area;
           int Z_pos;
  }

           parameters{
           vector[p_point] beta_p[n_point_var]; // coefficients
                    vector[p_area] beta_a[n_area_var];
                    row_vector<lower = 0>[n_w] Z_ban1;
 row_vector[n_w] Z_ban2;
 row_vector[n_w] Z_ban3;
           real<lower = 0> tau_sq[n_norm];
           positive_ordered[n_w] phi;
           matrix[n_w, n_point + n_grid + n_sample] noise;
           }

           transformed parameters{
           matrix[n_w, n_point + n_grid] w;
           matrix[n_w, n_point + n_grid] w_var;

           matrix[n_w, n_sample] wa;
           matrix[n_w, n_sample] wa_var;
           vector[n_area] wA[n_area_var];
           vector[n_neighbor] C_site_nei_C_nei_inv;
             vector[n_neighbor] C_site_nei_phi;
             matrix[n_neighbor,n_neighbor] C_nei_phi;
             row_vector[n_w] Z_1;
 row_vector[n_w] Z_2;
 row_vector[n_w] Z_3;
           Z_1 = Z_ban1;
 Z_2 = Z_ban2;
 Z_3 = Z_ban3;
           for (l in 1:n_w){
           w_var[l, 1] = 1;
           w[l, 1] = sqrt(w_var[l, 1]) * noise[l, 1];
           // for transforming w
           for (i in 2:(n_point + n_grid)) {
           int dim;
           matrix[ i < (n_neighbor + 1)? (i - 1) : n_neighbor, i < (n_neighbor + 1)? (i - 1): n_neighbor] sC_nei_phi;
           vector[ i < (n_neighbor + 1)? (i - 1) : n_neighbor] sC_site_nei_phi;
           vector[ i < (n_neighbor + 1)? (i - 1) : n_neighbor] sC_site_nei_C_nei_inv;
           dim = (i < (n_neighbor + 1))? (i-1) : n_neighbor;
           if(dim == 1){sC_nei_phi[1, 1] = 1;}
           else{
           for (j in 1:dim){
           for (k in j:dim){
           sC_nei_phi[j, k] = exp(- sC_nei[(i - 1)][j,k] / phi[l]);
           sC_nei_phi[k, j] = sC_nei_phi[j, k];
           }}}
           sC_site_nei_phi = exp(- sC_site_nei[(i - 1)][1:dim] / phi[l]);
           sC_site_nei_C_nei_inv = mdivide_left_spd(sC_nei_phi, sC_site_nei_phi);// m by m times m by n
           w_var[l, i] = (1 - dot_product(sC_site_nei_C_nei_inv, sC_site_nei_phi));
           w[l, i] = dot_product(sC_site_nei_C_nei_inv, w[l, nearind[i-1, 1:dim]]) + sqrt(w_var[l, i]) * noise[l, i]; // 1 by m, m by 1
           }

           // for transforming wa
           for (i in 1:n_sample) { // for each predicted location
           for (j in 1:n_neighbor){
           C_site_nei_phi[j] = exp(- C_site_nei[i][j]/phi[l]);
           }
           for (j in 1:n_neighbor){
           for (k in j:n_neighbor){
           C_nei_phi[j,k] = exp(- C_nei[i][j,k]/phi[l]);
           C_nei_phi[k,j] = C_nei_phi[j,k];
           }}
           C_site_nei_C_nei_inv = mdivide_left_spd(C_nei_phi, C_site_nei_phi);// m by m times m by n
           wa_var[l, i] = (1 - dot_product(C_site_nei_C_nei_inv, C_site_nei_phi));
           wa[l, i] = dot_product(C_site_nei_C_nei_inv, append_row(to_vector(w[l, ]), to_vector(wa[l, ]))[nearind_sample[i,]]) + sqrt(wa_var[l, i]) * noise[l, n_point + n_grid + i]; // 1 by m, m by 1

           }}
           wA[1] = log(A1 * to_vector(exp(Z_2 * wa))); // log link
}

           model{
           tau_sq ~ inv_gamma(2,1);
           phi ~ normal(1,10);
           for (i in 1:n_w){
           noise[i,] ~ normal(0, 1);
           }
  Z_1 ~ normal(1,1);
Z_2 ~ normal(1,1);
Z_3 ~ normal(1,1);
beta_p[1] ~ normal(0,10);
                      Y_point1 ~ normal(X_point * beta_p[1] + to_vector(Z_1 * w[,1:n_point]), sqrt(tau_sq[1]));
beta_a[1] ~ normal(0,10);
                      Y_area1 ~ poisson_log(X_area * beta_a[1] + wA[1]);
Y_pp1 ~ poisson_log(log(area) + to_vector(log(offset[,1])) + to_vector(Z_3 * w[,(n_point+1):(n_point + n_grid)]));
}
