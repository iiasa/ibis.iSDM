// Data Flags
int<lower=1> N;                    // total number of observations
int observed[N];                   // response variable
int<lower=1> K;                    // number of population-level effects
matrix[N, K] X;                    // population-level design matrix
vector[N] offsets;                 // Vector with exposure offset and other eventual offsets
// Generic parameters for modelling options
int<lower=0,upper=1> has_intercept;// 1 = yes
int<lower=0,upper=2> has_spatial;  // 0 = none, 1 = CAR, 2 = GRMF
// prior family: 0 = none, 1 = normal, 2 = student_t, 3 = hs, 4 = hs_plus,
//   5 = laplace, 6 = lasso, 7 = product_normal
// int<lower=0,upper=7> prior_dist;
// int<lower=0,upper=2> prior_dist_for_intercept;
// // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = exponential
// int<lower=0,upper=3> prior_dist_for_aux;
// // prior family: 0 = none, 1 = normal, 2 = student_t, 3 = exponential
// int<lower=0,upper=3> prior_dist_for_smooth;
// // Prediction flags
// int<lower=0,upper=1> do_prediction;// 1 = yes
// int<lower=1> Npred;                // Number of prediction observations
// matrix[Npred, K] Xpred;            // Prediction matrix
// vector[Npred] offset_pred_exposure;// Prediction offset
