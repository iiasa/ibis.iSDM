// Data Flags
int<lower=1> N;                    // total number of observations
int observed[N];                   // response variable
int<lower=1> K;                    // number of population-level effects
matrix[N, K] X;                    // population-level design matrix
vector[N] offsets;                 // Vector with exposure offset and other eventual offsets
// Generic parameters for modelling options
int<lower=0,upper=1> has_intercept;// 1 = yes
int<lower=0,upper=2> has_spatial;  // 0 = none, 1 = CAR, 2 = GRMF
