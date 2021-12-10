// Generic functions to use

# Exponential covariance structure through a Gaussian process
matrix GP(matrix x, real sigma_sq, real scale, real delta) {
    int N = dims(x)[1];
    matrix[N, N] K;
    for (i in 1:(N-1)) {
      K[i, i] = sigma_sq + delta;
      for (j in (i + 1):N) {
        K[i, j] = sigma_sq * exp(- x[i,j] / scale );
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sigma_sq + delta;
    return K;
}
