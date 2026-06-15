#include <Rcpp.h>
using namespace Rcpp;

//' Two-dimensional kernel density estimation
//'
//' Compute a bivariate Gaussian kernel density estimate.
//'
//' @param x Numeric vector of x coordinates.
//' @param y Numeric vector of y coordinates.
//' @param h Numeric vector of length 2 giving the bandwidths.
//' @param n Integer vector of length 2 giving the grid size.
//' @param lims Numeric vector of length 4 giving the limits
//'   \code{c(xmin, xmax, ymin, ymax)}.
//'
//' @return A list with components \code{x}, \code{y}, and \code{z}.
//'
//' @noRd
// [[Rcpp::export]]
List cpp_kde2d(NumericVector x, NumericVector y, NumericVector h,
               IntegerVector n, NumericVector lims) {

  int ndata = x.size();
  int n_gx  = n[0];
  int n_gy  = n[1];

  // Build the evaluation grid
  NumericVector gx(n_gx);
  NumericVector gy(n_gy);

  double step_x = (n_gx > 1) ? (lims[1] - lims[0]) / (n_gx - 1) : 0.0;
  double step_y = (n_gy > 1) ? (lims[3] - lims[2]) / (n_gy - 1) : 0.0;

  for (int i = 0; i < n_gx; i++) {
    gx[i] = lims[0] + i * step_x;
  }
  for (int i = 0; i < n_gy; i++) {
    gy[i] = lims[2] + i * step_y;
  }

  // Divide bandwidth by 4 for S compatibility (legacy)
  double hx = h[0] / 4.0;
  double hy = h[1] / 4.0;

  // Pre-compute dnorm values for x-dimension: dx[i * ndata + k]
  std::vector<double> dx(n_gx * ndata);
  for (int i = 0; i < n_gx; i++) {
    for (int k = 0; k < ndata; k++) {
      double val = (gx[i] - x[k]) / hx;
      dx[i * ndata + k] = R::dnorm(val, 0.0, 1.0, 0);
    }
  }

  // Pre-compute dnorm values for y-dimension: dy[j * ndata + k]
  std::vector<double> dy(n_gy * ndata);
  for (int j = 0; j < n_gy; j++) {
    for (int k = 0; k < ndata; k++) {
      double val = (gy[j] - y[k]) / hy;
      dy[j * ndata + k] = R::dnorm(val, 0.0, 1.0, 0);
    }
  }

  // z[i,j] = sum_k dx[i,k] * dy[j,k]  /  (ndata * hx * hy)
  // This is equivalent to tcrossprod(ax, ay) / (ndata * hx * hy)
  double denom = (double)ndata * hx * hy;
  NumericMatrix z(n_gx, n_gy);

  for (int i = 0; i < n_gx; i++) {
    for (int j = 0; j < n_gy; j++) {
      double sum = 0.0;
      for (int k = 0; k < ndata; k++) {
        sum += dx[i * ndata + k] * dy[j * ndata + k];
      }
      z(i, j) = sum / denom;
    }
  }

  return List::create(Named("x") = gx,
                      Named("y") = gy,
                      Named("z") = z);
}
