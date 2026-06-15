/* Regularized horseshoe prior for shared covariate effects.
 *
 * See Appendix C.1 in Piironen and Vehtari (2017):
 * https://projecteuclid.org/euclid.ejs/1513306866
 *
 * Args:
 *   z: standardized population-level coefficients
 *   lambda: local shrinkage parameters
 *   tau: global shrinkage parameter
 *   c2: slab regularization parameter
 *
 * Returns:
 *   population-level coefficients following the regularized horseshoe prior
 */
vector horseshoe(vector z, vector lambda, real tau, real c2) {
  int K = rows(z);
  vector[K] lambda2 = square(lambda);
  vector[K] lambda_tilde = sqrt(c2 * lambda2 ./ (c2 + tau^2 * lambda2));
  return z .* lambda_tilde * tau;
}
