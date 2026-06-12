/* Miscellaneous response-link helpers shared by Stan SDM templates. */

/* Inverse complementary log-log response transform.
 *
 * Args:
 *   eta: linear predictor on the cloglog link scale
 *
 * Returns:
 *   probability on the response scale
 */
real inv_cloglog_response(real eta) {
  return 1 - exp(-exp(eta));
}
