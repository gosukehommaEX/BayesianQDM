// Monte Carlo estimate of the upper tail probability of the difference of two
// independent location-scale t variables, P(T_t - T_c > q), using the random
// number generator of R.

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericVector cpp_ptdiff_mc_upper(Rcpp::NumericVector q,
                                        Rcpp::NumericVector mu_t,
                                        Rcpp::NumericVector scale_t,
                                        Rcpp::NumericVector df_t,
                                        Rcpp::NumericVector mu_c,
                                        Rcpp::NumericVector scale_c,
                                        Rcpp::NumericVector df_c,
                                        int nMC) {
  const R_xlen_t n = q.size();
  Rcpp::NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    // T_t - T_c > q  <=>  s_t X_t - s_c X_c > q - (mu_t - mu_c)
    const double threshold = q[i] - (mu_t[i] - mu_c[i]);
    long count = 0;
    for (int j = 0; j < nMC; ++j) {
      const double d = scale_t[i] * R::rt(df_t[i]) -
        scale_c[i] * R::rt(df_c[i]);
      if (d > threshold) {
        ++count;
      }
    }
    out[i] = static_cast<double>(count) / nMC;
    Rcpp::checkUserInterrupt();
  }
  return out;
}
