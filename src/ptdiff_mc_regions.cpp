// Monte Carlo estimates of the probabilities of the Go and NoGo regions of
// one continuous endpoint from the same draws of the difference of two
// independent location-scale t variables:
//   go = P(T_t - T_c > q_go),  nogo = P(T_t - T_c <= q_nogo).
// Because both probabilities use the same draws, go + nogo = 1 exactly when
// q_go = q_nogo (posterior predictive criterion) and go + nogo <= 1 when
// q_go > q_nogo (posterior criterion), as for the exact probabilities.

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_ptdiff_mc_regions(Rcpp::NumericVector q_go,
                                          Rcpp::NumericVector q_nogo,
                                          Rcpp::NumericVector mu_t,
                                          Rcpp::NumericVector scale_t,
                                          Rcpp::NumericVector df_t,
                                          Rcpp::NumericVector mu_c,
                                          Rcpp::NumericVector scale_c,
                                          Rcpp::NumericVector df_c,
                                          int nMC) {
  const int n = static_cast<int>(q_go.size());
  Rcpp::NumericMatrix out(n, 2);
  for (int i = 0; i < n; ++i) {
    // T_t - T_c compared with q  <=>  s_t X_t - s_c X_c compared with
    // q - (mu_t - mu_c)
    const double shift = mu_t[i] - mu_c[i];
    const double thr_go = q_go[i] - shift;
    const double thr_nogo = q_nogo[i] - shift;
    long count_go = 0;
    long count_nogo = 0;
    for (int j = 0; j < nMC; ++j) {
      const double d = scale_t[i] * R::rt(df_t[i]) -
        scale_c[i] * R::rt(df_c[i]);
      if (d > thr_go) {
        ++count_go;
      }
      if (d <= thr_nogo) {
        ++count_nogo;
      }
    }
    out(i, 0) = static_cast<double>(count_go) / nMC;
    out(i, 1) = static_cast<double>(count_nogo) / nMC;
    Rcpp::checkUserInterrupt();
  }
  return out;
}
