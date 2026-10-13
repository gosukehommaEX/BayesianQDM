// Upper tail probability of the difference of two independent Beta
// variables, P(pi_t - pi_c > q), computed by adaptive Gauss-Kronrod
// quadrature (QUADPACK routine dqags, as used by stats::integrate()).

#include <Rcpp.h>
#include <R_ext/Applic.h>

#include <algorithm>
#include <vector>

namespace {

struct BetaDiffArgs {
  double q;
  double a_t;
  double b_t;
  double a_c;
  double b_c;
};

// Integrand f_c(x) * P(pi_t > x + q). QUADPACK passes a vector of abscissae
// that is overwritten with the values of the integrand.
void beta_diff_integrand(double *x, int n, void *ex) {
  const BetaDiffArgs *p = static_cast<const BetaDiffArgs *>(ex);
  for (int i = 0; i < n; ++i) {
    const double xi = x[i];
    x[i] = R::dbeta(xi, p->a_c, p->b_c, 0) *
      R::pbeta(xi + p->q, p->a_t, p->b_t, 0, 0);
  }
}

// Integrates over [lo, hi] and records the largest error estimate and the
// largest QUADPACK error code.
double integrate_segment(BetaDiffArgs &args, double lo, double hi,
                         double &max_abserr, int &max_ier) {
  if (!(hi > lo)) {
    return 0.0;
  }
  double a = lo;
  double b = hi;
  double epsabs = 1e-13;
  double epsrel = 1e-11;
  double result = 0.0;
  double abserr = 0.0;
  int neval = 0;
  int ier = 0;
  int limit = 200;
  int lenw = 4 * limit;
  int last = 0;
  std::vector<int> iwork(limit);
  std::vector<double> work(lenw);
  Rdqags(beta_diff_integrand, &args, &a, &b, &epsabs, &epsrel, &result,
         &abserr, &neval, &ier, &limit, &lenw, &last, iwork.data(),
         work.data());
  max_abserr = std::max(max_abserr, abserr);
  max_ier = std::max(max_ier, ier);
  return result;
}

}  // namespace

// [[Rcpp::export]]
Rcpp::NumericVector cpp_pbetadiff_upper(Rcpp::NumericVector q,
                                        Rcpp::NumericVector a_t,
                                        Rcpp::NumericVector b_t,
                                        Rcpp::NumericVector a_c,
                                        Rcpp::NumericVector b_c) {
  const R_xlen_t n = q.size();
  Rcpp::NumericVector out(n);
  // Quantiles used to split the range of integration
  const double probs[3] = {1e-12, 0.5, 1.0 - 1e-12};
  double max_abserr = 0.0;
  int max_ier = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    BetaDiffArgs args = {q[i], a_t[i], b_t[i], a_c[i], b_c[i]};
    const double lo = std::max(0.0, -args.q);
    const double hi = std::min(1.0, 1.0 - args.q);
    double value = 0.0;
    if (args.q < 0.0) {
      // pi_c < -q implies pi_t - pi_c > q
      value += R::pbeta(-args.q, args.a_c, args.b_c, 1, 0);
    }
    if (hi > lo) {
      std::vector<double> cut;
      cut.push_back(lo);
      cut.push_back(hi);
      for (int k = 0; k < 3; ++k) {
        const double xc = R::qbeta(probs[k], args.a_c, args.b_c, 1, 0);
        const double xt = R::qbeta(probs[k], args.a_t, args.b_t, 1, 0) -
          args.q;
        if (xc > lo && xc < hi) {
          cut.push_back(xc);
        }
        if (xt > lo && xt < hi) {
          cut.push_back(xt);
        }
      }
      std::sort(cut.begin(), cut.end());
      cut.erase(std::unique(cut.begin(), cut.end()), cut.end());
      for (std::size_t k = 0; k + 1 < cut.size(); ++k) {
        value += integrate_segment(args, cut[k], cut[k + 1], max_abserr,
                                   max_ier);
      }
    }
    out[i] = std::min(1.0, std::max(0.0, value));
    if ((i + 1) % 1000 == 0) {
      Rcpp::checkUserInterrupt();
    }
  }
  if (max_ier > 0 && max_abserr > 1e-8) {
    Rcpp::warning("Numerical integration may be inaccurate "
                  "(estimated absolute error %g).", max_abserr);
  }
  return out;
}
