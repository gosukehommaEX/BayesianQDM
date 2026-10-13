// Upper tail probability of the difference of two independent location-scale
// t variables, P(T_t - T_c > q), computed by adaptive Gauss-Kronrod quadrature
// over the whole real line (QUADPACK routines dqags and dqagi, as used by
// stats::integrate()).
//
// With T_k = mu_k + s_k X_k, P(T_t - T_c > q) = E[P(X_t > a + b U)], where
// U ~ t(df_c), a = (mu_c + q - mu_t) / s_t and b = s_c / s_t.

#include <Rcpp.h>
#include <R_ext/Applic.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace {

struct TDiffArgs {
  double a;
  double b;
  double df_t;
  double df_c;
};

// Integrand f(u; df_c) * P(X_t > a + b u). QUADPACK passes a vector of
// abscissae that is overwritten with the values of the integrand.
void t_diff_integrand(double *x, int n, void *ex) {
  const TDiffArgs *p = static_cast<const TDiffArgs *>(ex);
  for (int i = 0; i < n; ++i) {
    const double u = x[i];
    x[i] = R::dt(u, p->df_c, 0) * R::pt(p->a + p->b * u, p->df_t, 0, 0);
  }
}

const double kEpsAbs = 1e-13;
const double kEpsRel = 1e-11;
const int kLimit = 200;

// Integral over the finite interval [lo, hi].
double integrate_finite(TDiffArgs &args, double lo, double hi,
                        double &max_abserr, int &max_ier) {
  if (!(hi > lo)) {
    return 0.0;
  }
  double a = lo;
  double b = hi;
  double epsabs = kEpsAbs;
  double epsrel = kEpsRel;
  double result = 0.0;
  double abserr = 0.0;
  int neval = 0;
  int ier = 0;
  int limit = kLimit;
  int lenw = 4 * limit;
  int last = 0;
  std::vector<int> iwork(limit);
  std::vector<double> work(lenw);
  Rdqags(t_diff_integrand, &args, &a, &b, &epsabs, &epsrel, &result, &abserr,
         &neval, &ier, &limit, &lenw, &last, iwork.data(), work.data());
  max_abserr = std::max(max_abserr, abserr);
  max_ier = std::max(max_ier, ier);
  return result;
}

// Integral over (-Inf, bound] (inf = -1) or [bound, Inf) (inf = 1).
double integrate_infinite(TDiffArgs &args, double bound, int inf,
                          double &max_abserr, int &max_ier) {
  double bnd = bound;
  int infinite = inf;
  double epsabs = kEpsAbs;
  double epsrel = kEpsRel;
  double result = 0.0;
  double abserr = 0.0;
  int neval = 0;
  int ier = 0;
  int limit = kLimit;
  int lenw = 4 * limit;
  int last = 0;
  std::vector<int> iwork(limit);
  std::vector<double> work(lenw);
  Rdqagi(t_diff_integrand, &args, &bnd, &infinite, &epsabs, &epsrel, &result,
         &abserr, &neval, &ier, &limit, &lenw, &last, iwork.data(),
         work.data());
  max_abserr = std::max(max_abserr, abserr);
  max_ier = std::max(max_ier, ier);
  return result;
}

}  // namespace

// [[Rcpp::export]]
Rcpp::NumericVector cpp_ptdiff_ni_upper(Rcpp::NumericVector q,
                                        Rcpp::NumericVector mu_t,
                                        Rcpp::NumericVector scale_t,
                                        Rcpp::NumericVector df_t,
                                        Rcpp::NumericVector mu_c,
                                        Rcpp::NumericVector scale_c,
                                        Rcpp::NumericVector df_c) {
  const R_xlen_t n = q.size();
  Rcpp::NumericVector out(n);
  double max_abserr = 0.0;
  int max_ier = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    TDiffArgs args = {(mu_c[i] + q[i] - mu_t[i]) / scale_t[i],
                      scale_c[i] / scale_t[i], df_t[i], df_c[i]};
    // Split at the peak of the control density (u = 0) and at the point
    // where P(X_t > a + b u) = 1/2 (u = -a / b).
    const double u_half = -args.a / args.b;
    const double lo = std::min(0.0, u_half);
    const double hi = std::max(0.0, u_half);
    double value = integrate_infinite(args, lo, -1, max_abserr, max_ier);
    value += integrate_finite(args, lo, hi, max_abserr, max_ier);
    value += integrate_infinite(args, hi, 1, max_abserr, max_ier);
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
