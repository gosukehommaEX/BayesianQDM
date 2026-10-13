// Upper tail probability of the difference of two independent Beta
// variables, P(pi_t - pi_c > q), computed by adaptive Gauss-Kronrod
// quadrature (QUADPACK routine dqags, as used by stats::integrate()).
//
// The probability is the integral of f_c(x) P(pi_t > x + q) over the values
// of x with 0 <= x + q <= 1, plus P(pi_c < -q) when q < 0. The range is split
// at quantiles of both distributions so that concentrated distributions are
// integrated accurately. When the control density is unbounded at 0
// (a_c < 1) or at 1 (b_c < 1), the integral on that side of the control
// median is computed in the variable s = x^a_c or s = (1 - x)^b_c, which
// removes the singularity. Without this change of variables, QUADPACK can
// return inaccurate results with small error estimates when a split point
// lies close to the singularity.

#include <Rcpp.h>
#include <R_ext/Applic.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace {

struct BetaDiffArgs {
  double q;
  double a_t;
  double b_t;
  double a_c;
  double b_c;
  double log_a_c;
  double log_b_c;
  double lbeta_c;
  // Variable of integration: 0 = x, 1 = s = x^a_c, 2 = s = (1 - x)^b_c
  int type;
};

// Integrand evaluated in place at the abscissae supplied by QUADPACK.
void beta_diff_integrand(double *x, int n, void *ex) {
  const BetaDiffArgs *p = static_cast<const BetaDiffArgs *>(ex);
  for (int i = 0; i < n; ++i) {
    double xi;
    double weight;
    if (p->type == 1) {
      // x = s^(1 / a_c) and f_c(x) dx = (1 - x)^(b_c - 1) / (a_c B) ds
      xi = std::pow(x[i], 1.0 / p->a_c);
      weight = std::exp((p->b_c - 1.0) * std::log1p(-xi) - p->log_a_c -
                        p->lbeta_c);
    } else if (p->type == 2) {
      // 1 - x = s^(1 / b_c) and f_c(x) |dx| = x^(a_c - 1) / (b_c B) ds
      xi = -std::expm1(std::log(x[i]) / p->b_c);
      weight = std::exp((p->a_c - 1.0) * std::log(xi) - p->log_b_c -
                        p->lbeta_c);
    } else {
      xi = x[i];
      weight = R::dbeta(xi, p->a_c, p->b_c, 0);
    }
    x[i] = weight * R::pbeta(xi + p->q, p->a_t, p->b_t, 0, 0);
  }
}

// Integrates over [lo, hi] (lo < hi), which lies on one side of the control
// median, and records the largest error estimate and QUADPACK error code.
double integrate_segment(BetaDiffArgs &args, double lo, double hi,
                         double median, double &max_abserr, int &max_ier) {
  const double mid = 0.5 * (lo + hi);
  double a;
  double b;
  if (mid <= median && args.a_c < 1.0) {
    args.type = 1;
    a = std::pow(lo, args.a_c);
    b = std::pow(hi, args.a_c);
  } else if (mid > median && args.b_c < 1.0) {
    args.type = 2;
    a = std::pow(1.0 - hi, args.b_c);
    b = std::pow(1.0 - lo, args.b_c);
  } else {
    args.type = 0;
    a = lo;
    b = hi;
  }
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
    BetaDiffArgs args = {q[i], a_t[i], b_t[i], a_c[i], b_c[i],
                         std::log(a_c[i]), std::log(b_c[i]),
                         R::lbeta(a_c[i], b_c[i]), 0};
    const double lo = std::max(0.0, -args.q);
    const double hi = std::min(1.0, 1.0 - args.q);
    const double median = R::qbeta(0.5, args.a_c, args.b_c, 1, 0);
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
        value += integrate_segment(args, cut[k], cut[k + 1], median,
                                   max_abserr, max_ier);
      }
    }
    out[i] = std::min(1.0, std::max(0.0, value));
    if ((i + 1) % 1000 == 0) {
      Rcpp::checkUserInterrupt();
    }
  }
  // Warning for extreme parameter values for which QUADPACK reports a
  // problem; not triggered by the tests (excluded from test coverage).
  // # nocov start
  if (max_ier > 0 && max_abserr > 1e-8) {
    Rcpp::warning("Numerical integration may be inaccurate "
                  "(estimated absolute error %g).", max_abserr);
  }
  // # nocov end
  return out;
}
