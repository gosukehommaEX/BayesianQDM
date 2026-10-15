// Upper tail probability of the difference of two milestone survival
// probabilities, P(S_t - S_c > q), where S_k = exp(-lambda_k t0) and the
// hazards lambda_t and lambda_c have independent Gamma distributions with
// shapes a_k and rates b_k. The probability is computed by adaptive
// Gauss-Kronrod quadrature (QUADPACK routine dqags, as used by
// stats::integrate()) over the control hazard.
//
// Given lambda_c, the event S_t - S_c > q is S_t > q + exp(-lambda_c t0),
// whose probability is 1 if q + exp(-lambda_c t0) <= 0, 0 if it is at least
// 1, and G_t(-log(q + exp(-lambda_c t0)) / t0) otherwise, where G_t is the
// distribution function of lambda_t. The probability is the integral of this
// conditional probability times the density of lambda_c over the values of
// lambda_c with 0 < q + exp(-lambda_c t0) < 1, plus P(S_c <= -q) when q < 0.
// The range is split at quantiles of the control hazard and at the values of
// lambda_c that correspond to quantiles of the treatment hazard, so that
// concentrated distributions are integrated accurately. When the control
// density is unbounded at 0 (a_c < 1), the integral below the control median
// is computed in the variable u = lambda_c^a_c, which removes the
// singularity.

#include <Rcpp.h>
#include <R_ext/Applic.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {

struct SurvDiffArgs {
  double q;
  double t0;
  double a_t;
  double b_t;
  double a_c;
  double b_c;
  // log of b_c^a_c / Gamma(a_c + 1), used with u = lambda_c^a_c
  double log_const;
  // Variable of integration: 0 = lambda_c, 1 = u = lambda_c^a_c
  int type;
};

// P(S_t > x) for the treatment hazard
inline double surv_upper(double x, const SurvDiffArgs *p) {
  if (x <= 0.0) {
    // Not reached: the range of integration excludes x <= 0 (a guard
    // against rounding at its end; excluded from test coverage)
    return 1.0;  // # nocov
  }
  if (x >= 1.0) {
    return 0.0;
  }
  return R::pgamma(-std::log(x) / p->t0, p->a_t, 1.0 / p->b_t, 1, 0);
}

// Integrand evaluated in place at the abscissae supplied by QUADPACK.
void surv_diff_integrand(double *x, int n, void *ex) {
  const SurvDiffArgs *p = static_cast<const SurvDiffArgs *>(ex);
  for (int i = 0; i < n; ++i) {
    double lambda;
    double weight;
    if (p->type == 1) {
      // lambda = u^(1 / a_c) and f_c(lambda) d lambda =
      // b_c^a_c exp(-b_c lambda) / Gamma(a_c + 1) du
      lambda = std::pow(x[i], 1.0 / p->a_c);
      weight = std::exp(p->log_const - p->b_c * lambda);
    } else {
      lambda = x[i];
      weight = R::dgamma(lambda, p->a_c, 1.0 / p->b_c, 0);
    }
    x[i] = weight * surv_upper(p->q + std::exp(-lambda * p->t0), p);
  }
}

// Integrates over [lo, hi] (lo < hi), which lies on one side of the control
// median, and records the largest error estimate and QUADPACK error code.
double integrate_segment(SurvDiffArgs &args, double lo, double hi,
                         double median, double &max_abserr, int &max_ier) {
  double a;
  double b;
  if (0.5 * (lo + hi) <= median && args.a_c < 1.0) {
    args.type = 1;
    a = std::pow(lo, args.a_c);
    b = std::pow(hi, args.a_c);
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
  Rdqags(surv_diff_integrand, &args, &a, &b, &epsabs, &epsrel, &result,
         &abserr, &neval, &ier, &limit, &lenw, &last, iwork.data(),
         work.data());
  max_abserr = std::max(max_abserr, abserr);
  max_ier = std::max(max_ier, ier);
  return result;
}

}  // namespace

// q, t0, a_t, b_t, a_c, b_c: vectors of the same length.
// Returns P(S_t - S_c > q) for each element.
// [[Rcpp::export]]
Rcpp::NumericVector cpp_psurvdiff_upper(Rcpp::NumericVector q,
                                        Rcpp::NumericVector t0,
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
    SurvDiffArgs args = {q[i], t0[i], a_t[i], b_t[i], a_c[i], b_c[i],
                         a_c[i] * std::log(b_c[i]) - R::lgammafn(a_c[i] + 1.0),
                         0};
    double value = 0.0;
    if (args.q <= -1.0) {
      value = 1.0;
    } else if (args.q < 1.0) {
      const double scale_c = 1.0 / args.b_c;
      // Range of lambda_c with 0 < q + exp(-lambda_c t0) < 1
      const double lo = args.q > 0.0 ? -std::log1p(-args.q) / args.t0 : 0.0;
      double hi = std::numeric_limits<double>::infinity();
      if (args.q < 0.0) {
        hi = -std::log(-args.q) / args.t0;
        // S_c <= -q, that is lambda_c >= hi, implies S_t - S_c > q
        value += R::pgamma(hi, args.a_c, scale_c, 0, 0);
      }
      // The control density beyond this quantile is negligible
      hi = std::min(hi, R::qgamma(1e-16, args.a_c, scale_c, 0, 0));
      if (hi > lo) {
        const double median = R::qgamma(0.5, args.a_c, scale_c, 1, 0);
        std::vector<double> cut;
        cut.push_back(lo);
        cut.push_back(hi);
        for (int k = 0; k < 3; ++k) {
          const double xc = R::qgamma(probs[k], args.a_c, scale_c, 1, 0);
          if (xc > lo && xc < hi) {
            cut.push_back(xc);
          }
          // Control hazard at which the conditional probability is at the
          // quantile of the treatment hazard
          const double st = std::exp(
            -R::qgamma(probs[k], args.a_t, 1.0 / args.b_t, 1, 0) * args.t0);
          if (st - args.q > 0.0 && st - args.q < 1.0) {
            const double xt = -std::log(st - args.q) / args.t0;
            if (xt > lo && xt < hi) {
              cut.push_back(xt);
            }
          }
        }
        std::sort(cut.begin(), cut.end());
        cut.erase(std::unique(cut.begin(), cut.end()), cut.end());
        for (std::size_t k = 0; k + 1 < cut.size(); ++k) {
          value += integrate_segment(args, cut[k], cut[k + 1], median,
                                     max_abserr, max_ier);
        }
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
