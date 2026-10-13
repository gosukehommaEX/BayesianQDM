// Upper tail probability of the difference of two independent location-scale
// t variables, P(T_t - T_c > q), computed by adaptive Gauss-Kronrod quadrature
// (QUADPACK routine dqags, as used by stats::integrate()).
//
// With T_k = mu_k + s_k X_k, the probability is written as an integral over
// the standardised variable of the group with the smaller scale:
//   s_c <= s_t: P = int f_c(u) P(X_t > a + b u) du,
//               a = (mu_c + q - mu_t) / s_t, b = s_c / s_t;
//   s_c >  s_t: P = int f_t(u) P(X_c < a + b u) du,
//               a = (mu_t - mu_c - q) / s_c, b = s_t / s_c;
// so that b <= 1 and the second factor changes no faster than the density.
// The real line is split at -1, 0 and 1, at the midpoint u_half = -a / b of
// the change of the second factor, at u_half - 1 / b and u_half + 1 / b, and
// at the mirror images of these three points (where |b u| is comparable to
// |a|). Segments outside (-1, 1) whose ends differ by more than a factor of 2
// are integrated in log(|u|), and the two tails beyond the outermost points c
// are integrated in w = c / u, which runs over (0, 1]. In these variables
// heavy tails and features far from the origin are integrated accurately.

#include <Rcpp.h>
#include <R_ext/Applic.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace {

struct TDiffArgs {
  double a;
  double b;
  double df_int;   // degrees of freedom of the variable of integration
  double df_oth;   // degrees of freedom of the other variable
  int lower;       // 1: P(X < a + b u); 0: P(X > a + b u)
  // Variable of integration: 0 = u, 1 = s with u = exp(s),
  // 2 = s with u = -exp(s), 3 = w with u = c / w
  int type;
  double c;
};

// Integrand evaluated in place at the abscissae supplied by QUADPACK. The
// density and the Jacobian are combined on the log scale, so that extreme
// values of the variable give 0 instead of 0 * Inf or NaN.
void t_diff_integrand(double *x, int n, void *ex) {
  const TDiffArgs *p = static_cast<const TDiffArgs *>(ex);
  for (int i = 0; i < n; ++i) {
    double u;
    double log_jac;
    if (p->type == 1) {
      u = std::exp(x[i]);
      log_jac = x[i];
    } else if (p->type == 2) {
      u = -std::exp(x[i]);
      log_jac = x[i];
    } else if (p->type == 3) {
      u = p->c / x[i];
      log_jac = std::log(std::fabs(p->c)) - 2.0 * std::log(x[i]);
    } else {
      u = x[i];
      log_jac = 0.0;
    }
    const double log_dens = R::dt(u, p->df_int, 1) + log_jac;
    x[i] = std::isfinite(log_dens) ?
      std::exp(log_dens) * R::pt(p->a + p->b * u, p->df_oth, p->lower, 0) :
      0.0;
  }
}

// Integral of the integrand in the variable given by args.type over [lo, hi].
double integrate_piece(TDiffArgs &args, double lo, double hi,
                       double &max_abserr, int &max_ier) {
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
  Rdqags(t_diff_integrand, &args, &a, &b, &epsabs, &epsrel, &result, &abserr,
         &neval, &ier, &limit, &lenw, &last, iwork.data(), work.data());
  max_abserr = std::max(max_abserr, abserr);
  max_ier = std::max(max_ier, ier);
  return result;
}

// Integral over the whole real line, split as described at the top.
double integrate_line(TDiffArgs &args, double &max_abserr, int &max_ier) {
  std::vector<double> cut;
  cut.push_back(-1.0);
  cut.push_back(0.0);
  cut.push_back(1.0);
  const double u_half = -args.a / args.b;
  const double centres[2] = {u_half, -u_half};
  const double offsets[3] = {0.0, -1.0 / args.b, 1.0 / args.b};
  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 3; ++j) {
      const double point = centres[k] + offsets[j];
      if (std::isfinite(point)) {
        cut.push_back(point);
      }
    }
  }
  std::sort(cut.begin(), cut.end());
  cut.erase(std::unique(cut.begin(), cut.end()), cut.end());
  double value = 0.0;
  // Tails beyond the outermost points, in w = c / u
  args.type = 3;
  args.c = cut.back();
  value += integrate_piece(args, 0.0, 1.0, max_abserr, max_ier);
  args.c = cut.front();
  value += integrate_piece(args, 0.0, 1.0, max_abserr, max_ier);
  for (std::size_t k = 0; k + 1 < cut.size(); ++k) {
    const double lo = cut[k];
    const double hi = cut[k + 1];
    if (lo >= 1.0 && hi > 2.0 * lo) {
      args.type = 1;
      value += integrate_piece(args, std::log(lo), std::log(hi), max_abserr,
                               max_ier);
    } else if (hi <= -1.0 && lo < 2.0 * hi) {
      args.type = 2;
      value += integrate_piece(args, std::log(-hi), std::log(-lo),
                               max_abserr, max_ier);
    } else {
      args.type = 0;
      value += integrate_piece(args, lo, hi, max_abserr, max_ier);
    }
  }
  return value;
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
    TDiffArgs args;
    if (scale_c[i] <= scale_t[i]) {
      args.a = (mu_c[i] + q[i] - mu_t[i]) / scale_t[i];
      args.b = scale_c[i] / scale_t[i];
      args.df_int = df_c[i];
      args.df_oth = df_t[i];
      args.lower = 0;
    } else {
      args.a = (mu_t[i] - mu_c[i] - q[i]) / scale_c[i];
      args.b = scale_t[i] / scale_c[i];
      args.df_int = df_t[i];
      args.df_oth = df_c[i];
      args.lower = 1;
    }
    args.type = 0;
    args.c = 0.0;
    const double value = integrate_line(args, max_abserr, max_ier);
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
