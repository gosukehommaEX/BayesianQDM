// Region probabilities of two continuous endpoints (Yamaguchi et al.,
// submitted, Appendix A.2) for many pairs of distributions of the treatment
// and control group means at once.
//
// The mean of each group has a bivariate t distribution t_nu(loc, V)
// (posterior distribution of the mean, or posterior predictive distribution
// of the mean of future subjects). The effect theta = mean_t - mean_c is
// classified into the nine regions of Table 1 of the article: region
// 3 c_1 + c_2 (0-based), where c_k is 0 if theta_k > TV_k, 1 if
// MAV_k < theta_k <= TV_k and 2 if theta_k <= MAV_k. The four regions of a
// predictive criterion are the regions 0, 2, 6 and 8 with TV_k = MAV_k =
// null value.
//
// Two methods are provided:
// - moment matching: theta is approximated by a bivariate t distribution
//   with the moments of Yamaguchi et al. (submitted, Moment_matching
//   article, Theorem 3), and the probabilities of the regions are obtained
//   from its distribution function, computed for real degrees of freedom by
//   Gauss-Legendre quadrature (see bvt_cdf());
// - Monte Carlo simulation: nMC draws from each distribution, independently
//   for each pair.

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace {

// Positive nodes and weights of the 48-point Gauss-Legendre rule on [-1, 1]
const double kGLx[24] = {
    0.032380170962869367,
    0.097004699209462697,
    0.16122235606889171,
    0.22476379039468905,
    0.28736248735545555,
    0.3487558862921607,
    0.40868648199071672,
    0.46690290475095841,
    0.523160974722233,
    0.57722472608397268,
    0.6288673967765136,
    0.67787237963266389,
    0.72403413092381463,
    0.76715903251574036,
    0.80706620402944262,
    0.84358826162439349,
    0.87657202027424785,
    0.90587913671556963,
    0.93138669070655433,
    0.9529877031604308,
    0.97059159254624727,
    0.98412458372282685,
    0.99353017226635076,
    0.99877100725242607};
const double kGLw[24] = {
    0.064737696812683654,
    0.064466164435949824,
    0.063924238584647866,
    0.063114192286253729,
    0.062039423159892422,
    0.060704439165893603,
    0.059114839698395344,
    0.057277292100402916,
    0.055199503699984033,
    0.052890189485193438,
    0.050359035553854202,
    0.047616658492490242,
    0.044674560856694231,
    0.041545082943464547,
    0.038241351065830501,
    0.034777222564770394,
    0.031167227832798097,
    0.027426509708357034,
    0.023570760839324047,
    0.019616160457356056,
    0.015579315722943226,
    0.011477234579234614,
    0.007327553901276135,
    0.0031533460523098414};

const double kHalfPi = 1.5707963267948966;

// Integral over phi in (a, b) of the integrand of bvt_cdf() with the
// 48-point Gauss-Legendre rule.
double bvt_piece(double a, double b, double k, double rho, double nu) {
  const double half = 0.5 * (b - a);
  const double mid = 0.5 * (b + a);
  const double s = std::sqrt((1.0 - rho * rho) / (nu + 1.0));
  double total = 0.0;
  for (int i = 0; i < 24; ++i) {
    for (int sign = -1; sign <= 1; sign += 2) {
      const double phi = mid + sign * half * kGLx[i];
      const double u = std::tan(phi);
      const double c = std::cos(phi);
      const double scale = s * std::sqrt(nu + u * u);
      const double num = k - rho * u;
      double cond;
      if (scale > 0.0) {
        cond = R::pt(num / scale, nu + 1.0, 1, 0);
      } else {
        // # nocov start
        // Degenerate distribution (|rho| = 1): not reached for positive
        // definite scale matrices
        cond = num > 0.0 ? 1.0 : 0.0;
        // # nocov end
      }
      total += kGLw[i] * R::dt(u, nu, 0) / (c * c) * cond;
    }
  }
  return half * total;
}

// P(X <= h, Y <= k) for the standard bivariate t distribution with
// correlation rho and nu (real, positive) degrees of freedom, for finite h
// and k. Conditionally on X = u, Y has a t distribution with nu + 1 degrees
// of freedom, location rho u and scale sqrt((nu + u^2) (1 - rho^2) /
// (nu + 1)), so the probability is the integral over u < h of the density
// of X times the conditional distribution function at k. With u = tan(phi)
// the range becomes (-pi / 2, atan(h)); it is split at 0 and at the point
// where rho u = k, and each part is integrated with the 48-point
// Gauss-Legendre rule (error below 1e-10 for nu > 4 and |rho| <= 0.99; the
// error grows to about 1e-8 at |rho| = 0.999 and 1e-6 at |rho| = 0.9999,
// where the conditional distribution function changes within a narrow
// range of u).
double bvt_cdf(double h, double k, double rho, double nu) {
  const double lo = -kHalfPi;
  const double hi = std::atan(h);
  double cuts[4];
  int n_cut = 0;
  cuts[n_cut++] = lo;
  if (rho != 0.0) {
    const double c = std::atan(k / rho);
    if (c > lo && c < hi) {
      cuts[n_cut++] = c;
    }
  }
  if (0.0 < hi) {
    cuts[n_cut++] = 0.0;
  }
  std::sort(cuts + 1, cuts + n_cut);
  cuts[n_cut++] = hi;
  double total = 0.0;
  for (int i = 0; i + 1 < n_cut; ++i) {
    if (cuts[i + 1] > cuts[i]) {
      total += bvt_piece(cuts[i], cuts[i + 1], k, rho, nu);
    }
  }
  return total;
}

// Coefficients of the distribution function F(x_i, y_j) on the grid
// x = (-Inf, MAV_1, TV_1, Inf), y = (-Inf, MAV_2, TV_2, Inf) for the union
// of the regions in a mask (bit r set: region r, 0-based).
void mask_coefficients(int mask, double coef[4][4]) {
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      coef[i][j] = 0.0;
    }
  }
  for (int r = 0; r < 9; ++r) {
    if (!((mask >> r) & 1)) {
      continue;
    }
    // Band from the bottom: 0 (<= MAV), 1 (MAV, TV], 2 (> TV)
    const int bx = 2 - r / 3;
    const int by = 2 - r % 3;
    coef[bx + 1][by + 1] += 1.0;
    coef[bx][by + 1] -= 1.0;
    coef[bx + 1][by] -= 1.0;
    coef[bx][by] += 1.0;
  }
}

}  // namespace

// loc_t, loc_c: locations (one row per pair, two columns).
// scale_t, scale_c: scale matrices (one row per pair; columns V11, V12, V22).
// df_t, df_c: degrees of freedom (one per pair, all greater than 4).
// tv, mav: the thresholds TV and MAV of the two endpoints.
// masks: one integer per union of regions (bit r: region r + 1).
// Returns the probabilities of the unions by moment matching, one row per
// pair and one column per union.
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_regions_continuous2_mm(Rcpp::NumericMatrix loc_t,
                                               Rcpp::NumericMatrix scale_t,
                                               Rcpp::NumericVector df_t,
                                               Rcpp::NumericMatrix loc_c,
                                               Rcpp::NumericMatrix scale_c,
                                               Rcpp::NumericVector df_c,
                                               Rcpp::NumericVector tv,
                                               Rcpp::NumericVector mav,
                                               Rcpp::IntegerVector masks) {
  const int n = loc_t.nrow();
  const int n_mask = static_cast<int>(masks.size());
  std::vector<double> coef(static_cast<std::size_t>(n_mask) * 16);
  bool needed[4][4] = {{false}};
  for (int m = 0; m < n_mask; ++m) {
    double c[4][4];
    mask_coefficients(masks[m], c);
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        coef[static_cast<std::size_t>(m) * 16 + i * 4 + j] = c[i][j];
        if (c[i][j] != 0.0) {
          needed[i][j] = true;
        }
      }
    }
  }
  // With TV = MAV (predictive criteria) the middle band is empty and the
  // second grid line equals the first
  const bool same_x = tv[0] == mav[0];
  const bool same_y = tv[1] == mav[1];
  Rcpp::NumericMatrix out(n, n_mask);
  for (int p = 0; p < n; ++p) {
    // Moment matching (Theorem 3 with dimension 2)
    const double n1 = df_t[p];
    const double n2 = df_c[p];
    const double b1 = n1 / (n1 - 2.0);
    const double b2 = n2 / (n2 - 2.0);
    const double v1[3] = {scale_t(p, 0), scale_t(p, 1), scale_t(p, 2)};
    const double v2[3] = {scale_c(p, 0), scale_c(p, 1), scale_c(p, 2)};
    const double s11 = b1 * v1[0] + b2 * v2[0];
    const double s12 = b1 * v1[1] + b2 * v2[1];
    const double s22 = b1 * v1[2] + b2 * v2[2];
    const double det = s11 * s22 - s12 * s12;
    // A = inverse of the summed covariance matrix
    const double a11 = s22 / det;
    const double a12 = -s12 / det;
    const double a22 = s11 / det;
    // P = A V1 and Q = A V2 (2 x 2, not symmetric)
    const double p11 = a11 * v1[0] + a12 * v1[1];
    const double p12 = a11 * v1[1] + a12 * v1[2];
    const double p21 = a12 * v1[0] + a22 * v1[1];
    const double p22 = a12 * v1[1] + a22 * v1[2];
    const double q11 = a11 * v2[0] + a12 * v2[1];
    const double q12 = a11 * v2[1] + a12 * v2[2];
    const double q21 = a12 * v2[0] + a22 * v2[1];
    const double q22 = a12 * v2[1] + a22 * v2[2];
    const double t1 = p11 + p22;
    const double t2 = q11 + q22;
    const double t11 = p11 * p11 + 2.0 * p12 * p21 + p22 * p22;
    const double t22 = q11 * q11 + 2.0 * q12 * q21 + q22 * q22;
    const double t12 = p11 * q11 + p12 * q21 + p21 * q12 + p22 * q22;
    const double qm =
      (n1 * n1 / ((n1 - 2.0) * (n1 - 4.0)) * (t1 * t1 + 2.0 * t11) +
       n2 * n2 / ((n2 - 2.0) * (n2 - 4.0)) * (t2 * t2 + 2.0 * t22) +
       2.0 * n1 * n2 / ((n1 - 2.0) * (n2 - 2.0)) * (t1 * t2 + 2.0 * t12)) /
      8.0;
    // qm > 1 for finite degrees of freedom; very large degrees of freedom
    // are capped (the distribution is then normal to double precision)
    double nu = 1e10;
    if (qm > 1.0) {
      nu = std::min((2.0 - 4.0 * qm) / (1.0 - qm), 1e10);
    }
    const double f = (nu - 2.0) / nu;
    const double sd1 = std::sqrt(s11 * f);
    const double sd2 = std::sqrt(s22 * f);
    double rho = s12 / std::sqrt(s11 * s22);
    rho = std::max(-1.0, std::min(1.0, rho));
    const double d1 = loc_t(p, 0) - loc_c(p, 0);
    const double d2 = loc_t(p, 1) - loc_c(p, 1);
    const double x[3] = {0.0, (mav[0] - d1) / sd1, (tv[0] - d1) / sd1};
    const double y[3] = {0.0, (mav[1] - d2) / sd2, (tv[1] - d2) / sd2};
    double g[4][4];
    for (int i = 0; i < 4; ++i) {
      g[i][0] = 0.0;
      g[0][i] = 0.0;
    }
    g[3][3] = 1.0;
    for (int i = 1; i <= 2; ++i) {
      g[i][3] = R::pt(x[i], nu, 1, 0);
      g[3][i] = R::pt(y[i], nu, 1, 0);
    }
    bool done[3][3] = {{false}};
    for (int i = 1; i <= 2; ++i) {
      for (int j = 1; j <= 2; ++j) {
        if (!needed[i][j]) {
          g[i][j] = 0.0;
          continue;
        }
        const int si = (same_x && i == 2) ? 1 : i;
        const int sj = (same_y && j == 2) ? 1 : j;
        if (!done[si][sj]) {
          g[si][sj] = bvt_cdf(x[si], y[sj], rho, nu);
          done[si][sj] = true;
        }
        g[i][j] = g[si][sj];
      }
    }
    for (int m = 0; m < n_mask; ++m) {
      double total = 0.0;
      const double *c = &coef[static_cast<std::size_t>(m) * 16];
      for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
          if (c[i * 4 + j] != 0.0) {
            total += c[i * 4 + j] * g[i][j];
          }
        }
      }
      out(p, m) = std::min(1.0, std::max(0.0, total));
    }
    if ((p + 1) % 1000 == 0) {
      Rcpp::checkUserInterrupt();
    }
  }
  return out;
}

// Same arguments as cpp_regions_continuous2_mm(); the probabilities are
// estimated from nMC draws of each distribution, drawn with the random
// number generator of R separately for each pair (treatment then control in
// each iteration).
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_regions_continuous2_mc(Rcpp::NumericMatrix loc_t,
                                               Rcpp::NumericMatrix scale_t,
                                               Rcpp::NumericVector df_t,
                                               Rcpp::NumericMatrix loc_c,
                                               Rcpp::NumericMatrix scale_c,
                                               Rcpp::NumericVector df_c,
                                               Rcpp::NumericVector tv,
                                               Rcpp::NumericVector mav,
                                               Rcpp::IntegerVector masks,
                                               int nMC) {
  const int n = loc_t.nrow();
  const int n_mask = static_cast<int>(masks.size());
  Rcpp::NumericMatrix out(n, n_mask);
  for (int p = 0; p < n; ++p) {
    // Cholesky factors of the scale matrices
    const double lt11 = std::sqrt(scale_t(p, 0));
    const double lt21 = scale_t(p, 1) / lt11;
    const double lt22 = std::sqrt(scale_t(p, 2) - lt21 * lt21);
    const double lc11 = std::sqrt(scale_c(p, 0));
    const double lc21 = scale_c(p, 1) / lc11;
    const double lc22 = std::sqrt(scale_c(p, 2) - lc21 * lc21);
    const double nt = df_t[p];
    const double nc = df_c[p];
    long count[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    for (int u = 0; u < nMC; ++u) {
      const double zt1 = R::norm_rand();
      const double zt2 = R::norm_rand();
      const double wt = std::sqrt(R::rchisq(nt) / nt);
      const double zc1 = R::norm_rand();
      const double zc2 = R::norm_rand();
      const double wc = std::sqrt(R::rchisq(nc) / nc);
      const double th1 = loc_t(p, 0) + lt11 * zt1 / wt -
        (loc_c(p, 0) + lc11 * zc1 / wc);
      const double th2 = loc_t(p, 1) + (lt21 * zt1 + lt22 * zt2) / wt -
        (loc_c(p, 1) + (lc21 * zc1 + lc22 * zc2) / wc);
      const int c1 = th1 > tv[0] ? 0 : (th1 > mav[0] ? 1 : 2);
      const int c2 = th2 > tv[1] ? 0 : (th2 > mav[1] ? 1 : 2);
      ++count[3 * c1 + c2];
    }
    for (int m = 0; m < n_mask; ++m) {
      long total = 0;
      for (int r = 0; r < 9; ++r) {
        if ((masks[m] >> r) & 1) {
          total += count[r];
        }
      }
      out(p, m) = static_cast<double>(total) / nMC;
    }
    Rcpp::checkUserInterrupt();
  }
  return out;
}
