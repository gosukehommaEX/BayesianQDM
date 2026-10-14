// Posterior probabilities of unions of the nine regions of two binary
// endpoints (Yamaguchi et al., submitted, Table 1) by Monte Carlo simulation,
// for many pairs of treatment and control posterior distributions at once.
//
// The joint response probabilities of each group have a Dirichlet posterior
// distribution with parameters (a_00, a_01, a_10, a_11). In each of nMC
// iterations, one draw is made from every treatment and every control
// posterior distribution (as normalized Gamma variables, using the random
// number generator of R), and the marginal response rates
// pi_1 = p_10 + p_11 and pi_2 = p_01 + p_11 are computed. Every requested
// pair of treatment and control distributions is then classified into one of
// the nine regions of theta_k = pi_tk - pi_ck (k = 1, 2), so all pairs share
// the same draws. The region of a pair is 3 c_1 + c_2 (0-based), where c_k is
// 0 if theta_k > TV_k, 1 if MAV_k < theta_k <= TV_k and 2 if
// theta_k <= MAV_k.

#include <Rcpp.h>

#include <vector>

namespace {

// Marginal response rates of one Dirichlet draw with parameters a[0..3] in
// the cell order 00, 01, 10, 11.
inline void draw_margins(const double *a, double &pi1, double &pi2) {
  double g[4];
  double total = 0.0;
  for (int k = 0; k < 4; ++k) {
    g[k] = R::rgamma(a[k], 1.0);
    total += g[k];
  }
  pi1 = (g[2] + g[3]) / total;
  pi2 = (g[1] + g[3]) / total;
}

inline int band(double theta, double tv, double mav) {
  return theta > tv ? 0 : (theta > mav ? 1 : 2);
}

}  // namespace

// alpha_t, alpha_c: Dirichlet parameters, one row per distribution.
// pair_t, pair_c: 0-based row indices of the pairs to evaluate.
// tv, mav: thresholds of the two endpoints.
// masks: one integer per union of regions; bit r (0 to 8) is set if region
//   r + 1 belongs to the union.
// Returns the estimated probabilities, one row per pair and one column per
// union.
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_post_regions_binary2(Rcpp::NumericMatrix alpha_t,
                                             Rcpp::NumericMatrix alpha_c,
                                             Rcpp::IntegerVector pair_t,
                                             Rcpp::IntegerVector pair_c,
                                             Rcpp::NumericVector tv,
                                             Rcpp::NumericVector mav,
                                             Rcpp::IntegerVector masks,
                                             int nMC) {
  const int n_t = alpha_t.nrow();
  const int n_c = alpha_c.nrow();
  const int n_pair = static_cast<int>(pair_t.size());
  const int n_mask = static_cast<int>(masks.size());
  // Parameters stored row by row for fast access
  std::vector<double> a_t(4 * static_cast<std::size_t>(n_t));
  std::vector<double> a_c(4 * static_cast<std::size_t>(n_c));
  for (int i = 0; i < n_t; ++i) {
    for (int k = 0; k < 4; ++k) {
      a_t[4 * static_cast<std::size_t>(i) + k] = alpha_t(i, k);
    }
  }
  for (int j = 0; j < n_c; ++j) {
    for (int k = 0; k < 4; ++k) {
      a_c[4 * static_cast<std::size_t>(j) + k] = alpha_c(j, k);
    }
  }
  // For each region, the unions that contain it
  std::vector<std::vector<int> > in_mask(9);
  for (int r = 0; r < 9; ++r) {
    for (int m = 0; m < n_mask; ++m) {
      if ((masks[m] >> r) & 1) {
        in_mask[r].push_back(m);
      }
    }
  }
  std::vector<double> pt1(n_t), pt2(n_t), pc1(n_c), pc2(n_c);
  std::vector<int> count(static_cast<std::size_t>(n_pair) * n_mask, 0);
  for (int u = 0; u < nMC; ++u) {
    for (int i = 0; i < n_t; ++i) {
      draw_margins(&a_t[4 * static_cast<std::size_t>(i)], pt1[i], pt2[i]);
    }
    for (int j = 0; j < n_c; ++j) {
      draw_margins(&a_c[4 * static_cast<std::size_t>(j)], pc1[j], pc2[j]);
    }
    for (int p = 0; p < n_pair; ++p) {
      const int i = pair_t[p];
      const int j = pair_c[p];
      const int r = 3 * band(pt1[i] - pc1[j], tv[0], mav[0]) +
        band(pt2[i] - pc2[j], tv[1], mav[1]);
      const std::vector<int> &ms = in_mask[r];
      for (std::size_t q = 0; q < ms.size(); ++q) {
        ++count[static_cast<std::size_t>(ms[q]) * n_pair + p];
      }
    }
    if ((u + 1) % 100 == 0) {
      Rcpp::checkUserInterrupt();
    }
  }
  Rcpp::NumericMatrix out(n_pair, n_mask);
  for (int m = 0; m < n_mask; ++m) {
    for (int p = 0; p < n_pair; ++p) {
      out(p, m) = static_cast<double>(
        count[static_cast<std::size_t>(m) * n_pair + p]) / nMC;
    }
  }
  return out;
}
