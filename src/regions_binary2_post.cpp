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
//
// The draws are stored for blocks of iterations and every pair is then
// classified for the whole block, so that its counts stay in local variables;
// the draws and the counts are the same as when the pairs are classified
// after each iteration.
//
// With a fixed control (no control distributions), the marginal response
// rates of the control group are the fixed values, and only the treatment
// distributions are drawn.

#include <Rcpp.h>

#include <algorithm>
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

// Band of theta: 0 if theta > tv, 1 if mav < theta <= tv and 2 if
// theta <= mav (tv > mav, as required by criterion_posterior()), computed
// without branches.
inline int band(double theta, double tv, double mav) {
  return 2 - static_cast<int>(theta > mav) - static_cast<int>(theta > tv);
}

}  // namespace

// alpha_t, alpha_c: Dirichlet parameters, one row per distribution.
// pair_t, pair_c: 0-based row indices of the pairs to evaluate.
// tv, mav: thresholds of the two endpoints.
// masks: one integer per union of regions; bit r (0 to 8) is set if region
//   r + 1 belongs to the union.
// value_c: marginal response rates of a fixed control, used when alpha_c has
//   no rows (pair_c is then not used).
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
                                             int nMC,
                                             Rcpp::NumericVector value_c) {
  const int n_t = alpha_t.nrow();
  const int n_c = alpha_c.nrow();
  const bool fixed = n_c == 0;
  const double vc1 = fixed ? value_c[0] : 0.0;
  const double vc2 = fixed ? value_c[1] : 0.0;
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
  // Margins of the draws of a block of iterations, stored by distribution
  // (row i holds the draws of distribution i in the iterations of the block)
  const int block = 64;
  std::vector<double> pt1(static_cast<std::size_t>(n_t) * block);
  std::vector<double> pt2(static_cast<std::size_t>(n_t) * block);
  std::vector<double> pc1(static_cast<std::size_t>(n_c) * block);
  std::vector<double> pc2(static_cast<std::size_t>(n_c) * block);
  // Counts by pair (row) and union (column)
  std::vector<int> count(static_cast<std::size_t>(n_pair) * n_mask, 0);
  const double tv1 = tv[0];
  const double tv2 = tv[1];
  const double mav1 = mav[0];
  const double mav2 = mav[1];
  for (int u0 = 0; u0 < nMC; u0 += block) {
    const int nb = std::min(block, nMC - u0);
    // Same order of draws as one iteration at a time: all treatment, then
    // all control distributions
    for (int b = 0; b < nb; ++b) {
      for (int i = 0; i < n_t; ++i) {
        const std::size_t k = static_cast<std::size_t>(i) * block + b;
        draw_margins(&a_t[4 * static_cast<std::size_t>(i)], pt1[k], pt2[k]);
      }
      for (int j = 0; j < n_c; ++j) {
        const std::size_t k = static_cast<std::size_t>(j) * block + b;
        draw_margins(&a_c[4 * static_cast<std::size_t>(j)], pc1[k], pc2[k]);
      }
    }
    for (int p = 0; p < n_pair; ++p) {
      const std::size_t i = static_cast<std::size_t>(pair_t[p]) * block;
      int hist[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
      if (fixed) {
        for (int b = 0; b < nb; ++b) {
          ++hist[3 * band(pt1[i + b] - vc1, tv1, mav1) +
                 band(pt2[i + b] - vc2, tv2, mav2)];
        }
      } else {
        const std::size_t j = static_cast<std::size_t>(pair_c[p]) * block;
        for (int b = 0; b < nb; ++b) {
          ++hist[3 * band(pt1[i + b] - pc1[j + b], tv1, mav1) +
                 band(pt2[i + b] - pc2[j + b], tv2, mav2)];
        }
      }
      int *cnt = &count[static_cast<std::size_t>(p) * n_mask];
      for (int r = 0; r < 9; ++r) {
        const std::vector<int> &ms = in_mask[r];
        for (std::size_t q = 0; q < ms.size(); ++q) {
          cnt[ms[q]] += hist[r];
        }
      }
    }
    Rcpp::checkUserInterrupt();
  }
  Rcpp::NumericMatrix out(n_pair, n_mask);
  for (int m = 0; m < n_mask; ++m) {
    for (int p = 0; p < n_pair; ++p) {
      out(p, m) = static_cast<double>(
        count[static_cast<std::size_t>(p) * n_mask + m]) / nMC;
    }
  }
  return out;
}
