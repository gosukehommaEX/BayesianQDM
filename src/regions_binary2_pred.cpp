// Posterior predictive probabilities of unions of the four regions of two
// binary endpoints (Yamaguchi et al., submitted, Table 2), computed exactly,
// for many pairs of treatment and control posterior distributions at once.
//
// Given Dirichlet parameters a, the future counts of the four cells among m
// subjects follow the Dirichlet-multinomial distribution. The future effects
// depend only on the numbers of responders S_1 = X_10 + X_11 and
// S_2 = X_01 + X_11, so the joint distribution of (S_1, S_2) is computed
// once for every distribution by summing the Dirichlet-multinomial
// probabilities over X_11. The future effect of endpoint k is
// S_tk / m_t - S_ck / m_c; it exceeds the null value q_k if and only if
// S_ck <= K_k(S_tk), where K_k(s) is the largest integer with
// (m_c s - m_t K_k(s)) / (m_t m_c) > q_k (this comparison of the exact
// difference with q_k is the same as in pbetabinomdiff(), so that a future
// effect equal to q_k is not counted as larger). The probabilities of the
// regions of a pair are then sums over (S_t1, S_t2) of the treatment
// probabilities times values of the bivariate distribution function of the
// control counts.
//
// Regions (0-based): 0 = both effects above the null values, 1 = endpoint 1
// above and endpoint 2 not, 2 = endpoint 2 above and endpoint 1 not,
// 3 = neither.
//
// With a fixed control (no control distributions), the future effect of
// endpoint k is S_tk / m_t - v_k for the fixed response rate v_k. It exceeds
// the null value q_k if S_tk / m_t - v_k - q_k > 1e-10, so that a future
// effect that differs from q_k only by rounding errors is not counted as
// larger.

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace {

// Joint probabilities of (S_1, S_2) for m future subjects, stored as
// out[s1 * (m + 1) + s2].
void margin_probs(const double *a, int m, std::vector<double> &out) {
  const int w = m + 1;
  std::fill(out.begin(), out.end(), 0.0);
  // lg[k][c] = log Gamma(a_k + c) - log Gamma(a_k) - log Gamma(c + 1)
  std::vector<double> lg(4 * static_cast<std::size_t>(w));
  double a_sum = 0.0;
  for (int k = 0; k < 4; ++k) {
    a_sum += a[k];
    const double base = R::lgammafn(a[k]);
    for (int c = 0; c <= m; ++c) {
      lg[k * static_cast<std::size_t>(w) + c] =
        R::lgammafn(a[k] + c) - base - R::lgammafn(c + 1.0);
    }
  }
  const double lconst = R::lgammafn(m + 1.0) + R::lgammafn(a_sum) -
    R::lgammafn(a_sum + m);
  for (int s1 = 0; s1 <= m; ++s1) {
    for (int s2 = 0; s2 <= m; ++s2) {
      const int lo = std::max(0, s1 + s2 - m);
      const int hi = std::min(s1, s2);
      double total = 0.0;
      for (int c11 = lo; c11 <= hi; ++c11) {
        const int c10 = s1 - c11;
        const int c01 = s2 - c11;
        const int c00 = m - c10 - c01 - c11;
        total += std::exp(lconst + lg[c00] + lg[w + c01] + lg[2 * w + c10] +
                          lg[3 * w + c11]);
      }
      out[static_cast<std::size_t>(s1) * w + s2] = total;
    }
  }
}

// K(s) for s = 0, ..., m_t: the largest integer k in [-1, m_c] with
// (m_c s - m_t k) / (m_t m_c) > q.
std::vector<int> cut_points(int m_t, int m_c, double q) {
  // The differences lie in [-1, 1], so values of q outside [-2, 2] give the
  // same cut points as -2 or 2 (and would lose precision below)
  q = std::min(std::max(q, -2.0), 2.0);
  std::vector<int> out(m_t + 1);
  const double den = static_cast<double>(m_t) * m_c;
  for (int s = 0; s <= m_t; ++s) {
    const double num = static_cast<double>(m_c) * s;
    const double k_est = std::floor((num - q * den) / m_t);
    double k_max = -1.0;
    for (int shift = -1; shift <= 1; ++shift) {
      const double cand = k_est + shift;
      if ((num - m_t * cand) / den > q && cand > k_max) {
        k_max = cand;
      }
    }
    k_max = std::min(std::max(k_max, -1.0), static_cast<double>(m_c));
    out[s] = static_cast<int>(k_max);
  }
  return out;
}

}  // namespace

// alpha_t, alpha_c: Dirichlet parameters, one row per distribution.
// pair_t, pair_c: 0-based row indices of the pairs to evaluate.
// m_t, m_c: numbers of future subjects; q: null values of the two endpoints.
// masks: one integer per union of regions; bit r (0 to 3) is set if region
//   r + 1 belongs to the union.
// value_c: marginal response rates of a fixed control, used when alpha_c has
//   no rows (pair_c and m_c are then not used).
// Returns the probabilities, one row per pair and one column per union.
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_pred_regions_binary2(Rcpp::NumericMatrix alpha_t,
                                             Rcpp::NumericMatrix alpha_c,
                                             Rcpp::IntegerVector pair_t,
                                             Rcpp::IntegerVector pair_c,
                                             int m_t, int m_c,
                                             Rcpp::NumericVector q,
                                             Rcpp::IntegerVector masks,
                                             Rcpp::NumericVector value_c) {
  const int n_t = alpha_t.nrow();
  const int n_c = alpha_c.nrow();
  const bool fixed = n_c == 0;
  const int n_pair = static_cast<int>(pair_t.size());
  const int n_mask = static_cast<int>(masks.size());
  const int wt = m_t + 1;
  const int wc = m_c + 2;  // control distribution function, index shift 1
  // Treatment probabilities of (S_t1, S_t2)
  std::vector<double> prob_t(static_cast<std::size_t>(n_t) * wt * wt);
  std::vector<double> buf_t(static_cast<std::size_t>(wt) * wt);
  double a[4];
  for (int i = 0; i < n_t; ++i) {
    for (int k = 0; k < 4; ++k) {
      a[k] = alpha_t(i, k);
    }
    margin_probs(a, m_t, buf_t);
    std::copy(buf_t.begin(), buf_t.end(),
              prob_t.begin() + static_cast<std::size_t>(i) * wt * wt);
  }
  // Control distribution functions F(x, y) = P(S_c1 <= x, S_c2 <= y) for
  // x, y = -1, ..., m_c, stored as cdf_c[(x + 1) * wc + (y + 1)]
  std::vector<double> cdf_c(static_cast<std::size_t>(n_c) * wc * wc, 0.0);
  std::vector<double> buf_c(static_cast<std::size_t>(m_c + 1) * (m_c + 1));
  for (int j = 0; j < n_c; ++j) {
    for (int k = 0; k < 4; ++k) {
      a[k] = alpha_c(j, k);
    }
    margin_probs(a, m_c, buf_c);
    double *f = &cdf_c[static_cast<std::size_t>(j) * wc * wc];
    for (int x = 0; x <= m_c; ++x) {
      double row = 0.0;
      for (int y = 0; y <= m_c; ++y) {
        row += buf_c[static_cast<std::size_t>(x) * (m_c + 1) + y];
        f[(x + 1) * wc + (y + 1)] = f[x * wc + (y + 1)] + row;
      }
    }
  }
  const std::vector<int> k1 = cut_points(m_t, m_c, q[0]);
  const std::vector<int> k2 = cut_points(m_t, m_c, q[1]);
  // Fixed control: whether S_tk / m_t - v_k exceeds q_k
  std::vector<int> up1(m_t + 1, 0);
  std::vector<int> up2(m_t + 1, 0);
  if (fixed) {
    for (int s = 0; s <= m_t; ++s) {
      const double rate = static_cast<double>(s) / m_t;
      up1[s] = (rate - value_c[0]) - q[0] > 1e-10;
      up2[s] = (rate - value_c[1]) - q[1] > 1e-10;
    }
  }
  Rcpp::NumericMatrix out(n_pair, n_mask);
  for (int p = 0; p < n_pair; ++p) {
    const double *pt = &prob_t[static_cast<std::size_t>(pair_t[p]) * wt * wt];
    double both = 0.0;
    double above1 = 0.0;
    double above2 = 0.0;
    if (fixed) {
      for (int s1 = 0; s1 <= m_t; ++s1) {
        for (int s2 = 0; s2 <= m_t; ++s2) {
          const double w = pt[s1 * wt + s2];
          both += w * up1[s1] * up2[s2];
          above1 += w * up1[s1];
          above2 += w * up2[s2];
        }
      }
    } else {
      const double *f =
        &cdf_c[static_cast<std::size_t>(pair_c[p]) * wc * wc];
      for (int s1 = 0; s1 <= m_t; ++s1) {
        const int x = k1[s1] + 1;
        for (int s2 = 0; s2 <= m_t; ++s2) {
          const double w = pt[s1 * wt + s2];
          const int y = k2[s2] + 1;
          both += w * f[x * wc + y];
          above1 += w * f[x * wc + (m_c + 1)];
          above2 += w * f[(m_c + 1) * wc + y];
        }
      }
    }
    double region[4];
    region[0] = both;
    region[1] = above1 - both;
    region[2] = above2 - both;
    region[3] = 1.0 - above1 - above2 + both;
    for (int m = 0; m < n_mask; ++m) {
      double total = 0.0;
      for (int r = 0; r < 4; ++r) {
        if ((masks[m] >> r) & 1) {
          total += region[r];
        }
      }
      out(p, m) = std::min(1.0, std::max(0.0, total));
    }
    if ((p + 1) % 10000 == 0) {
      Rcpp::checkUserInterrupt();
    }
  }
  return out;
}
