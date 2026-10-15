#' Upper tail of the difference of two t variables by moment matching
#'
#' Approximates \eqn{\Pr(T_t - T_c > q)} by a single scaled t distribution
#' whose second and fourth moments match those of the difference (Yamaguchi et
#' al., submitted, Theorem 1). The matched degrees of freedom are
#' \deqn{\nu^* = 4 + \frac{(v_t + v_c)^2}{v_t^2 / (\nu_t - 4) +
#'   v_c^2 / (\nu_c - 4)},}
#' where \eqn{v_k = s_k^2 \nu_k / (\nu_k - 2)} is the variance of group
#' \eqn{k}, and the matched squared scale is
#' \eqn{(v_t + v_c)(\nu^* - 2) / \nu^*}. Infinite degrees of freedom (normal
#' distributions) are allowed.
#'
#' @param q,mu_t,scale_t,df_t,mu_c,scale_c,df_c Numeric vectors of a common
#'   length; all degrees of freedom must be greater than 4.
#'
#' @return A numeric vector of upper tail probabilities.
#' @noRd
ptdiff_mm_upper <- function(q, mu_t, scale_t, df_t, mu_c, scale_c, df_c) {
  var_t <- ifelse(is.finite(df_t), scale_t ^ 2 * df_t / (df_t - 2),
                  scale_t ^ 2)
  var_c <- ifelse(is.finite(df_c), scale_c ^ 2 * df_c / (df_c - 2),
                  scale_c ^ 2)
  kurt_t <- ifelse(is.finite(df_t), var_t ^ 2 / (df_t - 4), 0)
  kurt_c <- ifelse(is.finite(df_c), var_c ^ 2 / (df_c - 4), 0)
  kurt <- kurt_t + kurt_c
  nu_star <- ifelse(kurt > 0, 4 + (var_t + var_c) ^ 2 / kurt, Inf)
  scale_star <- sqrt((var_t + var_c) *
                       ifelse(is.finite(nu_star), (nu_star - 2) / nu_star, 1))
  stats::pt((q - (mu_t - mu_c)) / scale_star, df = nu_star,
            lower.tail = FALSE)
}
