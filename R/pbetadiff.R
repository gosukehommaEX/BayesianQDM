#' Distribution function of the difference of two Beta variables
#'
#' Computes \eqn{\Pr(\pi_t - \pi_c \le q)} (or the upper tail
#' \eqn{\Pr(\pi_t - \pi_c > q)}) for independent
#' \eqn{\pi_t \sim \textrm{Beta}(a_t, b_t)} and
#' \eqn{\pi_c \sim \textrm{Beta}(a_c, b_c)}. This is the posterior probability
#' used for a single binary endpoint.
#'
#' The upper tail is computed as
#' \deqn{\Pr(\pi_t - \pi_c > q) = \Pr(\pi_c < -q) +
#'   \int f_c(x) \Pr(\pi_t > x + q) \, dx,}
#' where the integral runs over the values of \eqn{x} with
#' \eqn{0 \le x + q \le 1}, and the first term is zero for
#' \eqn{q \ge 0}. The integral is evaluated by adaptive Gauss-Kronrod
#' quadrature (the QUADPACK routine used by \code{\link[stats]{integrate}}),
#' with the range split at quantiles of both Beta distributions so that
#' concentrated distributions are integrated accurately. When the density of
#' \eqn{\pi_c} is unbounded at 0 (\eqn{a_c < 1}) or at 1 (\eqn{b_c < 1}), the
#' integral on that side of the median of \eqn{\pi_c} is evaluated in the
#' variable \eqn{s = x^{a_c}} or \eqn{s = (1 - x)^{b_c}}, in which the
#' integrand is bounded. The lower tail is
#' computed directly from the same formula with the roles of the two groups
#' exchanged, so both tails are accurate near 0.
#'
#' @param q Numeric vector of quantiles.
#' @param a_t,b_t Positive shape parameters of the treatment group.
#' @param a_c,b_c Positive shape parameters of the control group.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are
#'   \eqn{\Pr(\pi_t - \pi_c \le q)}, otherwise \eqn{\Pr(\pi_t - \pi_c > q)}.
#'
#' @return A numeric vector of probabilities. All arguments except
#'   \code{lower.tail} are recycled to a common length; each must have length
#'   1 or the common length.
#'
#' @family distributions
#' @examples
#' # Posterior probability that the response rate difference exceeds 0.1,
#' # with Jeffreys priors and 7/12 versus 4/12 responders
#' pbetadiff(0.1, a_t = 7.5, b_t = 5.5, a_c = 4.5, b_c = 8.5,
#'           lower.tail = FALSE)
#' @export
pbetadiff <- function(q, a_t, b_t, a_c, b_c, lower.tail = TRUE) {
  check_numeric(q, "q", len = NULL)
  check_numeric(a_t, "a_t", len = NULL, lower = 0, lower_open = TRUE)
  check_numeric(b_t, "b_t", len = NULL, lower = 0, lower_open = TRUE)
  check_numeric(a_c, "a_c", len = NULL, lower = 0, lower_open = TRUE)
  check_numeric(b_c, "b_c", len = NULL, lower = 0, lower_open = TRUE)
  check_flag(lower.tail, "lower.tail")
  args <- recycle_common(list(q = q, a_t = a_t, b_t = b_t, a_c = a_c,
                              b_c = b_c))
  if (lower.tail) {
    # P(pi_t - pi_c <= q) = P(pi_c - pi_t >= -q), the upper tail with the
    # groups exchanged (the distribution is continuous).
    cpp_pbetadiff_upper(-args$q, args$a_c, args$b_c, args$a_t, args$b_t)
  } else {
    cpp_pbetadiff_upper(args$q, args$a_t, args$b_t, args$a_c, args$b_c)
  }
}
