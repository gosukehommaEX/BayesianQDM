#' Distribution function of the difference of two Beta-binomial proportions
#'
#' Computes \eqn{\Pr(Y_t / m_t - Y_c / m_c \le q)} (or the upper tail
#' \eqn{\Pr(Y_t / m_t - Y_c / m_c > q)}) for independent
#' \eqn{Y_t \sim \textrm{Beta-binomial}(m_t, a_t, b_t)} and
#' \eqn{Y_c \sim \textrm{Beta-binomial}(m_c, a_c, b_c)}. This is the posterior
#' predictive probability for the difference in response proportions of a
#' future study with \eqn{m_t} and \eqn{m_c} subjects.
#'
#' The probability is an exact finite sum. The probability mass functions are
#' computed on the log scale, so large values of \eqn{m_t} and \eqn{m_c} do not
#' cause overflow. Each difference is compared with \code{q} after a single
#' division of the exact integer numerator \eqn{m_c Y_t - m_t Y_c} by
#' \eqn{m_t m_c}, so differences that are mathematically equal to a decimal
#' \code{q} are classified correctly (they belong to the lower tail).
#'
#' @param q Numeric vector of quantiles.
#' @param m_t,m_c Numbers of future subjects in the treatment and control
#'   groups (whole numbers of at least 1).
#' @param a_t,b_t Positive shape parameters of the treatment group.
#' @param a_c,b_c Positive shape parameters of the control group.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are
#'   \eqn{\Pr(Y_t / m_t - Y_c / m_c \le q)}, otherwise
#'   \eqn{\Pr(Y_t / m_t - Y_c / m_c > q)}.
#'
#' @return A numeric vector of probabilities. The arguments \code{q},
#'   \code{a_t}, \code{b_t}, \code{a_c} and \code{b_c} are recycled to a common
#'   length; each must have length 1 or the common length.
#'
#' @family distributions
#' @examples
#' # Predictive probability that the future difference exceeds 0, with
#' # Jeffreys priors, 7/12 versus 4/12 responders and 20 future subjects per
#' # group
#' pbetabinomdiff(0, m_t = 20, a_t = 7.5, b_t = 5.5, m_c = 20, a_c = 4.5,
#'                b_c = 8.5, lower.tail = FALSE)
#' @export
pbetabinomdiff <- function(q, m_t, a_t, b_t, m_c, a_c, b_c,
                           lower.tail = TRUE) {
  check_numeric(q, "q", len = NULL)
  check_numeric(m_t, "m_t", lower = 1, integer = TRUE)
  check_numeric(m_c, "m_c", lower = 1, integer = TRUE)
  check_numeric(a_t, "a_t", len = NULL, lower = 0, lower_open = TRUE)
  check_numeric(b_t, "b_t", len = NULL, lower = 0, lower_open = TRUE)
  check_numeric(a_c, "a_c", len = NULL, lower = 0, lower_open = TRUE)
  check_numeric(b_c, "b_c", len = NULL, lower = 0, lower_open = TRUE)
  check_flag(lower.tail, "lower.tail")
  args <- recycle_common(list(q = q, a_t = a_t, b_t = b_t, a_c = a_c,
                              b_c = b_c))
  m_t <- round(m_t)
  m_c <- round(m_c)
  k_t <- 0:m_t
  k_c <- 0:m_c
  # Numerators m_c * k_t of the differences on the common denominator
  # m_t * m_c (exact in double precision)
  num_t <- m_c * k_t
  den <- m_t * m_c
  out <- numeric(length(args$q))
  for (i in seq_along(out)) {
    q_i <- args$q[i]
    # For each k_t, the largest k_c with (m_c * k_t - m_t * k_c) / den > q_i.
    # The difference decreases in k_c, so the estimate below is within one of
    # the exact value; the three neighbouring candidates are checked with the
    # exact comparison.
    k_est <- floor((num_t - q_i * den) / m_t)
    k_max <- rep(-1, length(k_t))
    for (shift in c(-1, 0, 1)) {
      cand <- k_est + shift
      ok <- (num_t - m_t * cand) / den > q_i
      k_max <- ifelse(ok & cand > k_max, cand, k_max)
    }
    k_max <- pmin(pmax(k_max, -1), m_c)
    p_t <- dbetabinom(k_t, m_t, args$a_t[i], args$b_t[i])
    p_c <- dbetabinom(k_c, m_c, args$a_c[i], args$b_c[i])
    if (lower.tail) {
      # P(Y_c >= k_max + 1), with P(Y_c >= m_c + 1) = 0
      tail_c <- c(rev(cumsum(rev(p_c))), 0)
    } else {
      # P(Y_c <= k_max), with P(Y_c <= -1) = 0
      tail_c <- c(0, cumsum(p_c))
    }
    out[i] <- sum(p_t * tail_c[k_max + 2])
  }
  pmin(pmax(out, 0), 1)
}
