#' Cumulative Distribution Function of the Difference Between Two
#' Independent Beta-Binomial Proportions
#'
#' Calculates the cumulative distribution function (CDF) of the difference
#' between two independent Beta-Binomial proportions by exact enumeration.
#' Specifically, computes \eqn{P((Y_t / m_t) - (Y_c / m_c) \le q)} or
#' \eqn{P((Y_t / m_t) - (Y_c / m_c) > q)}, where
#' \eqn{Y_j \sim \mathrm{BetaBinomial}(m_j, \alpha_j, \beta_j)} for
#' \eqn{j \in \{t, c\}}.
#'
#' @param q A numeric scalar representing the quantile threshold for the
#'        difference in proportions.
#' @param m_t A positive integer giving the number of future patients in the
#'        treatment group.
#' @param m_c A positive integer giving the number of future patients in the
#'        control group.
#' @param alpha_t A positive numeric scalar giving the first shape parameter
#'        of the Beta mixing distribution for the treatment group.
#' @param alpha_c A positive numeric scalar giving the first shape parameter
#'        of the Beta mixing distribution for the control group.
#' @param beta_t A positive numeric scalar giving the second shape parameter
#'        of the Beta mixing distribution for the treatment group.
#' @param beta_c A positive numeric scalar giving the second shape parameter
#'        of the Beta mixing distribution for the control group.
#' @param lower.tail A logical scalar; if \code{TRUE} (default), the function
#'        returns \eqn{P((Y_t / m_t) - (Y_c / m_c) \le q)}, otherwise
#'        \eqn{P((Y_t / m_t) - (Y_c / m_c) > q)}.
#'
#' @return A numeric scalar in \code{[0, 1]}.
#'
#' @details
#' The probability mass function of
#' \eqn{Y_j \sim \mathrm{BetaBinomial}(m_j, \alpha_j, \beta_j)} is:
#' \deqn{P(Y_j = k) = \binom{m_j}{k}
#'   \frac{B(k + \alpha_j,\; m_j - k + \beta_j)}{B(\alpha_j, \beta_j)},
#'   \quad k = 0, \ldots, m_j}
#' where \eqn{B(\cdot, \cdot)} is the Beta function.
#'
#' The exact CDF is obtained by enumerating all
#' \eqn{(m_t + 1)(m_c + 1)} outcome combinations and summing the joint
#' probabilities for which the proportion difference satisfies the specified
#' condition. Computation time therefore grows quadratically in \eqn{m_t} and
#' \eqn{m_c}; for large future sample sizes consider a normal approximation.
#'
#' The Beta-Binomial distribution arises when the success probability in a
#' Binomial model follows a Beta prior, making it appropriate for
#' posterior predictive calculations in Bayesian binary-endpoint trials.
#'
#' @examples
#' # P((Y_t/12) - (Y_c/12) > 0.2) with symmetric Beta(0.5, 0.5) priors
#' pbetabinomdiff(0.2, 12, 12, 0.5, 0.5, 0.5, 0.5, lower.tail = FALSE)
#'
#' # P((Y_t/20) - (Y_c/15) > 0.1) with different future sample sizes
#' pbetabinomdiff(0.1, 20, 15, 1, 1, 1, 1, lower.tail = FALSE)
#'
#' # P((Y_t/10) - (Y_c/10) > 0) with informative priors
#' pbetabinomdiff(0, 10, 10, 2, 3, 3, 2, lower.tail = FALSE)
#'
#' # Lower tail: P((Y_t/15) - (Y_c/15) <= 0.05) with vague priors
#' pbetabinomdiff(0.05, 15, 15, 1, 1, 1, 1, lower.tail = TRUE)
#'
#' @export
pbetabinomdiff <- function(q, m_t, m_c, alpha_t, alpha_c, beta_t, beta_c,
                           lower.tail = TRUE) {

  # --- Input validation ---
  if (!is.numeric(q) || length(q) != 1L || is.na(q)) {
    stop("'q' must be a single numeric value")
  }

  if (!is.numeric(m_t) || length(m_t) != 1L || is.na(m_t) ||
      m_t != floor(m_t) || m_t < 1L) {
    stop("'m_t' must be a single positive integer")
  }

  if (!is.numeric(m_c) || length(m_c) != 1L || is.na(m_c) ||
      m_c != floor(m_c) || m_c < 1L) {
    stop("'m_c' must be a single positive integer")
  }

  if (!is.numeric(alpha_t) || length(alpha_t) != 1L || is.na(alpha_t) || alpha_t <= 0) {
    stop("'alpha_t' must be a single positive numeric value")
  }

  if (!is.numeric(alpha_c) || length(alpha_c) != 1L || is.na(alpha_c) || alpha_c <= 0) {
    stop("'alpha_c' must be a single positive numeric value")
  }

  if (!is.numeric(beta_t) || length(beta_t) != 1L || is.na(beta_t) || beta_t <= 0) {
    stop("'beta_t' must be a single positive numeric value")
  }

  if (!is.numeric(beta_c) || length(beta_c) != 1L || is.na(beta_c) || beta_c <= 0) {
    stop("'beta_c' must be a single positive numeric value")
  }

  if (!is.logical(lower.tail) || length(lower.tail) != 1L || is.na(lower.tail)) {
    stop("'lower.tail' must be a single logical value (TRUE or FALSE)")
  }

  # --- Compute Beta-Binomial PMFs ---
  # PMF for Y_t ~ BetaBinomial(m_t, alpha_t, beta_t)
  k_t <- 0:m_t
  dbetabinom_t <- choose(m_t, k_t) *
    beta(k_t + alpha_t, m_t - k_t + beta_t) / beta(alpha_t, beta_t)

  # PMF for Y_c ~ BetaBinomial(m_c, alpha_c, beta_c)
  k_c <- 0:m_c
  dbetabinom_c <- choose(m_c, k_c) *
    beta(k_c + alpha_c, m_c - k_c + beta_c) / beta(alpha_c, beta_c)

  # --- Enumerate all (m_t+1) x (m_c+1) outcome combinations ---
  # Proportion difference: (k_t/m_t) - (k_c/m_c) = (m_c*k_t - m_t*k_c) / (m_t*m_c)
  # Use integer arithmetic via outer() to avoid floating-point comparison issues.
  diff_numer <- outer(m_c * k_t, m_t * k_c, "-")  # numerator matrix (m_c*k_t - m_t*k_c)
  denom      <- m_t * m_c

  if (lower.tail) {
    I <- (diff_numer / denom <= q)
  } else {
    I <- (diff_numer / denom > q)
  }

  # Sum joint probabilities over outcome pairs satisfying the condition
  result <- as.numeric(crossprod(dbetabinom_t[row(I)[I]], dbetabinom_c[col(I)[I]]))

  return(result)
}
