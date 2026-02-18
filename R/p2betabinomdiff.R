#' Cumulative Distribution Function of the Difference Between Two Beta-Binomial Proportions
#'
#' Calculates the cumulative distribution function (CDF) of the difference
#' between two independent Beta-Binomial proportions by exact enumeration.
#' Specifically, computes \eqn{P((Y_1/m_1) - (Y_2/m_2) \le q)} or
#' \eqn{P((Y_1/m_1) - (Y_2/m_2) > q)}, where
#' \eqn{Y_1 \sim \mathrm{BetaBinomial}(m_1, \alpha_1, \beta_1)} and
#' \eqn{Y_2 \sim \mathrm{BetaBinomial}(m_2, \alpha_2, \beta_2)}.
#'
#' @param q A numeric scalar representing the quantile threshold for the
#'        difference in proportions.
#' @param m1 A positive integer giving the number of trials for the first
#'        Beta-Binomial distribution.
#' @param m2 A positive integer giving the number of trials for the second
#'        Beta-Binomial distribution.
#' @param alpha1 A positive numeric scalar giving the first shape parameter of
#'        the first Beta mixing distribution.
#' @param alpha2 A positive numeric scalar giving the first shape parameter of
#'        the second Beta mixing distribution.
#' @param beta1 A positive numeric scalar giving the second shape parameter of
#'        the first Beta mixing distribution.
#' @param beta2 A positive numeric scalar giving the second shape parameter of
#'        the second Beta mixing distribution.
#' @param lower.tail A logical scalar; if \code{TRUE} (default), the function
#'        returns \eqn{P((Y_1/m_1) - (Y_2/m_2) \le q)}, otherwise
#'        \eqn{P((Y_1/m_1) - (Y_2/m_2) > q)}.
#'
#' @return A numeric scalar in \code{[0, 1]}.
#'
#' @details
#' The probability mass function of \eqn{Y \sim \mathrm{BetaBinomial}(m, \alpha, \beta)} is:
#' \deqn{P(Y = k) = \binom{m}{k}
#'   \frac{B(k + \alpha,\; m - k + \beta)}{B(\alpha, \beta)}, \quad k = 0, \ldots, m}
#' where \eqn{B(\cdot, \cdot)} is the Beta function.
#'
#' The exact CDF is obtained by enumerating all \eqn{(m_1 + 1)(m_2 + 1)} outcome
#' combinations and summing the joint probabilities for which the proportion
#' difference satisfies the specified condition. Computation time therefore grows
#' quadratically in \eqn{m_1} and \eqn{m_2}; for large trial sizes consider
#' a normal approximation.
#'
#' The Beta-Binomial distribution arises when the success probability in a Binomial
#' model follows a Beta prior, making it appropriate for overdispersed count data
#' and for posterior predictive calculations in Bayesian binary-endpoint trials.
#'
#' @examples
#' # P((Y1/12) - (Y2/12) > 0.2) with symmetric Beta(0.5, 0.5) priors
#' p2betabinomdiff(0.2, 12, 12, 0.5, 0.5, 0.5, 0.5, lower.tail = FALSE)
#'
#' # P((Y1/20) - (Y2/15) > 0.1) with different sample sizes
#' p2betabinomdiff(0.1, 20, 15, 1, 1, 1, 1, lower.tail = FALSE)
#'
#' # P((Y1/10) - (Y2/10) > 0) with informative priors
#' p2betabinomdiff(0, 10, 10, 2, 3, 3, 2, lower.tail = FALSE)
#'
#' # Lower tail: P((Y1/15) - (Y2/15) <= 0.05)
#' p2betabinomdiff(0.05, 15, 15, 1, 1, 1, 1, lower.tail = TRUE)
#'
#' @export
p2betabinomdiff <- function(q, m1, m2, alpha1, alpha2, beta1, beta2,
                            lower.tail = TRUE) {

  # --- Input validation ---
  if (!is.numeric(q) || length(q) != 1L || is.na(q)) {
    stop("'q' must be a single numeric value")
  }

  for (arg in list(m1, m2)) {
    nm <- deparse(substitute(arg))
    if (!is.numeric(arg) || length(arg) != 1L || is.na(arg) ||
        arg != floor(arg) || arg < 1L) {
      stop(paste0("'", nm, "' must be a single positive integer"))
    }
  }

  if (!is.numeric(alpha1) || length(alpha1) != 1L || is.na(alpha1) || alpha1 <= 0) {
    stop("'alpha1' must be a single positive numeric value")
  }

  if (!is.numeric(alpha2) || length(alpha2) != 1L || is.na(alpha2) || alpha2 <= 0) {
    stop("'alpha2' must be a single positive numeric value")
  }

  if (!is.numeric(beta1) || length(beta1) != 1L || is.na(beta1) || beta1 <= 0) {
    stop("'beta1' must be a single positive numeric value")
  }

  if (!is.numeric(beta2) || length(beta2) != 1L || is.na(beta2) || beta2 <= 0) {
    stop("'beta2' must be a single positive numeric value")
  }

  if (!is.logical(lower.tail) || length(lower.tail) != 1L || is.na(lower.tail)) {
    stop("'lower.tail' must be a single logical value (TRUE or FALSE)")
  }

  # --- Compute Beta-Binomial PMFs ---
  # PMF for Y1 ~ BetaBinomial(m1, alpha1, beta1)
  k1 <- 0:m1
  dbetabinom1 <- choose(m1, k1) * beta(k1 + alpha1, m1 - k1 + beta1) / beta(alpha1, beta1)

  # PMF for Y2 ~ BetaBinomial(m2, alpha2, beta2)
  k2 <- 0:m2
  dbetabinom2 <- choose(m2, k2) * beta(k2 + alpha2, m2 - k2 + beta2) / beta(alpha2, beta2)

  # --- Enumerate all (m1+1) x (m2+1) outcome combinations ---
  # Proportion difference: (k1/m1) - (k2/m2) = (m2*k1 - m1*k2) / (m1*m2)
  # Use integer arithmetic via outer() to avoid floating-point comparison issues.
  diff_numer <- outer(m2 * k1, m1 * k2, "-")  # numerator matrix (m2*k1 - m1*k2)
  denom      <- m1 * m2

  if (lower.tail) {
    I <- (diff_numer / denom <= q)
  } else {
    I <- (diff_numer / denom > q)
  }

  # Sum joint probabilities over outcome pairs satisfying the condition
  result <- as.numeric(crossprod(dbetabinom1[row(I)[I]], dbetabinom2[col(I)[I]]))

  return(result)
}
