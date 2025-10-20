#' Cumulative Distribution Function of the Difference Between Two Beta-Binomial Variables
#'
#' This function calculates the cumulative distribution function (CDF) of the difference
#' between two independent beta-binomial random variables. The computation uses exact
#' enumeration over all possible combinations of outcomes. Specifically, it computes
#' P((Y1/m1) - (Y2/m2) ≤ q) or P((Y1/m1) - (Y2/m2) > q) where Y1 ~ BetaBinomial(m1, alpha1, beta1)
#' and Y2 ~ BetaBinomial(m2, alpha2, beta2).
#'
#' @param q A numeric value representing the quantile threshold for the difference
#'        in proportions.
#' @param m1 A positive integer representing the number of trials for the first
#'        beta-binomial distribution.
#' @param m2 A positive integer representing the number of trials for the second
#'        beta-binomial distribution.
#' @param alpha1 A positive numeric value representing the first shape parameter
#'        of the first beta distribution.
#' @param alpha2 A positive numeric value representing the first shape parameter
#'        of the second beta distribution.
#' @param beta1 A positive numeric value representing the second shape parameter
#'        of the first beta distribution.
#' @param beta2 A positive numeric value representing the second shape parameter
#'        of the second beta distribution.
#' @param lower.tail A logical value; if TRUE (default), probabilities are
#'        P((Y1/m1) - (Y2/m2) ≤ q), otherwise P((Y1/m1) - (Y2/m2) > q).
#'
#' @return A numeric value in \code{[0, 1]} representing the cumulative probability that
#'         the difference in proportions between the two beta-binomial variables
#'         is below (if lower.tail = TRUE) or exceeds (if lower.tail = FALSE) the quantile q.
#'
#' @details
#' The function computes the exact probability by enumerating all possible combinations
#' of outcomes for both beta-binomial distributions and summing the joint probabilities
#' for which the proportion difference satisfies the specified condition.
#'
#' The probability mass function for BetaBinomial(m, α, β) is:
#' \deqn{P(Y = k) = \binom{m}{k} \frac{B(k + \alpha, m - k + \beta)}{B(\alpha, \beta)}}
#'
#' where B(·,·) is the beta function.
#'
#' The beta-binomial distribution arises when the success probability in a binomial
#' distribution follows a beta distribution, making it useful for modeling overdispersed
#' count data. For large values of m1 and m2, the computation may become slow as it
#' requires evaluating (m1 + 1) × (m2 + 1) combinations.
#'
#' @examples
#' # Calculate P((Y1/12) - (Y2/12) > 0.2) for symmetric beta priors
#' p2betabinomdiff(0.2, 12, 12, 0.5, 0.5, 0.5, 0.5, lower.tail = FALSE)
#'
#' # Calculate P((Y1/20) - (Y2/15) > 0.1) for different sample sizes
#' p2betabinomdiff(0.1, 20, 15, 1, 1, 1, 1, lower.tail = FALSE)
#'
#' # Calculate P((Y1/10) - (Y2/10) > 0) for informative priors
#' p2betabinomdiff(0, 10, 10, 2, 3, 3, 2, lower.tail = FALSE)
#'
#' # Calculate lower tail probability P((Y1/15) - (Y2/15) ≤ 0.05)
#' p2betabinomdiff(0.05, 15, 15, 1, 1, 1, 1, lower.tail = TRUE)
#'
#' @export
p2betabinomdiff <- function(q, m1, m2, alpha1, alpha2, beta1, beta2, lower.tail = TRUE) {
  # Calculate probability mass functions for both beta-binomial distributions
  # PMF for first beta-binomial: Y1 ~ BetaBinomial(m1, alpha1, beta1)
  dbetabinom1 <- choose(m1, 0:m1) * beta(0:m1 + alpha1, m1 - (0:m1) + beta1) / beta(alpha1, beta1)

  # PMF for second beta-binomial: Y2 ~ BetaBinomial(m2, alpha2, beta2)
  dbetabinom2 <- choose(m2, 0:m2) * beta(0:m2 + alpha2, m2 - (0:m2) + beta2) / beta(alpha2, beta2)

  # Create indicator matrix for the condition (Y1/m1) - (Y2/m2) compared to q
  # This matrix compares all possible combinations of outcomes
  if(lower.tail) {
    # For lower tail: P((Y1/m1) - (Y2/m2) ≤ q)
    I <- (outer(m2 * 0:m1, m1 * 0:m2, '-') / (m1 * m2) <= q)
  } else {
    # For upper tail: P((Y1/m1) - (Y2/m2) > q)
    I <- (outer(m2 * 0:m1, m1 * 0:m2, '-') / (m1 * m2) > q)
  }

  # Calculate the probability by summing over all valid combinations
  # Use efficient matrix operations: multiply PMFs for outcomes that satisfy the condition
  results <- as.numeric(crossprod(dbetabinom1[row(I)[I]], dbetabinom2[col(I)[I]]))

  return(results)
}
