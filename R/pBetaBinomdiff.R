#' Cumulative Distribution Function of the Difference Between Two Beta-Binomial Variables
#'
#' This function calculates the cumulative distribution function (CDF) of the difference
#' between two independent beta-binomial random variables. Specifically, it computes
#' P((Y1/m1) - (Y2/m2) > q) where Y1 ~ BetaBinomial(m1, alpha1, beta1) and
#' Y2 ~ BetaBinomial(m2, alpha2, beta2).
#'
#' @param q A numeric value representing the quantile threshold for the proportion difference.
#' @param m1 A positive integer representing the number of trials for the first beta-binomial distribution.
#' @param m2 A positive integer representing the number of trials for the second beta-binomial distribution.
#' @param alpha1 A positive numeric value representing the first shape parameter of the first beta distribution.
#' @param alpha2 A positive numeric value representing the first shape parameter of the second beta distribution.
#' @param beta1 A positive numeric value representing the second shape parameter of the first beta distribution.
#' @param beta2 A positive numeric value representing the second shape parameter of the second beta distribution.
#'
#' @return A numeric value representing P((Y1/m1) - (Y2/m2) > q), the probability that the
#'         difference in proportions between the two beta-binomial variables exceeds the quantile q.
#'
#' @details
#' The function uses the probability mass functions of beta-binomial distributions to compute
#' the exact probability. The beta-binomial distribution arises when the success probability
#' in a binomial distribution follows a beta distribution, making it useful for modeling
#' overdispersed count data.
#'
#' The probability mass function for BetaBinomial(m, α, β) is:
#' \deqn{P(Y = k) = \binom{m}{k} \frac{B(k + \alpha, m - k + \beta)}{B(\alpha, \beta)}}
#'
#' where B(·,·) is the beta function.
#'
#' @examples
#' # Calculate P((Y1/12) - (Y2/12) > 0.2) for symmetric beta priors
#' pBetaBinomdiff(0.2, 12, 12, 0.5, 0.5, 0.5, 0.5)
#'
#' # Calculate P((Y1/20) - (Y2/15) > 0.1) for different sample sizes
#' pBetaBinomdiff(0.1, 20, 15, 1, 1, 1, 1)
#'
#' # Calculate P((Y1/10) - (Y2/10) > 0) for informative priors
#' pBetaBinomdiff(0, 10, 10, 2, 3, 3, 2)
#'
#' @export
pBetaBinomdiff <- function(q, m1, m2, alpha1, alpha2, beta1, beta2) {
  # Calculate probability mass functions of beta-binomial distributions
  # PMF for first beta-binomial: BetaBinomial(m1, alpha1, beta1)
  dbetabinom1 <- choose(m1, 0:m1) * beta(0:m1 + alpha1, m1 - (0:m1) + beta1) / beta(alpha1, beta1)

  # PMF for second beta-binomial: BetaBinomial(m2, alpha2, beta2)
  dbetabinom2 <- choose(m2, 0:m2) * beta(0:m2 + alpha2, m2 - (0:m2) + beta2) / beta(alpha2, beta2)

  # Create indicator matrix for the condition (Y1/m1) - (Y2/m2) > q
  # This compares all possible combinations of outcomes
  I <- (outer(m2 * 0:m1, m1 * 0:m2, '-') / (m1 * m2) > q)

  # Calculate the posterior predictive probability by summing over all valid combinations
  # Use matrix operations for efficient computation
  results <- as.numeric(crossprod(dbetabinom1[row(I)[I]], dbetabinom2[col(I)[I]]))

  return(results)
}
