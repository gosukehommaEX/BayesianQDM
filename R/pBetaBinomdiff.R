#' Cumulative distribution function (CDF) of the difference between two beta-binomial variables
#'
#' This function calculates CDF(s) of difference between two beta-binomial variables.
#'
#' @param q Quantile value.
#' @param m1 A number of patients in group 1 for the future trial data.
#' @param m2 A number of patients in group 2 for the future trial data.
#' @param alpha1 The first non-negative numeric value of beta function for group 1.
#' @param alpha2 The first non-negative numeric value of beta function for group 2.
#' @param beta1 The second non-negative numeric value of beta function for group 1.
#' @param beta2 The second non-negative numeric value of beta function for group 2.
#'
#' @return The \code{pBetaBinomdiff} gives the numeric value of the density.
#'
#' @examples
#' pBetaBinomdiff(0.2, 12, 12, 0.5, 0.5, 0.5, 0.5)
#'
#' @export
pBetaBinomdiff <- function(q, m1, m2, alpha1, alpha2, beta1, beta2) {
  # Probability mass functions of beta-binomial distribution
  dbetabinom1 <- choose(m1, 0:m1) * beta(0:m1 + alpha1, m1 - (0:m1) + beta1) / beta(alpha1, beta1)
  dbetabinom2 <- choose(m2, 0:m2) * beta(0:m2 + alpha2, m2 - (0:m2) + beta2) / beta(alpha2, beta2)
  # A posterior predictive probability
  I <- (outer(m2 * 0:m1, m1 * 0:m2, '-') / (m1 * m2) > q)
  results <- as.numeric(crossprod(dbetabinom1[row(I)[I]], dbetabinom2[col(I)[I]]))
  return(results)
}
