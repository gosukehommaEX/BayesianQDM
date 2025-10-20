#' Cumulative Distribution Function of the Difference Between Two Beta Variables
#'
#' This function calculates the cumulative distribution function (CDF) of the difference
#' between two independent beta-distributed random variables using numerical integration
#' of the probability density function. The density function uses Appell's first
#' hypergeometric function for computation. Specifically, it computes P(X1 - X2 ≤ q) or
#' P(X1 - X2 > q) where X1 ~ Beta(alpha1, beta1) and X2 ~ Beta(alpha2, beta2).
#'
#' @param q A numeric value in the range \code{[-1, 1]} representing the quantile threshold.
#' @param alpha1 A positive numeric value representing the first shape parameter
#'        of the first beta distribution.
#' @param alpha2 A positive numeric value representing the first shape parameter
#'        of the second beta distribution.
#' @param beta1 A positive numeric value representing the second shape parameter
#'        of the first beta distribution.
#' @param beta2 A positive numeric value representing the second shape parameter
#'        of the second beta distribution.
#' @param lower.tail A logical value; if TRUE (default), probabilities are P(X1 - X2 ≤ q),
#'        otherwise P(X1 - X2 > q).
#'
#' @return A numeric value in \code{[0, 1]} representing the cumulative probability that
#'         the difference between the two beta variables is below (if lower.tail = TRUE)
#'         or exceeds (if lower.tail = FALSE) the quantile q.
#'
#' @details
#' The function computes the CDF by numerically integrating the probability density
#' function of the difference using adaptive quadrature. The density is evaluated
#' using piecewise formulations based on Appell's F1 hypergeometric function:
#' \itemize{
#'   \item For π ∈ \code{[-1, 0)}: Uses Appell's F1 with parameters adjusted for negative differences
#'   \item For π ∈ \code{[0, 1)}: Uses Appell's F1 with parameters adjusted for positive differences
#'   \item Special case handling when π = 0 and certain parameter conditions are met
#' }
#'
#' The upper tail probability P(X1 - X2 > q) is computed as:
#' \deqn{P(X1 - X2 > q) = \int_q^1 f(π) dπ}
#'
#' where f(π) is the probability density function of the difference.
#'
#' @examples
#' # Calculate P(Beta(0.5, 0.5) - Beta(0.5, 0.5) > 0.2)
#' p2betadiff(0.2, 0.5, 0.5, 0.5, 0.5, lower.tail = FALSE)
#'
#' # Calculate P(Beta(2, 3) - Beta(1, 4) > -0.1)
#' p2betadiff(-0.1, 2, 3, 1, 4, lower.tail = FALSE)
#'
#' # Calculate P(Beta(1, 1) - Beta(1, 1) > 0) - should be approximately 0.5
#' p2betadiff(0, 1, 1, 1, 1, lower.tail = FALSE)
#'
#' # Calculate lower tail probability P(Beta(2, 2) - Beta(2, 2) ≤ 0.1)
#' p2betadiff(0.1, 2, 2, 2, 2, lower.tail = TRUE)
#'
#' @importFrom stats integrate
#' @export
p2betadiff <- function(q, alpha1, alpha2, beta1, beta2, lower.tail = TRUE) {
  # Calculate the normalization constant using beta functions
  k <- 1 / (beta(alpha1, beta1) * beta(alpha2, beta2))

  # Calculate gamma parameter for Appell's function
  gamma <- alpha1 + alpha2 + beta1 + beta2 - 2

  # Compute the upper tail probability P(X1 - X2 > q) using numerical integration
  # Integrate the density function from q to 1
  results <- integrate(
    Vectorize(function(pi) {
      density <- d2betadiff(pi, alpha1, alpha2, beta1, beta2)
      return(density)
    }),
    lower = q,
    upper = 1
  )[['value']]

  # Convert to lower tail probability if requested
  # If lower.tail = TRUE: return 1 - P(X1 - X2 > q) = P(X1 - X2 ≤ q)
  # If lower.tail = FALSE: return P(X1 - X2 > q)
  results <- lower.tail + c(1, -1)[lower.tail + 1] * results

  return(results)
}
