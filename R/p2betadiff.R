#' Cumulative Distribution Function of the Difference Between Two Beta Variables
#'
#' This function calculates the cumulative distribution function (CDF) of the difference
#' between two independent beta-distributed random variables. Specifically, it computes
#' P(X1 - X2 <= q) or P(X1 - X2 > q) where X1 ~ Beta(alpha1, beta1) and
#' X2 ~ Beta(alpha2, beta2).
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
#' @param lower.tail A logical value; if TRUE (default), probabilities are P(X1 - X2 <= q),
#'        otherwise P(X1 - X2 > q).
#'
#' @return A numeric value in \code{[0, 1]} representing the cumulative probability that
#'         the difference between the two beta variables is below (if lower.tail = TRUE)
#'         or exceeds (if lower.tail = FALSE) the quantile q.
#'
#' @details
#' The upper tail probability P(X1 - X2 > q) is expressed as the convolution:
#'
#' P(X1 - X2 > q) = integral_0^1 P(X2 < x - q) * f1(x) dx
#'                = integral_0^1 pbeta(x - q, alpha2, beta2) * dbeta(x, alpha1, beta1) dx
#'
#' where f1(x) = dbeta(x, alpha1, beta1) is the density of X1, and
#' pbeta(x - q, alpha2, beta2) = P(X2 < x - q) is the CDF of X2 evaluated at x - q.
#' The argument x - q is implicitly clipped to [0, 1] by pbeta (returns 0 for negative
#' values and 1 for values exceeding 1).
#'
#' This formulation replaces the original double integration structure
#' (outer integral over the difference, inner integral for Appell's F1 function)
#' with a single numerical integration. The two approaches are mathematically
#' equivalent by the law of total probability. The integrand uses only
#' R's native pbeta and dbeta functions, which are implemented in C and are
#' substantially faster than the recursive numerical integration required by
#' Appell's F1 hypergeometric function.
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
#' # Calculate lower tail probability P(Beta(2, 2) - Beta(2, 2) <= 0.1)
#' p2betadiff(0.1, 2, 2, 2, 2, lower.tail = TRUE)
#'
#' @importFrom stats integrate dbeta pbeta
#' @export
p2betadiff <- function(q, alpha1, alpha2, beta1, beta2, lower.tail = TRUE) {

  # Compute upper tail probability P(X1 - X2 > q) via convolution:
  #   P(X1 - X2 > q) = integral_0^1 P(X2 < x - q) * f1(x) dx
  #                  = integral_0^1 pbeta(x - q, alpha2, beta2) * dbeta(x, alpha1, beta1) dx
  #
  # pbeta handles boundary cases automatically:
  #   pbeta(x - q, ...) = 0 when x - q <= 0  (i.e., x <= q)
  #   pbeta(x - q, ...) = 1 when x - q >= 1  (i.e., x >= 1 + q)
  # So the effective integration range is (q, 1), but integrating over [0, 1]
  # is safe and avoids manual boundary handling.
  results <- integrate(
    function(x) pbeta(x - q, alpha2, beta2) * dbeta(x, alpha1, beta1),
    lower = 0,
    upper = 1
  )[['value']]

  # Convert to lower tail if requested:
  #   lower.tail = TRUE  -> return P(X1 - X2 <= q) = 1 - P(X1 - X2 > q)
  #   lower.tail = FALSE -> return P(X1 - X2 > q)
  results <- lower.tail + c(1, -1)[lower.tail + 1] * results

  return(results)
}
