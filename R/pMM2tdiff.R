#' Cumulative Distribution Function of the Difference of Two t-Distributed Variables
#' Using the Moment-Matching Approximation
#'
#' This function calculates the cumulative distribution function (CDF) of the difference
#' between two independent t-distributed random variables using the Moment-Matching Approximation.
#' This method provides a fast and reasonably accurate approximation,
#' particularly useful for unequal variances and degrees of freedom. Specifically,
#' it computes P(T1 - T2 <= q) or P(T1 - T2 > q) where T1 and T2 follow t-distributions
#' with potentially different parameters.
#'
#' @param q A numeric value representing the quantile threshold.
#' @param mu.t1 A numeric value representing the location parameter of the first
#'        t-distribution.
#' @param mu.t2 A numeric value representing the location parameter of the second
#'        t-distribution.
#' @param sd.t1 A positive numeric value representing the scale parameter of the
#'        first t-distribution.
#' @param sd.t2 A positive numeric value representing the scale parameter of the
#'        second t-distribution.
#' @param nu.t1 A positive numeric value representing the degrees of freedom of
#'        the first t-distribution. Must be > 2 for finite variance.
#' @param nu.t2 A positive numeric value representing the degrees of freedom of
#'        the second t-distribution. Must be > 2 for finite variance.
#' @param lower.tail A logical value; if TRUE (default), probabilities are P(T1 - T2 <= q),
#'        otherwise P(T1 - T2 > q).
#'
#' @return A numeric value in \code{[0, 1]} representing the cumulative probability that the
#'         difference between the two t-distributed variables is below (if lower.tail = TRUE)
#'         or exceeds (if lower.tail = FALSE) the quantile q.
#'
#' @details
#' For two independent non-standardized t-distributed random variables:
#' \itemize{
#'   \item X1 ~ t(mu.t1, sd.t1^2, nu.t1)
#'   \item X2 ~ t(mu.t2, sd.t2^2, nu.t2)
#' }
#' The difference Z = X1 - X2 is approximated as: Z ~ t(mu.t1 - mu.t2, sigma_star^2, nu_star)
#' where the effective parameters are computed through moment matching:
#' \itemize{
#'   \item sigma_star is derived from the second moment matching
#'   \item nu_star is derived from the fourth moment matching
#' }
#' The method requires nu.t1 > 4 and nu.t2 > 4 for the existence of fourth moments. The approximation
#' quality improves as degrees of freedom increase and approaches exactness as
#' nu approaches infinity (normal limit).
#'
#' **Advantages**:
#' \itemize{
#'   \item Computationally fast and efficient
#'   \item Works well for unequal variances and degrees of freedom
#'   \item Provides good approximation in most practical scenarios
#'   \item Particularly accurate when degrees of freedom are not too small
#' }
#'
#' **Computational considerations**:
#' \itemize{
#'   \item Much faster than numerical integration methods
#'   \item Suitable for exploratory analyses and large-scale simulations
#'   \item For final analyses requiring highest accuracy, consider using
#'         numerical integration (pNI2tdiff)
#' }
#'
#' @examples
#' # Calculate P(t1 - t2 > 3) for equal parameters
#' pMM2tdiff(q = 3, mu.t1 = 2, mu.t2 = 0, sd.t1 = 1, sd.t2 = 1,
#'          nu.t1 = 17, nu.t2 = 17, lower.tail = FALSE)
#'
#' # Calculate P(t1 - t2 > 1) for unequal variances
#' pMM2tdiff(q = 1, mu.t1 = 5, mu.t2 = 3, sd.t1 = 2, sd.t2 = 1.5,
#'          nu.t1 = 10, nu.t2 = 15, lower.tail = FALSE)
#'
#' # Calculate P(t1 - t2 > 0) for different degrees of freedom
#' pMM2tdiff(q = 0, mu.t1 = 1, mu.t2 = 1, sd.t1 = 1, sd.t2 = 1,
#'          nu.t1 = 5, nu.t2 = 20, lower.tail = FALSE)
#'
#' # Calculate lower tail probability P(t1 - t2 <= 2)
#' pMM2tdiff(q = 2, mu.t1 = 3, mu.t2 = 0, sd.t1 = 1.5, sd.t2 = 1.2,
#'          nu.t1 = 12, nu.t2 = 15, lower.tail = TRUE)
#'
#' # Vectorized usage: calculate for multiple parameter sets simultaneously
#' pMM2tdiff(q = 1, mu.t1 = c(2, 3, 4), mu.t2 = c(0, 1, 2),
#'          sd.t1 = c(1, 1.2, 1.5), sd.t2 = c(1, 1.1, 1.3),
#'          nu.t1 = 10, nu.t2 = 10, lower.tail = FALSE)
#'
#' @importFrom stats pt
#' @export
pMM2tdiff <- function(q, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2, lower.tail = TRUE) {

  # Step 1: Calculate mean difference
  mu.diff <- mu.t1 - mu.t2

  # Step 2: Obtain sigma^{star}
  var1 <- sd.t1 ^ 2 * nu.t1 / (nu.t1 - 2)
  var2 <- sd.t2 ^ 2 * nu.t2 / (nu.t2 - 2)
  Qu_star <- (var1 + var2) ^ 2
  term1 <- (sd.t1 ^ 2) ^ 2 * nu.t1 ^ 2 / ((nu.t1 - 2) * (nu.t1 - 4))
  term2 <- (sd.t2 ^ 2) ^ 2 * nu.t2 ^ 2 / ((nu.t2 - 2) * (nu.t2 - 4))
  term3 <- 2 * var1 * var2
  Qu <- term1 + term2 + term3
  nu.star <- (2 * Qu_star - 4 * Qu) / (Qu_star - Qu)
  sigma.star <- sqrt((var1 + var2) * (nu.star - 2) / nu.star)

  # Step 3: Standardized statistic using moment-matching approximation
  standardized <- (q - mu.diff) / sigma.star
  results <- pt(standardized, df = nu.star, lower.tail = FALSE)

  # Convert to lower tail probability if requested
  # If lower.tail = TRUE: return 1 - P(T1 - T2 > q) = P(T1 - T2 <= q)
  # If lower.tail = FALSE: return P(T1 - T2 > q)
  results <- lower.tail + c(1, -1)[lower.tail + 1] * results

  return(results)
}
