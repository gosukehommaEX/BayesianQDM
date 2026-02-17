#' Cumulative Distribution Function of the Difference of Two t-Distributed Variables
#' by Numerical Integration
#'
#' This function calculates the cumulative distribution function (CDF) of the difference
#' between two independent t-distributed random variables using exact numerical integration
#' via convolution. This method provides the most accurate results (within numerical
#' precision) compared to approximation methods. Specifically, it computes P(T1 - T2 <= q)
#' or P(T1 - T2 > q) where T1 and T2 follow t-distributions with potentially different
#' location, scale, and degrees of freedom parameters.
#'
#' @param q A numeric value representing the quantile threshold.
#' @param mu.t1 A numeric value or vector representing the location parameter (mu) of the
#'        first t-distribution.
#' @param mu.t2 A numeric value or vector representing the location parameter (mu) of the
#'        second t-distribution.
#' @param sd.t1 A positive numeric value or vector representing the scale parameter (sigma)
#'        of the first t-distribution.
#' @param sd.t2 A positive numeric value or vector representing the scale parameter (sigma)
#'        of the second t-distribution.
#' @param nu.t1 A positive numeric value representing the degrees of freedom (nu) of
#'        the first t-distribution. Must be > 2 for finite variance.
#' @param nu.t2 A positive numeric value representing the degrees of freedom (nu) of
#'        the second t-distribution. Must be > 2 for finite variance.
#' @param lower.tail A logical value; if TRUE (default), probabilities are P(T1 - T2 <= q),
#'        otherwise P(T1 - T2 > q).
#'
#' @return A numeric value or vector in \code{[0, 1]} representing the cumulative probability
#'         that the difference between the two t-distributed variables is below
#'         (if lower.tail = TRUE) or exceeds (if lower.tail = FALSE) the quantile q.
#'         If mu.t1, mu.t2, sd.t1, or sd.t2 are vectors, returns a vector of the same length.
#'
#' @details
#' This function uses the exact convolution approach to compute the distribution of the
#' difference between two t-distributed variables. The method involves:
#' \itemize{
#'   \item Using the convolution formula:
#'         P(T_1 - T_2 > q) = integral from -infinity to infinity of f_1(x) * F_2(x - q) dx
#'   \item Where f_1(x) is the probability density function (PDF) of the first
#'         t-distribution
#'   \item And F_2(x - q) is the cumulative distribution function (CDF) of the
#'         second t-distribution evaluated at (x - q)
#'   \item Adaptive integration bounds based on the distribution characteristics
#'         (approximately plus/minus 8 standard deviations from the mean)
#'   \item High-precision numerical integration with relative tolerance 1e-6 and
#'         absolute tolerance 1e-8
#' }
#'
#' When mu.t1, mu.t2, sd.t1, or sd.t2 are vectors, mapply is used to apply the
#' scalar integration function across all parameter sets efficiently.
#'
#' **Advantages**:
#' \itemize{
#'   \item Provides exact results within numerical precision
#'   \item Handles arbitrary combinations of parameters
#'   \item No approximations required
#' }
#'
#' **Computational considerations**:
#' \itemize{
#'   \item More computationally intensive than approximation methods
#'         (e.g., Welch-Satterthwaite)
#'   \item Recommended for final analyses where accuracy is critical
#'   \item For exploratory analyses, consider using faster approximation methods
#' }
#'
#' @examples
#' # Calculate P(t1 - t2 > 3) for equal parameters
#' pNI2tdiff(q = 3, mu.t1 = 2, mu.t2 = 0, sd.t1 = 1, sd.t2 = 1,
#'          nu.t1 = 17, nu.t2 = 17, lower.tail = FALSE)
#'
#' # Calculate P(t1 - t2 > 1) for unequal variances
#' pNI2tdiff(q = 1, mu.t1 = 5, mu.t2 = 3, sd.t1 = 2, sd.t2 = 1.5,
#'          nu.t1 = 10, nu.t2 = 15, lower.tail = FALSE)
#'
#' # Calculate P(t1 - t2 > 0) for different degrees of freedom
#' pNI2tdiff(q = 0, mu.t1 = 1, mu.t2 = 1, sd.t1 = 1, sd.t2 = 1,
#'          nu.t1 = 5, nu.t2 = 20, lower.tail = FALSE)
#'
#' # Calculate lower tail probability P(t1 - t2 <= 2)
#' pNI2tdiff(q = 2, mu.t1 = 3, mu.t2 = 0, sd.t1 = 1.5, sd.t2 = 1.2,
#'          nu.t1 = 12, nu.t2 = 15, lower.tail = TRUE)
#'
#' # Vectorized usage: calculate for multiple parameter sets simultaneously
#' pNI2tdiff(q = 1, mu.t1 = c(2, 3, 4), mu.t2 = c(0, 1, 2),
#'          sd.t1 = c(1, 1.2, 1.5), sd.t2 = c(1, 1.1, 1.3),
#'          nu.t1 = 10, nu.t2 = 10, lower.tail = FALSE)
#'
#' @importFrom stats dt pt integrate
#' @export
pNI2tdiff <- function(q, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2, lower.tail = TRUE) {

  # Define scalar integration function for a single parameter set
  .integrate_one <- function(mu1, mu2, sd1, sd2) {
    # Adaptive integration bounds: center +/- 8 * spread of the first distribution
    spread1 <- sd1 * sqrt((nu.t1 + 1) / (nu.t1 - 2))

    integrate(
      function(x) {
        # PDF of the first t-distribution at x
        f1 <- dt((x - mu1) / sd1, df = nu.t1) / sd1

        # CDF of the second t-distribution evaluated at (x - q)
        F2 <- pt(((x - q) - mu2) / sd2, df = nu.t2)

        f1 * F2
      },
      lower = mu1 - 8 * spread1,
      upper = mu1 + 8 * spread1,
      rel.tol = 1e-6,
      abs.tol = 1e-8
    )[['value']]
  }

  # Apply integration across all parameter sets using mapply
  # For scalar inputs, mapply reduces to a single call (no overhead)
  results <- mapply(.integrate_one, mu1 = mu.t1, mu2 = mu.t2, sd1 = sd.t1, sd2 = sd.t2)

  # Convert to lower tail if requested:
  #   lower.tail = TRUE  -> return P(T1 - T2 <= q) = 1 - P(T1 - T2 > q)
  #   lower.tail = FALSE -> return P(T1 - T2 > q)
  results <- lower.tail + c(1, -1)[lower.tail + 1] * results

  return(results)
}
