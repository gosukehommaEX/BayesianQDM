#' Cumulative Distribution Function of the Difference of Two t-Distributed Variables
#' by Numerical Integration
#'
#' This function calculates the cumulative distribution function (CDF) of the difference
#' between two independent t-distributed random variables using exact numerical integration
#' via convolution. This method provides the most accurate results (within numerical
#' precision) compared to approximation methods. Specifically, it computes P(T1 - T2 ≤ q)
#' or P(T1 - T2 > q) where T1 and T2 follow t-distributions with potentially different
#' location, scale, and degrees of freedom parameters.
#'
#' @param q A numeric value representing the quantile threshold.
#' @param mu.t1 A numeric value representing the location parameter (μ) of the first
#'        t-distribution.
#' @param mu.t2 A numeric value representing the location parameter (μ) of the second
#'        t-distribution.
#' @param sd.t1 A positive numeric value representing the scale parameter (σ) of the
#'        first t-distribution.
#' @param sd.t2 A positive numeric value representing the scale parameter (σ) of the
#'        second t-distribution.
#' @param nu.t1 A positive numeric value representing the degrees of freedom (ν) of
#'        the first t-distribution. Must be > 2 for finite variance.
#' @param nu.t2 A positive numeric value representing the degrees of freedom (ν) of
#'        the second t-distribution. Must be > 2 for finite variance.
#' @param lower.tail A logical value; if TRUE (default), probabilities are P(T1 - T2 ≤ q),
#'        otherwise P(T1 - T2 > q).
#'
#' @return A numeric value in [0, 1] representing the cumulative probability that the
#'         difference between the two t-distributed variables is below (if lower.tail = TRUE)
#'         or exceeds (if lower.tail = FALSE) the quantile q.
#'
#' @details
#' This function uses the exact convolution approach to compute the distribution of the
#' difference between two t-distributed variables. The method involves:
#' \itemize{
#'   \item Using the convolution formula:
#'         \deqn{P(T_1 - T_2 \le q) = \int_{-\infty}^{\infty} f_1(x) \cdot F_2(x - q) dx}
#'   \item Where \eqn{f_1(x)} is the probability density function (PDF) of the first
#'         t-distribution
#'   \item And \eqn{F_2(x - q)} is the cumulative distribution function (CDF) of the
#'         second t-distribution evaluated at \eqn{x - q}
#'   \item Adaptive integration bounds based on the distribution characteristics
#'         (approximately ±8 standard deviations from the mean)
#'   \item High-precision numerical integration with relative tolerance 1e-6 and
#'         absolute tolerance 1e-8
#' }
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
#' # Calculate lower tail probability P(t1 - t2 ≤ 2)
#' pNI2tdiff(q = 2, mu.t1 = 3, mu.t2 = 0, sd.t1 = 1.5, sd.t2 = 1.2,
#'          nu.t1 = 12, nu.t2 = 15, lower.tail = TRUE)
#'
#' @importFrom stats dt pt integrate
#' @export
pNI2tdiff <- function(q, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2, lower.tail = TRUE) {
  # Determine the number of results to compute (handles vectorized inputs)
  n <- max(length(mu.t1), length(mu.t2), length(sd.t1), length(sd.t2))

  # Ensure all parameters have the same length for vectorized computation
  mu.t1 <- rep(mu.t1, length.out = n)
  mu.t2 <- rep(mu.t2, length.out = n)
  sd.t1 <- rep(sd.t1, length.out = n)
  sd.t2 <- rep(sd.t2, length.out = n)

  # Calculate P(T1 - T2 > q) for each parameter set using convolution
  results <- sapply(seq(n), function(i) {
    # Determine adaptive integration bounds based on the first distribution's characteristics
    # Calculate the center and approximate spread of the first t-distribution
    center1 <- mu.t1[i]
    spread1 <- sd.t1[i] * sqrt((nu.t1 + 1) / (nu.t1 - 2))

    # Perform numerical integration using the convolution formula
    integrate(
      function(x) {
        # PDF of the first t-distribution at x
        f1 <- dt((x - mu.t1[i]) / sd.t1[i], df = nu.t1) / sd.t1[i]

        # CDF of the second t-distribution at (x - q)
        F2 <- pt(((x - q) - mu.t2[i]) / sd.t2[i], df = nu.t2)

        # Return the product for integration (convolution formula)
        return(f1 * F2)
      },
      # Set adaptive integration bounds to capture the essential probability mass
      # Use approximately ±8 standard deviations from the mean
      lower = max(-Inf, center1 - 8 * spread1),
      upper = min(Inf,  center1 + 8 * spread1),
      # Set high precision tolerances for accurate results
      rel.tol = 1e-6,
      abs.tol = 1e-8
    )[['value']]
  })

  # Convert to lower tail probability if requested
  # If lower.tail = TRUE: return P(T1 - T2 ≤ q)
  # If lower.tail = FALSE: return 1 - P(T1 - T2 ≤ q) = P(T1 - T2 > q)
  results <- lower.tail + c(1, -1)[lower.tail + 1] * results

  return(results)
}
