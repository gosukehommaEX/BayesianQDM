#' Cumulative Distribution Function of the Difference of Two t-Distributed Variables by Numerical Integration
#'
#' This function calculates the cumulative distribution function (CDF) of the difference
#' between two independent t-distributed random variables using numerical integration
#' via convolution. Specifically, it computes P(T1 - T2 ≥ q) where T1 and T2 follow
#' t-distributions with potentially different location, scale, and degrees of freedom parameters.
#'
#' @param q A numeric value representing the quantile threshold.
#' @param mu.t1 A numeric value representing the location parameter of the first t-distribution.
#' @param mu.t2 A numeric value representing the location parameter of the second t-distribution.
#' @param sd.t1 A positive numeric value representing the scale parameter of the first t-distribution.
#' @param sd.t2 A positive numeric value representing the scale parameter of the second t-distribution.
#' @param nu.t1 A positive numeric value representing the degrees of freedom of the first t-distribution.
#' @param nu.t2 A positive numeric value representing the degrees of freedom of the second t-distribution.
#'
#' @return A numeric vector representing P(T1 - T2 ≥ q), the probability that the difference
#'         between the two t-distributed variables is greater than or equal to the quantile q.
#'
#' @details
#' This function uses the exact convolution approach to compute the distribution of the
#' difference between two t-distributed variables. The method involves:
#' \itemize{
#'   \item Using the convolution formula: \deqn{P(T_1 - T_2 ≥ q) = \int_{-\infty}^{\infty} f_1(x) \cdot F_2(x - q) dx}
#'   \item Where \eqn{f_1(x)} is the PDF of the first t-distribution
#'   \item And \eqn{F_2(x - q)} is the CDF of the second t-distribution evaluated at \eqn{x - q}
#'   \item Adaptive integration bounds based on the distribution characteristics
#'   \item High-precision numerical integration with specified tolerances
#' }
#'
#' This approach provides exact results (within numerical precision) but is more
#' computationally intensive than approximation methods like Welch-Satterthwaite.
#'
#' @examples
#' # Calculate P(t1 - t2 ≥ 3) for equal parameters
#' pNIdifft(q = 3, mu.t1 = 2, mu.t2 = 0, sd.t1 = 1, sd.t2 = 1, nu.t1 = 17, nu.t2 = 17)
#'
#' # Calculate P(t1 - t2 ≥ 1) for unequal variances
#' pNIdifft(q = 1, mu.t1 = 5, mu.t2 = 3, sd.t1 = 2, sd.t2 = 1.5, nu.t1 = 10, nu.t2 = 15)
#'
#' @importFrom stats dt pt integrate
#' @export
pNIdifft <- function(q, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2) {
  # Determine the number of results to compute
  n <- max(length(mu.t1), length(mu.t2), length(sd.t1), length(sd.t2))

  # Ensure all parameters have the same length
  mu.t1 <- rep(mu.t1, length.out = n)
  mu.t2 <- rep(mu.t2, length.out = n)
  sd.t1 <- rep(sd.t1, length.out = n)
  sd.t2 <- rep(sd.t2, length.out = n)

  # Calculate P(T1 - T2 ≥ q) for each parameter set
  results <- sapply(seq(n), function(i) {
    # Determine integration bounds based on the first distribution's characteristics
    # Calculate the center and spread of the first t-distribution
    center1 <- mu.t1[i]
    spread1 <- sd.t1[i] * sqrt((nu.t1 + 1) / (nu.t1 - 2))

    # Perform numerical integration using the convolution formula
    integrate(
      function(x) {
        # PDF of the first t-distribution at x
        f1 <- dt((x - mu.t1[i]) / sd.t1[i], df = nu.t1) / sd.t1[i]

        # CDF of the second t-distribution at (x - q)
        F2 <- pt(((x - q) - mu.t2[i]) / sd.t2[i], df = nu.t2)

        # Return the product for integration
        return(f1 * F2)
      },
      # Set adaptive integration bounds to capture the essential mass
      lower = max(-Inf, center1 - 8 * spread1),
      upper = min(Inf,  center1 + 8 * spread1),
      # Set high precision tolerances for accurate results
      rel.tol = 1e-6,
      abs.tol = 1e-8
    )[['value']]
  })

  return(results)
}
