#' Cumulative Distribution Function of the Difference of Two t-Distributed Variables
#' Using the Welch-Satterthwaite Approximation
#'
#' This function calculates the cumulative distribution function (CDF) of the difference
#' between two independent t-distributed random variables using the Welch-Satterthwaite
#' approximation. This method provides a fast and reasonably accurate approximation,
#' particularly useful for unequal variances and degrees of freedom. Specifically,
#' it computes P(T1 - T2 ≤ q) or P(T1 - T2 > q) where T1 and T2 follow t-distributions
#' with potentially different parameters.
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
#' @return A numeric value in \code{[0, 1]} representing the cumulative probability that the
#'         difference between the two t-distributed variables is below (if lower.tail = TRUE)
#'         or exceeds (if lower.tail = FALSE) the quantile q.
#'
#' @details
#' The Welch-Satterthwaite approximation approximates the distribution of the difference
#' between two t-distributed variables using a single t-distribution with adjusted
#' parameters. The method involves:
#' \itemize{
#'   \item **Step 1**: Calculate the pooled variance (sum of squared scale parameters):
#'         \deqn{\sigma_{\text{pooled}}^2 = \sigma_1^2 + \sigma_2^2}
#'   \item **Step 2**: Calculate the effective degrees of freedom using the
#'         Welch-Satterthwaite formula:
#'         \deqn{\nu^* = \frac{(\sigma_1^2 + \sigma_2^2)^2}{\frac{\sigma_1^4}{\nu_1} + \frac{\sigma_2^4}{\nu_2}}}
#'   \item **Step 3**: Standardize the quantile:
#'         \deqn{q_{\text{standardized}} = \frac{q - (\mu_1 - \mu_2)}{\sqrt{\sigma_1^2 + \sigma_2^2}}}
#'   \item **Step 4**: Compute the CDF using a single t-distribution with the
#'         effective degrees of freedom ν*
#' }
#'
#' **Advantages**:
#' \itemize{
#'   \item Computationally fast and efficient
#'   \item Works well for unequal variances and degrees of freedom
#'   \item Provides good approximation in most practical scenarios
#'   \item Particularly accurate when degrees of freedom are not too small (ν > 5)
#' }
#'
#' **Computational considerations**:
#' \itemize{
#'   \item Much faster than numerical integration methods
#'   \item Suitable for exploratory analyses and large-scale simulations
#'   \item For final analyses requiring highest accuracy, consider using
#'         numerical integration (pNIdifft)
#' }
#'
#' @examples
#' # Calculate P(t1 - t2 > 3) for equal parameters
#' pWS2tdiff(q = 3, mu.t1 = 2, mu.t2 = 0, sd.t1 = 1, sd.t2 = 1,
#'          nu.t1 = 17, nu.t2 = 17, lower.tail = FALSE)
#'
#' # Calculate P(t1 - t2 > 1) for unequal variances
#' pWS2tdiff(q = 1, mu.t1 = 5, mu.t2 = 3, sd.t1 = 2, sd.t2 = 1.5,
#'          nu.t1 = 10, nu.t2 = 15, lower.tail = FALSE)
#'
#' # Calculate P(t1 - t2 > 0) for different degrees of freedom - should be ~0.5
#' pWS2tdiff(q = 0, mu.t1 = 1, mu.t2 = 1, sd.t1 = 1, sd.t2 = 1,
#'          nu.t1 = 5, nu.t2 = 20, lower.tail = FALSE)
#'
#' # Calculate lower tail probability P(t1 - t2 ≤ 2)
#' pWS2tdiff(q = 2, mu.t1 = 3, mu.t2 = 0, sd.t1 = 1.5, sd.t2 = 1.2,
#'          nu.t1 = 12, nu.t2 = 15, lower.tail = TRUE)
#'
#' @importFrom stats pt
#' @export
pWS2tdiff <- function(q, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2, lower.tail = TRUE) {
  # Step 1: Calculate the adjusted scale parameter for the difference
  # This is the square root of the pooled variance (sum of individual variances)
  adj.scale <- sqrt(sd.t1 ^ 2 + sd.t2 ^ 2)

  # Step 2: Calculate effective degrees of freedom using Welch-Satterthwaite formula
  # This accounts for the combination of two t-distributions with different parameters
  # The formula weights each distribution's contribution by its variance
  nu.star <- (sd.t1 ^ 2 + sd.t2 ^ 2) ^ 2 / ((sd.t1 ^ 4 / nu.t1) + (sd.t2 ^ 4 / nu.t2))

  # Step 3: Calculate the location parameter for the difference T1 - T2
  # This is simply the difference of the two location parameters
  location.param <- mu.t1 - mu.t2

  # Step 4: Standardize the quantile using the adjusted scale and location
  # This transforms q to the standard form for the approximating t-distribution
  standardized.q <- (q - location.param) / adj.scale

  # Step 5: Calculate the probability using the approximating t-distribution
  # with the effective degrees of freedom
  results <- pt(standardized.q, df = nu.star, lower.tail = FALSE)

  # Convert to lower tail probability if requested
  # If lower.tail = TRUE: return 1 - P(T1 - T2 > q) = P(T1 - T2 ≤ q)
  # If lower.tail = FALSE: return P(T1 - T2 > q)
  results <- lower.tail + c(1, -1)[lower.tail + 1] * results

  return(results)
}
