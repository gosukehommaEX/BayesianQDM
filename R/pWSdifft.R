#' Cumulative Distribution Function of the Welch-Satterthwaite Approximated Difference of Two t-Distributed Variables
#'
#' This function calculates the cumulative distribution function (CDF) of the difference
#' between two independent t-distributed random variables using the Welch-Satterthwaite
#' approximation. Specifically, it computes P(T1 - T2 > q) where T1 and T2 follow
#' t-distributions with potentially different degrees of freedom and scale parameters.
#'
#' @param q A numeric value representing the quantile threshold.
#' @param mu.t1 A numeric value representing the location parameter of the first t-distribution.
#' @param mu.t2 A numeric value representing the location parameter of the second t-distribution.
#' @param sd.t1 A positive numeric value representing the scale parameter of the first t-distribution.
#' @param sd.t2 A positive numeric value representing the scale parameter of the second t-distribution.
#' @param nu.t1 A positive numeric value representing the degrees of freedom of the first t-distribution.
#' @param nu.t2 A positive numeric value representing the degrees of freedom of the second t-distribution.
#'
#' @return A numeric value representing P(T1 - T2 > q), the probability that the difference
#'         between the two t-distributed variables exceeds the quantile q.
#'
#' @details
#' The Welch-Satterthwaite approximation is used to approximate the distribution of the
#' difference between two t-distributed variables. The method involves:
#' \itemize{
#'   \item Calculating the effective degrees of freedom using the Welch-Satterthwaite formula:
#'         \deqn{\nu^* = \frac{(s_1^2 + s_2^2)^2}{\frac{s_1^4}{\nu_1} + \frac{s_2^4}{\nu_2}}}
#'   \item Using a single t-distribution with adjusted parameters to approximate the difference
#'   \item Computing the standardized quantile and applying the t-distribution CDF
#' }
#'
#' This approximation is particularly useful for unequal variances and degrees of freedom,
#' providing a computationally efficient alternative to exact convolution methods.
#'
#' @examples
#' # Calculate P(t1 - t2 > 3) for equal parameters
#' pWSdifft(q = 3, mu.t1 = 2, mu.t2 = 0, sd.t1 = 1, sd.t2 = 1, nu.t1 = 17, nu.t2 = 17)
#'
#' # Calculate P(t1 - t2 > 1) for unequal variances
#' pWSdifft(q = 1, mu.t1 = 5, mu.t2 = 3, sd.t1 = 2, sd.t2 = 1.5, nu.t1 = 10, nu.t2 = 15)
#'
#' # Calculate P(t1 - t2 > 0) for different degrees of freedom
#' pWSdifft(q = 0, mu.t1 = 1, mu.t2 = 1, sd.t1 = 1, sd.t2 = 1, nu.t1 = 5, nu.t2 = 20)
#'
#' @importFrom stats pt
#' @export
pWSdifft <- function(q, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2) {
  # Calculate the adjusted scale parameter for the difference
  # This combines the variances of both t-distributions
  adj.scale <- sqrt(sd.t1 ^ 2 + sd.t2 ^ 2)

  # Calculate effective degrees of freedom using Welch-Satterthwaite formula
  # This accounts for the combination of two t-distributions with different parameters
  nu.star <- (sd.t1 ^ 2 + sd.t2 ^ 2) ^ 2 / ((sd.t1 ^ 4 / nu.t1) + (sd.t2 ^ 4 / nu.t2))

  # Calculate the location parameter for the difference T1 - T2
  location.param <- mu.t1 - mu.t2

  # Standardize the quantile using the adjusted scale and location
  standardized.q <- (q - location.param) / adj.scale

  # Return the upper tail probability P(T1 - T2 > q) using the approximated t-distribution
  results <- pt(standardized.q, df = nu.star, lower.tail = FALSE)

  return(results)
}
