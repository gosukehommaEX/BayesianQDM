#' Cumulative Distribution Function of the Difference of Two t-Distributed Variables by Monte Carlo Simulation
#'
#' This function calculates the cumulative distribution function (CDF) of the difference
#' between two independent t-distributed random variables using Monte Carlo simulation.
#' Specifically, it computes P(T1 - T2 <= q) or P(T1 - T2 > q) where T1 and T2 follow t-distributions
#' with potentially different location, scale, and degrees of freedom parameters.
#'
#' @param nMC A positive integer representing the number of Monte Carlo iterations for simulation.
#' @param q A numeric value representing the quantile threshold.
#' @param mu.t1 A numeric value representing the location parameter of the first t-distribution.
#' @param mu.t2 A numeric value representing the location parameter of the second t-distribution.
#' @param sd.t1 A positive numeric value representing the scale parameter of the first t-distribution.
#' @param sd.t2 A positive numeric value representing the scale parameter of the second t-distribution.
#' @param nu.t1 A positive numeric value representing the degrees of freedom of the first t-distribution.
#' @param nu.t2 A positive numeric value representing the degrees of freedom of the second t-distribution.
#' @param lower.tail logical; if TRUE (default), probabilities are P(T1 - T2 <= q), otherwise, P(T1 - T2 > q).
#'
#' @return A numeric value representing the probability that the difference
#'         between the two t-distributed variables below or exceeds the quantile q.
#'
#' @details
#' This function uses Monte Carlo simulation to approximate the distribution of the
#' difference between two t-distributed variables. The method involves:
#' \itemize{
#'   \item Generating nMC random samples from each t-distribution
#'   \item Computing the difference for each pair of samples
#'   \item Calculating the proportion of differences that exceed the threshold q
#'   \item The accuracy increases with larger nMC values
#' }
#'
#' Monte Carlo simulation provides a flexible and intuitive approach that can handle
#' any combination of parameters, but the accuracy depends on the number of simulations.
#' Larger nMC values yield more precise estimates but require more computational time.
#'
#' @examples
#' # Calculate P(t1 - t2 > 3) for equal parameters
#' pMCdifft(nMC = 1e+5, q = 3, mu.t1 = 2, mu.t2 = 0, sd.t1 = 1, sd.t2 = 1, nu.t1 = 17, nu.t2 = 17, lower.tail = FALSE)
#'
#' # Calculate P(t1 - t2 > 1) for unequal variances
#' pMCdifft(nMC = 1e+5, q = 1, mu.t1 = 5, mu.t2 = 3, sd.t1 = 2, sd.t2 = 1.5, nu.t1 = 10, nu.t2 = 15, lower.tail = FALSE)
#'
#' # Calculate P(t1 - t2 > 0) for different degrees of freedom
#' pMCdifft(nMC = 1e+5, q = 0, mu.t1 = 1, mu.t2 = 1, sd.t1 = 1, sd.t2 = 1, nu.t1 = 5, nu.t2 = 20, lower.tail = FALSE)
#'
#' @importFrom stats rt
#' @export
pMCdifft <- function(nMC, q, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2, lower.tail = TRUE) {
  # Determine the number of results to compute
  n = max(length(mu.t1), length(mu.t2), length(sd.t1), length(sd.t2))

  # Ensure all parameters have the same length
  mu.t1 <- rep(mu.t1, length.out = n)
  mu.t2 <- rep(mu.t2, length.out = n)
  sd.t1 <- rep(sd.t1, length.out = n)
  sd.t2 <- rep(sd.t2, length.out = n)

  # Calculate P(T1 - T2 > q) for each parameter set using Monte Carlo simulation
  results <- sapply(seq(n), function(i) {
    # Generate nMC random samples from the first t-distribution
    # Scale and shift the standard t-distribution to match desired parameters
    rand.t1 <- rt(nMC, df = nu.t1) * sqrt(sd.t1[i] ^ 2) + mu.t1[i]

    # Generate nMC random samples from the second t-distribution
    # Scale and shift the standard t-distribution to match desired parameters
    rand.t2 <- rt(nMC, df = nu.t2) * sqrt(sd.t2[i] ^ 2) + mu.t2[i]

    # Calculate the empirical probability P(T1 - T2 > q)
    # Count how many differences exceed q and divide by total simulations
    sum(rand.t1 - rand.t2 > q) / nMC
  })

  # Return the final results
  results <- lower.tail + c(1, -1)[lower.tail + 1] * results

  return(results)
}
