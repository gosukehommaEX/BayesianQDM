#' Cumulative Distribution Function of the Difference of Two t-Distributed Variables
#' by Monte Carlo Simulation
#'
#' This function calculates the cumulative distribution function (CDF) of the difference
#' between two independent t-distributed random variables using Monte Carlo simulation.
#' This method provides a flexible approach that can handle any combination of parameters,
#' with accuracy improving as the number of simulations increases. Specifically, it computes
#' P(T1 - T2 <= q) or P(T1 - T2 > q) where T1 and T2 follow t-distributions with potentially
#' different location, scale, and degrees of freedom parameters.
#'
#' @param nMC A positive integer representing the number of Monte Carlo iterations for
#'        simulation. Typical values range from 10,000 (for quick estimates) to 100,000+
#'        (for high precision). Larger values yield more accurate results but require
#'        more computational time.
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
#' @return A numeric value or vector in \code{[0, 1]} representing the estimated cumulative
#'         probability that the difference between the two t-distributed variables is below
#'         (if lower.tail = TRUE) or exceeds (if lower.tail = FALSE) the quantile q.
#'         If mu.t1, mu.t2, sd.t1, or sd.t2 are vectors of length n, returns a vector of
#'         length n. The estimate is subject to Monte Carlo error that decreases as
#'         sqrt(1/nMC).
#'
#' @details
#' This function uses Monte Carlo simulation to approximate the distribution of the
#' difference between two t-distributed variables. The algorithm consists of:
#' \itemize{
#'   \item **Step 1**: Generate an nMC x n matrix of random samples from the first
#'         t-distribution T1 ~ t(nu_1), then scale each column by sd.t1[i] and shift
#'         by mu.t1[i] simultaneously
#'   \item **Step 2**: Generate an nMC x n matrix of random samples from the second
#'         t-distribution T2 ~ t(nu_2), then scale and shift similarly
#'   \item **Step 3**: Compute the difference matrix D = T1 - T2
#'   \item **Step 4**: Calculate the empirical probability as the column-wise proportion
#'         of differences that satisfy the condition (D > q or D <= q)
#' }
#'
#' When mu.t1, mu.t2, sd.t1, or sd.t2 are vectors of length n, all n parameter sets are
#' processed simultaneously using matrix operations instead of a loop, which is
#' substantially faster for large n (e.g., n = nsim in pGNGsinglecontinuous).
#'
#' **Monte Carlo error**:
#' \itemize{
#'   \item The standard error of the estimate is approximately sqrt(p(1-p)/nMC),
#'         where p is the true probability
#'   \item For a probability near 0.5, the standard error is roughly 0.5/sqrt(nMC)
#'   \item For nMC = 10,000: SE approximately 0.005 (plus/minus 0.01 with 95 percent confidence)
#'   \item For nMC = 100,000: SE approximately 0.0016 (plus/minus 0.003 with 95 percent confidence)
#' }
#'
#' **Advantages**:
#' \itemize{
#'   \item Highly flexible - works with any parameter combination
#'   \item Easy to understand and implement
#'   \item Naturally handles complex scenarios
#'   \item Provides approximate confidence intervals via bootstrap
#' }
#'
#' **Computational considerations**:
#' \itemize{
#'   \item Computational time scales linearly with nMC
#'   \item Suitable for moderate precision requirements
#'   \item For highest precision, consider numerical integration (pNI2tdiff)
#'   \item Results are stochastic - different runs yield slightly different estimates
#' }
#'
#' @examples
#' # Calculate P(t1 - t2 > 3) for equal parameters
#' pMC2tdiff(nMC = 1e5, q = 3, mu.t1 = 2, mu.t2 = 0, sd.t1 = 1, sd.t2 = 1,
#'          nu.t1 = 17, nu.t2 = 17, lower.tail = FALSE)
#'
#' # Calculate P(t1 - t2 > 1) for unequal variances
#' pMC2tdiff(nMC = 1e5, q = 1, mu.t1 = 5, mu.t2 = 3, sd.t1 = 2, sd.t2 = 1.5,
#'          nu.t1 = 10, nu.t2 = 15, lower.tail = FALSE)
#'
#' # Calculate P(t1 - t2 > 0) for different degrees of freedom
#' pMC2tdiff(nMC = 1e5, q = 0, mu.t1 = 1, mu.t2 = 1, sd.t1 = 1, sd.t2 = 1,
#'          nu.t1 = 5, nu.t2 = 20, lower.tail = FALSE)
#'
#' # Calculate lower tail probability P(t1 - t2 <= 2)
#' pMC2tdiff(nMC = 1e5, q = 2, mu.t1 = 3, mu.t2 = 0, sd.t1 = 1.5, sd.t2 = 1.2,
#'          nu.t1 = 12, nu.t2 = 15, lower.tail = TRUE)
#'
#' # Vectorized usage: calculate for multiple parameter sets simultaneously
#' pMC2tdiff(nMC = 1e5, q = 1, mu.t1 = c(2, 3, 4), mu.t2 = c(0, 1, 2),
#'          sd.t1 = c(1, 1.2, 1.5), sd.t2 = c(1, 1.1, 1.3),
#'          nu.t1 = 10, nu.t2 = 10, lower.tail = FALSE)
#'
#' @importFrom stats rt
#' @export
pMC2tdiff <- function(nMC, q, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2, lower.tail = TRUE) {

  # Determine number of parameter sets
  n <- max(length(mu.t1), length(mu.t2), length(sd.t1), length(sd.t2))

  # Recycle all location/scale parameters to length n
  mu.t1 <- rep(mu.t1, length.out = n)
  mu.t2 <- rep(mu.t2, length.out = n)
  sd.t1 <- rep(sd.t1, length.out = n)
  sd.t2 <- rep(sd.t2, length.out = n)

  # Generate nMC x n matrices of standard t random samples for both groups.
  # Each column corresponds to one parameter set; columns are scaled and shifted
  # simultaneously using rep(..., each = nMC) to broadcast location/scale vectors.
  rand.t1 <- matrix(rt(nMC * n, df = nu.t1), nrow = nMC, ncol = n) *
    rep(sd.t1, each = nMC) + rep(mu.t1, each = nMC)

  rand.t2 <- matrix(rt(nMC * n, df = nu.t2), nrow = nMC, ncol = n) *
    rep(sd.t2, each = nMC) + rep(mu.t2, each = nMC)

  # Compute column-wise proportion of differences exceeding q:
  # P(T1 - T2 > q) for each parameter set
  results <- colMeans(rand.t1 - rand.t2 > q)

  # Convert to lower tail if requested:
  #   lower.tail = TRUE  -> return P(T1 - T2 <= q) = 1 - P(T1 - T2 > q)
  #   lower.tail = FALSE -> return P(T1 - T2 > q)
  results <- lower.tail + c(1, -1)[lower.tail + 1] * results

  return(results)
}
