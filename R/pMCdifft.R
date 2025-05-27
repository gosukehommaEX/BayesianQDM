#' Cumulative distribution function (CDF) of the difference of two t-distributed outcomes by Monte Carlo simulation.
#'
#' This function calculates CDF(s) of difference between two variables following t-distribution
#' using a Monte Carlo simulation (1,000,000 iterations).
#'
#' @param q Quantile value.
#' @param mu.t1 Location parameter of t-distribution for group 1.
#' @param mu.t2 Location parameter of t-distribution for group 2.
#' @param sd.t1 Scale parameter of t-distribution for group 1.
#' @param sd.t2 Scale parameter of t-distribution for group 2.
#' @param nu.t1 Degrees of freedom of t-distribution for group 1.
#' @param nu.t2 Degrees of freedom of t-distribution for group 2.
#'
#' @return Cumulative probability(s)
#'
#' @examples
#' pMCdifft(q = 3, mu.t1 = 2, mu.t2 = 0, sd.t1 = 1, sd.t2 = 1, nu.t1 = 17, nu.t2 = 17)
#'
#' @importFrom stats rt
#' @export
pMCdifft = function(q, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2) {
  # Set the number of results
  n = max(length(mu.t1), length(mu.t2), length(sd.t1), length(sd.t2))
  mu.t1 = rep(mu.t1, length.out = n)
  mu.t2 = rep(mu.t2, length.out = n)
  sd.t1 = rep(sd.t1, length.out = n)
  sd.t2 = rep(sd.t2, length.out = n)
  # Calculate Pr(t1 - t2 >= q)
  results = sapply(seq(n), function(i) {
    # t-distributed random numbers for group 1
    rand.t1 = rt(1e+6, df = nu.t1) * sqrt(sd.t1[i] ^ 2) + mu.t1[i]
    # t-distributed random numbers for group 1
    rand.t2 = rt(1e+6, df = nu.t2) * sqrt(sd.t2[i] ^ 2) + mu.t2[i]
    # CDF
    sum(rand.t1 - rand.t2 > q) / 1e+6
  })
  return(results)
}
