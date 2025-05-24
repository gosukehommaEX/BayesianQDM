#' Cumulative distribution function (CDF) of the Welch-Satterthwaite approximated difference of two t-distributed outcomes.
#'
#' This function calculates CDF(s) of difference between two variables following t-distribution
#' using Welch-Satterthwaite approximation.
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
#' pWSdifft(q = 3, mu.t1 = 2, mu.t2 = 0, sd.t1 = 1, sd.t2 = 1, nu.t1 = 17, nu.t2 = 17)
#'
#' @importFrom stats pt
#' @export
pWSdifft = function(q, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2) {
  # Calculate adjusted scale parameter
  adj.scale = sqrt(sd.t1 ^ 2 + sd.t2 ^ 2)
  # Calculate effective degrees of freedom using Welch-Satterthwaite formula
  nu.star = (sd.t1 ^ 2 + sd.t2 ^ 2) ^ 2 / ((sd.t1 ^ 4 / nu.t1) + (sd.t2 ^ 4 / nu.t2))
  # Location parameter for the difference
  location.param = mu.t1 - mu.t2
  # Standardize the quantile
  standardized.q = (q - location.param) / adj.scale
  # Return CDF of t-distribution with effective degrees of freedom
  return(pt(standardized.q, df = nu.star, lower.tail = FALSE))
}
