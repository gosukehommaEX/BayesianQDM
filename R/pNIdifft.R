#' Cumulative distribution function (CDF) of the difference of two t-distributed outcomes by convolution of a difference of the t-distributions.
#'
#' This function calculates CDF(s) of difference between two variables following t-distribution
#' using a numerical integration.
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
#' pNIdifft(q = 3, mu.t1 = 2, mu.t2 = 0, sd.t1 = 1, sd.t2 = 1, nu.t1 = 17, nu.t2 = 17)
#'
#' @importFrom stats dt pt integrate
#' @export
pNIdifft <- function(q, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2) {
  # Set the number of results
  n <- max(length(mu.t1), length(mu.t2), length(sd.t1), length(sd.t2))
  mu.t1 <- rep(mu.t1, length.out = n)
  mu.t2 <- rep(mu.t2, length.out = n)
  sd.t1 <- rep(sd.t1, length.out = n)
  sd.t2 <- rep(sd.t2, length.out = n)
  # Calculate Pr(t1 - t2 >= q)
  results <- sapply(seq(n), function(i) {
    ## Set a range of integral
    # Based on the center of distribution for group 1
    center1 <- mu.t1[i]
    spread1 <- sd.t1[i] * sqrt((nu.t1 + 1) / (nu.t1 - 2))
    # Run numerical integration
    integrate(
      function(x) {
        f1 <- dt((x - mu.t1[i]) / sd.t1[i], df = nu.t1) / sd.t1[i]
        F2 <- pt(((x - q) - mu.t2[i]) / sd.t2[i], df = nu.t2)
        return(f1 * F2)
      },
      # Adjust the range of integral dynamically
      lower = max(-Inf, center1 - 8 * spread1),
      upper = min(Inf,  center1 + 8 * spread1),
      rel.tol = 1e-6,
      abs.tol = 1e-8
    )[['value']]
  })
  return(results)
}
