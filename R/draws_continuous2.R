#' Standard random numbers for the simulation of studies with two continuous
#' endpoints
#'
#' Draws the random numbers from which the sufficient statistics of
#' \code{nsim} studies are computed under any scenario (see
#' \code{stats_continuous2()}), so that all scenarios use the same random
#' numbers. They are drawn in the order \code{z1}, \code{z2}, \code{x1},
#' \code{x2} and \code{z12}.
#'
#' @param nsim Number of simulated studies.
#' @param n Number of subjects (at least 3).
#'
#' @return A list with elements \code{z1}, \code{z2}, \code{z12} (standard
#'   normal), \code{x1} (chi-squared with \code{n - 1} degrees of freedom)
#'   and \code{x2} (chi-squared with \code{n - 2} degrees of freedom), each
#'   of length \code{nsim}.
#' @noRd
draws_continuous2 <- function(nsim, n) {
  z1 <- stats::rnorm(nsim)
  z2 <- stats::rnorm(nsim)
  x1 <- stats::rchisq(nsim, n - 1)
  x2 <- stats::rchisq(nsim, n - 2)
  z12 <- stats::rnorm(nsim)
  list(z1 = z1, z2 = z2, x1 = x1, x2 = x2, z12 = z12)
}
