#' Sufficient statistics of simulated studies with two continuous endpoints
#'
#' With \eqn{\Sigma = L L^T} (Cholesky factor of the covariance matrix of the
#' scenario), the sample mean is \eqn{\bar y = \mu + L z / \sqrt{n}} with
#' \eqn{z = (z_1, z_2)}, and the matrix of sums of squares and
#' cross-products is \eqn{S = L A A^T L^T} with the Bartlett factor
#' \eqn{A = ((\sqrt{x_1}, 0), (z_{12}, \sqrt{x_2}))}, so that \eqn{S} has the
#' Wishart distribution with \eqn{n - 1} degrees of freedom and scale matrix
#' \eqn{\Sigma}, independently of \eqn{\bar y}.
#'
#' @param draws Random numbers from \code{draws_continuous2()}.
#' @param n Number of subjects.
#' @param mu,sd Means and standard deviations of the two endpoints.
#' @param rho Correlation of the two endpoints.
#'
#' @return A list with elements \code{mean} (matrix with two columns) and
#'   \code{cov} (sample covariance matrices, matrix with the columns
#'   \code{c11}, \code{c12} and \code{c22}).
#' @noRd
stats_continuous2 <- function(draws, n, mu, sd, rho) {
  l11 <- sd[1L]
  l21 <- rho * sd[2L]
  l22 <- sd[2L] * sqrt(1 - rho ^ 2)
  mean <- cbind(mu[1L] + l11 * draws$z1 / sqrt(n),
                mu[2L] + (l21 * draws$z1 + l22 * draws$z2) / sqrt(n))
  w11 <- draws$x1
  w12 <- sqrt(draws$x1) * draws$z12
  w22 <- draws$z12 ^ 2 + draws$x2
  cov <- cbind(c11 = l11 ^ 2 * w11,
               c12 = l11 * (l21 * w11 + l22 * w12),
               c22 = l21 ^ 2 * w11 + 2 * l21 * l22 * w12 + l22 ^ 2 * w22) /
    (n - 1)
  list(mean = mean, cov = cov)
}
