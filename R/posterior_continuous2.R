#' Posterior parameters of two continuous endpoints
#'
#' Updates the Normal-Inverse-Wishart parameters
#' \eqn{(\mu, \kappa, \nu, \Lambda)} with the external data (weight
#' \eqn{a_0}) and then with the current data (Yamaguchi et al., submitted,
#' equation A.10 and Section A.2.5). The vague prior
#' \eqn{|\Sigma|^{-3/2}} corresponds to \eqn{\kappa = 0}, \eqn{\nu = -1}
#' and \eqn{\Lambda = 0}. The current data can be given for many simulated
#' studies at once.
#'
#' @param prior An object of class \code{"qdm_prior_niw"} or
#'   \code{"qdm_prior_vague"}.
#' @param n Number of subjects of the current data (0 for none).
#' @param mean Sample means: a vector of length 2 or a matrix with two
#'   columns (one row per study).
#' @param cov Sample covariance matrices: a 2 x 2 matrix, or a matrix with
#'   the columns \code{c11}, \code{c12} and \code{c22} (one row per study).
#' @param external An object of class \code{"qdm_borrow"} with data of two
#'   continuous endpoints, or \code{NULL}.
#'
#' @return A list with elements \code{mu} (matrix with two columns),
#'   \code{kappa}, \code{nu} (numbers) and \code{Lambda} (matrix with the
#'   columns \code{L11}, \code{L12} and \code{L22}), one row per study.
#' @noRd
posterior_continuous2 <- function(prior, n, mean = c(0, 0),
                                  cov = matrix(0, 2L, 2L), external = NULL) {
  sym3 <- function(m) {
    if (is.matrix(m) && nrow(m) == 2L && ncol(m) == 2L) {
      matrix(c(m[1L, 1L], m[1L, 2L], m[2L, 2L]), 1L)
    } else {
      matrix(m, ncol = 3L)
    }
  }
  mean <- matrix(mean, ncol = 2L)
  cov <- sym3(cov)
  if (inherits(prior, "qdm_prior_vague")) {
    mu <- matrix(0, 1L, 2L)
    kappa <- 0
    nu <- -1
    lam <- matrix(0, 1L, 3L)
  } else {
    mu <- matrix(prior$mu0, 1L, 2L)
    kappa <- prior$kappa0
    nu <- prior$nu0
    lam <- sym3(prior$Lambda0)
  }
  updates <- list()
  if (!is.null(external) && external$a0 > 0) {
    updates <- list(list(w = external$a0, m = external$data$n,
                         ybar = matrix(external$data$mean, 1L),
                         s = sym3(external$data$cov)))
  }
  if (n > 0) {
    updates <- c(updates, list(list(w = 1, m = n, ybar = mean, s = cov)))
  }
  for (u in updates) {
    kappa_new <- kappa + u$w * u$m
    nr <- max(nrow(mu), nrow(u$ybar))
    d <- u$ybar[rep_len(seq_len(nrow(u$ybar)), nr), , drop = FALSE] -
      mu[rep_len(seq_len(nrow(mu)), nr), , drop = FALSE]
    f <- u$w * u$m * kappa / kappa_new
    lam <- lam[rep_len(seq_len(nrow(lam)), nr), , drop = FALSE] +
      u$w * (u$m - 1) * u$s[rep_len(seq_len(nrow(u$s)), nr), , drop = FALSE] +
      f * cbind(d[, 1L] ^ 2, d[, 1L] * d[, 2L], d[, 2L] ^ 2)
    mu <- mu[rep_len(seq_len(nrow(mu)), nr), , drop = FALSE] +
      u$w * u$m / kappa_new * d
    nu <- nu + u$w * u$m
    kappa <- kappa_new
  }
  dimnames(mu) <- NULL
  dimnames(lam) <- list(NULL, c("L11", "L12", "L22"))
  list(mu = mu, kappa = kappa, nu = nu, Lambda = lam)
}
