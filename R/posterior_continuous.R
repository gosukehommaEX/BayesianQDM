#' Posterior parameters for a continuous endpoint
#'
#' Combines a Normal-Inverse-chi-squared prior (or the vague prior), external
#' data borrowed through a power prior and the current data into the
#' parameters of the Normal-Inverse-chi-squared posterior distribution of the
#' mean and variance of one group. Data with \eqn{m} subjects, sample mean
#' \eqn{\bar y} and sample standard deviation \eqn{s}, entered with weight
#' \eqn{w} (\eqn{w = a_0} for external data and 1 for current data), update
#' the parameters as
#' \deqn{\kappa' = \kappa + w m, \quad \nu' = \nu + w m, \quad
#'   \mu' = (\kappa \mu + w m \bar y) / \kappa',}
#' \deqn{\nu' \sigma'^2 = \nu \sigma^2 + w (m - 1) s^2 +
#'   w m \kappa (\bar y - \mu)^2 / \kappa'.}
#' The vague prior \eqn{p(\mu, \sigma^2) \propto \sigma^{-2}} corresponds to
#' \eqn{\kappa = 0}, \eqn{\nu = -1} and \eqn{\nu \sigma^2 = 0}, which gives the
#' posterior distribution of Huang et al. (submitted, Theorem 4) with external
#' data and \eqn{\mu \mid D \sim t_{n - 1}(\bar y, s^2 / n)} without. The
#' marginal posterior distribution of the mean is
#' \eqn{t_{\nu}(\mu, \sigma^2 / \kappa)}.
#'
#' @param prior An object of class \code{"qdm_prior_nix"} or
#'   \code{"qdm_prior_vague"}.
#' @param n Number of subjects in the current data (0 if there are none).
#' @param mean,sd Numeric vectors of sample means and sample standard
#'   deviations of the current data (not used when \code{n} is 0).
#' @param external An object of class \code{"qdm_borrow"} or \code{NULL}.
#'
#' @return A list with elements \code{mu} and \code{sigma2} (numeric vectors
#'   of the length of \code{mean}, or of length 1 when \code{n} is 0),
#'   \code{kappa} and \code{nu} (numbers).
#' @noRd
posterior_continuous <- function(prior, n, mean = 0, sd = 0,
                                 external = NULL) {
  if (inherits(prior, "qdm_prior_vague")) {
    mu <- 0
    kappa <- 0
    nu <- -1
    ss <- 0
  } else {
    mu <- prior$mu0
    kappa <- prior$kappa0
    nu <- prior$nu0
    ss <- prior$nu0 * prior$sigma0 ^ 2
  }
  updates <- list()
  if (!is.null(external) && external$a0 > 0) {
    updates <- list(list(w = external$a0, m = external$data$n,
                         ybar = external$data$mean, s = external$data$sd))
  }
  if (n > 0) {
    updates <- c(updates, list(list(w = 1, m = n, ybar = mean, s = sd)))
  }
  for (u in updates) {
    kappa_new <- kappa + u$w * u$m
    ss <- ss + u$w * (u$m - 1) * u$s ^ 2 +
      u$w * u$m * kappa * (u$ybar - mu) ^ 2 / kappa_new
    mu <- (kappa * mu + u$w * u$m * u$ybar) / kappa_new
    nu <- nu + u$w * u$m
    kappa <- kappa_new
  }
  list(mu = mu, kappa = kappa, nu = nu, sigma2 = ss / nu)
}
