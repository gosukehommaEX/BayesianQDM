#' Normal-Inverse-Wishart prior for two continuous endpoints
#'
#' Creates the conjugate Normal-Inverse-Wishart prior for the mean vector
#' \eqn{\mu} and covariance matrix \eqn{\Sigma} of a bivariate normal endpoint
#' in one group:
#' \deqn{\mu \mid \Sigma \sim N_2(\mu_0, \Sigma / \kappa_0), \quad
#'   \Sigma \sim \textrm{Inv-Wishart}(\nu_0, \Lambda_0).}
#'
#' @param mu0 Numeric vector of length 2, the prior mean of \eqn{\mu}.
#' @param kappa0 Positive prior sample size for \eqn{\mu}.
#' @param nu0 Prior degrees of freedom of the inverse Wishart distribution;
#'   must be greater than 1.
#' @param Lambda0 A 2 x 2 symmetric positive definite scale matrix.
#'
#' @return An object of class \code{c("qdm_prior_niw", "qdm_prior")} with
#'   elements \code{mu0}, \code{kappa0}, \code{nu0} and \code{Lambda0}.
#'
#' @family priors
#' @examples
#' prior_niw(mu0 = c(0, 0), kappa0 = 1, nu0 = 4, Lambda0 = diag(2))
#' @export
prior_niw <- function(mu0, kappa0, nu0, Lambda0) {
  check_numeric(mu0, "mu0", len = 2L)
  check_numeric(kappa0, "kappa0", lower = 0, lower_open = TRUE)
  check_numeric(nu0, "nu0", lower = 1, lower_open = TRUE)
  check_pd_matrix(Lambda0, "Lambda0")
  structure(list(mu0 = as.numeric(mu0), kappa0 = kappa0, nu0 = nu0,
                 Lambda0 = unname(Lambda0)),
            class = c("qdm_prior_niw", "qdm_prior"))
}
