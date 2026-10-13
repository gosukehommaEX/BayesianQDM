#' Normal-Inverse-chi-squared prior for a continuous endpoint
#'
#' Creates the conjugate Normal-Inverse-chi-squared prior for the mean
#' \eqn{\mu} and variance \eqn{\sigma^2} of a normally distributed endpoint in
#' one group:
#' \deqn{\mu \mid \sigma^2 \sim N(\mu_0, \sigma^2 / \kappa_0), \quad
#'   \sigma^2 \sim \textrm{Inv-}\chi^2(\nu_0, \sigma_0^2).}
#'
#' @param mu0 Prior mean of \eqn{\mu}.
#' @param kappa0 Positive prior sample size for \eqn{\mu}.
#' @param nu0 Positive prior degrees of freedom for \eqn{\sigma^2}.
#' @param sigma0 Positive prior scale of \eqn{\sigma}, on the standard
#'   deviation scale (the scaled inverse chi-squared distribution has scale
#'   parameter \code{sigma0^2}).
#'
#' @return An object of class \code{c("qdm_prior_nix", "qdm_prior")} with
#'   elements \code{mu0}, \code{kappa0}, \code{nu0} and \code{sigma0}.
#'
#' @family priors
#' @examples
#' prior_nix(mu0 = 0, kappa0 = 1, nu0 = 3, sigma0 = 2)
#' @export
prior_nix <- function(mu0, kappa0, nu0, sigma0) {
  check_numeric(mu0, "mu0")
  check_numeric(kappa0, "kappa0", lower = 0, lower_open = TRUE)
  check_numeric(nu0, "nu0", lower = 0, lower_open = TRUE)
  check_numeric(sigma0, "sigma0", lower = 0, lower_open = TRUE)
  structure(list(mu0 = mu0, kappa0 = kappa0, nu0 = nu0, sigma0 = sigma0),
            class = c("qdm_prior_nix", "qdm_prior"))
}
