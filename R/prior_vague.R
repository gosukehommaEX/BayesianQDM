#' Vague prior for continuous endpoints
#'
#' Creates the vague (reference) prior for the parameters of a normally
#' distributed endpoint. For one endpoint the prior density is
#' \eqn{p(\mu, \sigma^2) \propto \sigma^{-2}}; for two endpoints it is
#' \eqn{p(\mu, \Sigma) \propto |\Sigma|^{-3/2}}. With \eqn{n} observations the
#' posterior distribution of the mean is a t distribution with \eqn{n - 1}
#' degrees of freedom for one endpoint and a bivariate t distribution with
#' \eqn{n - 2} degrees of freedom for two endpoints.
#'
#' @return An object of class \code{c("qdm_prior_vague", "qdm_prior")}.
#'
#' @family priors
#' @examples
#' prior_vague()
#' @export
prior_vague <- function() {
  structure(list(), class = c("qdm_prior_vague", "qdm_prior"))
}
