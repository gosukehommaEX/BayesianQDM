#' Beta prior for a response rate
#'
#' Creates a Beta prior distribution for the response rate of one group in an
#' analysis of a single binary endpoint.
#'
#' @param a,b Positive shape parameters. \code{prior_beta(0.5, 0.5)} is the
#'   Jeffreys prior and \code{prior_beta(1, 1)} is the uniform prior.
#'
#' @return An object of class \code{c("qdm_prior_beta", "qdm_prior")} with
#'   elements \code{a} and \code{b}.
#'
#' @family priors
#' @examples
#' prior_beta(0.5, 0.5)
#' @export
prior_beta <- function(a, b) {
  check_numeric(a, "a", lower = 0, lower_open = TRUE)
  check_numeric(b, "b", lower = 0, lower_open = TRUE)
  structure(list(a = a, b = b), class = c("qdm_prior_beta", "qdm_prior"))
}
