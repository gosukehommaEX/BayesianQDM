#' Print a prior specification
#'
#' @param x An object of class \code{"qdm_prior"}.
#' @param ... Passed to \code{\link{format.qdm_prior}}.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' print(prior_nix(mu0 = 0, kappa0 = 1, nu0 = 3, sigma0 = 2))
#' @export
print.qdm_prior <- function(x, ...) {
  cat("Prior:", format(x, ...), "\n")
  invisible(x)
}
