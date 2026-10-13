#' Format a prior specification
#'
#' @param x An object of class \code{"qdm_prior"}.
#' @param digits Number of significant digits for numeric values.
#' @param ... Not used.
#'
#' @return A character string describing the prior distribution.
#'
#' @examples
#' format(prior_beta(0.5, 0.5))
#' @export
format.qdm_prior <- function(x, digits = 4L, ...) {
  switch(
    class(x)[1L],
    qdm_prior_beta = sprintf("Beta(a = %s, b = %s)",
                             format(x$a, digits = digits),
                             format(x$b, digits = digits)),
    qdm_prior_dirichlet = sprintf(
      "Dirichlet(alpha = (%s))",
      format_num(x$alpha, digits)
    ),
    qdm_prior_nix = sprintf(
      "N-Inv-Chisq(mu0 = %s, kappa0 = %s, nu0 = %s, sigma0 = %s)",
      format(x$mu0, digits = digits), format(x$kappa0, digits = digits),
      format(x$nu0, digits = digits), format(x$sigma0, digits = digits)
    ),
    qdm_prior_niw = sprintf(
      "N-Inv-Wishart(mu0 = (%s), kappa0 = %s, nu0 = %s, Lambda0 = (%s; %s))",
      format_num(x$mu0, digits),
      format(x$kappa0, digits = digits), format(x$nu0, digits = digits),
      format_num(x$Lambda0[1L, ], digits),
      format_num(x$Lambda0[2L, ], digits)
    ),
    qdm_prior_vague = "Vague",
    stop("Unknown prior class.", call. = FALSE)
  )
}
