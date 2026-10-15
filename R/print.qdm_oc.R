#' Print operating characteristics
#'
#' @param x An object of class \code{"qdm_oc"}.
#' @param digits Number of significant digits for the probabilities.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' des <- design_binary(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                   gamma_go = 0.5, gamma_nogo = 0.8)
#' )
#' print(qdm_oc(des, data.frame(pi_t = c(0.4, 0.6), pi_c = 0.4)))
#' @export
print.qdm_oc <- function(x, digits = 4L, ...) {
  method <- attr(x, "method")
  if (!is.null(method)) {
    cat(sprintf("Operating characteristics (%s)\n", method))
  }
  y <- x
  class(y) <- "data.frame"
  attr(y, "design") <- NULL
  attr(y, "method") <- NULL
  print(y, digits = digits, row.names = FALSE)
  invisible(x)
}
