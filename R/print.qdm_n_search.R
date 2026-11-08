#' Print a search for the sample size
#'
#' @param x An object of class \code{"qdm_n_search"}.
#' @param digits Number of significant digits.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#' @examples
#' des <- design_binary(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                   gamma_go = 0.5, gamma_nogo = 0.8)
#' )
#' qdm_find_n(des, data.frame(pi_t = 0.4, pi_c = 0.4),
#'            data.frame(pi_t = 0.7, pi_c = 0.4), target_go = 0.1,
#'            target_nogo = 0.1, n_t = seq(6, 20, by = 2))
#' @export
print.qdm_n_search <- function(x, digits = 4L, ...) {
  cat(format(x, digits = digits), sep = "\n")
  cat("Probabilities for the candidate sample sizes\n")
  print(x$curve, digits = digits, row.names = FALSE)
  if (!is.null(x$oc)) {
    cat("Decision probabilities with the chosen sample size\n")
    print(x$oc, digits = digits, row.names = FALSE)
  }
  invisible(x)
}
