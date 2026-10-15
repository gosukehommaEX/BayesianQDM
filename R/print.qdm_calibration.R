#' Print a calibration
#'
#' Prints the chosen thresholds and the decision probabilities of the design
#' with these thresholds under the two scenarios.
#'
#' @param x An object of class \code{"qdm_calibration"}.
#' @param digits Number of significant digits for numeric values.
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
#' print(qdm_calibrate(des, data.frame(pi_t = 0.4, pi_c = 0.4),
#'                     data.frame(pi_t = 0.7, pi_c = 0.4),
#'                     target_go = 0.05, target_nogo = 0.2))
#' @export
print.qdm_calibration <- function(x, digits = 4L, ...) {
  cat(format(x, digits = digits), sep = "\n")
  if (!is.null(x$oc)) {
    cat("Decision probabilities with the chosen thresholds\n")
    print(x$oc, digits = digits, row.names = FALSE)
  }
  invisible(x)
}
