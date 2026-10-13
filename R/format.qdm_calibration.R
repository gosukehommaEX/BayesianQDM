#' Format a calibration
#'
#' @param x An object of class \code{"qdm_calibration"}.
#' @param digits Number of significant digits for numeric values.
#' @param ... Not used.
#'
#' @return A character vector with one element per line of the description
#'   (without the table of decision probabilities).
#'
#' @examples
#' des <- design_binary(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                   gamma_go = 0.5, gamma_nogo = 0.8)
#' )
#' format(qdm_calibrate(des, data.frame(pi_t = 0.4, pi_c = 0.4),
#'                      data.frame(pi_t = 0.7, pi_c = 0.4),
#'                      target_go = 0.05, target_nogo = 0.2))
#' @export
format.qdm_calibration <- function(x, digits = 4L, ...) {
  describe <- function(name, gamma, prob, event, scenario, target) {
    if (is.na(gamma)) {
      sprintf("  %s: no value in the grid meets the target %s", name,
              format_num(target, digits))
    } else {
      sprintf("  %s = %s: P(%s criterion met | %s) = %s (target %s)", name,
              format_num(gamma, digits), event, scenario,
              format_num(prob, digits), format_num(target, digits))
    }
  }
  c(sprintf("Calibration of the decision thresholds (%s)", x$method),
    describe("gamma_go  ", x$gamma_go, x$prob_go, "Go", "scenario_go",
             x$target_go),
    describe("gamma_nogo", x$gamma_nogo, x$prob_nogo, "NoGo",
             "scenario_nogo", x$target_nogo))
}
