#' Calibration as a data frame
#'
#' Converts a result of \code{\link{qdm_calibrate}} to a data frame with one
#' row: the chosen thresholds, the probabilities of meeting the criteria at
#' these thresholds and the targets. The probabilities for all candidate
#' thresholds are in the element \code{curve} of the result, which is a data
#' frame.
#'
#' @param x An object of class \code{"qdm_calibration"}.
#' @param row.names \code{NULL} (the default) or the name of the row.
#' @param optional Not used.
#' @param ... Not used.
#'
#' @return A data frame with one row and the columns \code{gamma_go},
#'   \code{prob_go}, \code{target_go}, \code{gamma_nogo}, \code{prob_nogo}
#'   and \code{target_nogo} (see \code{\link{qdm_calibrate}}).
#'
#' @examples
#' des <- design_binary(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                   gamma_go = 0.5, gamma_nogo = 0.8)
#' )
#' cal <- qdm_calibrate(des,
#'                      scenario_go = data.frame(pi_t = 0.4, pi_c = 0.4),
#'                      scenario_nogo = data.frame(pi_t = 0.7, pi_c = 0.4),
#'                      target_go = 0.05, target_nogo = 0.2)
#' as.data.frame(cal)
#' @export
as.data.frame.qdm_calibration <- function(x, row.names = NULL,
                                          optional = FALSE, ...) {
  out <- data.frame(gamma_go = x$gamma_go, prob_go = x$prob_go,
                    target_go = x$target_go, gamma_nogo = x$gamma_nogo,
                    prob_nogo = x$prob_nogo, target_nogo = x$target_nogo)
  if (!is.null(row.names)) {
    rownames(out) <- row.names
  }
  out
}
