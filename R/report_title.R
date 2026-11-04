#' Default title of a report
#'
#' @param x The object of the report (see \code{show_report()}).
#'
#' @return A character string.
#' @noRd
report_title <- function(x) {
  if (inherits(x, "qdm_oc")) {
    "Decision criteria and operating characteristics"
  } else if (inherits(x, "qdm_calibration")) {
    "Decision criteria and calibration of the thresholds"
  } else if (inherits(x, "qdm_boundary")) {
    "Decision criteria and decision boundaries"
  } else if (inherits(x, "qdm_decision")) {
    "Decision criteria and decision for the observed data"
  } else {
    "Decision criteria"
  }
}
