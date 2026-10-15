#' Format the conduct of a study
#'
#' @param x An object of class \code{"qdm_conduct"}.
#' @param digits Number of significant digits for numeric values.
#' @param ... Not used.
#'
#' @return A character vector with one element per line of the description.
#'
#' @examples
#' format(conduct_tte(accrual_time = c(0, 12), cut_events = 20))
#' @export
format.qdm_conduct <- function(x, digits = 4L, ...) {
  num <- function(v) format_num(v, digits)
  enroll <- if (length(x$accrual_rate) == 1L) {
    sprintf("Enrollment: uniform from 0 to %s", num(x$accrual_time[2L]))
  } else {
    sprintf("Enrollment: piecewise uniform on (%s) with rates (%s)",
            num(x$accrual_time), num(x$accrual_rate))
  }
  dropout <- if (all(x$dropout_hazard == 0)) {
    "Dropout: none"
  } else if (length(x$dropout_hazard) == 1L) {
    sprintf("Dropout: hazard %s", num(x$dropout_hazard))
  } else {
    sprintf("Dropout: hazards (%s) from times (%s)", num(x$dropout_hazard),
            num(x$dropout_time))
  }
  analysis <- if (is.null(x$cut_time)) {
    sprintf("Analysis: at event %d", x$cut_events)
  } else {
    sprintf("Analysis: at time %s", num(x$cut_time))
  }
  c(enroll, dropout, analysis)
}
