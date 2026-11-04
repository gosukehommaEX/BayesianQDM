#' Conduct of a time-to-event study in a report
#'
#' The description of \code{format.qdm_conduct()}, with each number formatted
#' separately (see \code{report_num()}), so that the numbers of a vector do
#' not share the number of decimals.
#'
#' @param conduct Conduct of a study, created by \code{conduct_tte()}.
#'
#' @return A character vector: the enrollment, the dropout and the time of
#'   the analysis.
#' @noRd
report_conduct <- function(conduct) {
  num <- report_num
  enroll <- if (length(conduct$accrual_rate) == 1L) {
    sprintf("Enrollment: uniform from 0 to %s", num(conduct$accrual_time[2L]))
  } else {
    sprintf("Enrollment: piecewise uniform on (%s) with rates (%s)",
            num(conduct$accrual_time), num(conduct$accrual_rate))
  }
  dropout <- if (all(conduct$dropout_hazard == 0)) {
    "Dropout: none"
  } else if (length(conduct$dropout_hazard) == 1L) {
    sprintf("Dropout: hazard %s", num(conduct$dropout_hazard))
  } else {
    sprintf("Dropout: hazards (%s) from times (%s)",
            num(conduct$dropout_hazard), num(conduct$dropout_time))
  }
  analysis <- if (is.null(conduct$cut_time)) {
    sprintf("Analysis: at event %d", conduct$cut_events)
  } else {
    sprintf("Analysis: at time %s", num(conduct$cut_time))
  }
  c(enroll, dropout, analysis)
}
