#' Labels of the probabilities of a false Go and a false NoGo
#'
#' The plots and reports of a calibration (\code{qdm_calibrate()}) and of a
#' search for the sample size (\code{qdm_find_n()}) show the probability that
#' the Go criterion is met under \code{scenario_go}, a scenario under which a
#' Go decision would be wrong, and the probability that the NoGo criterion is
#' met under \code{scenario_nogo}, a scenario under which a NoGo decision
#' would be wrong. The labels describe them as the probabilities of a false
#' Go and a false NoGo, without the argument names.
#'
#' @return A list with \code{legend}, the legend labels of the plots, and
#'   \code{scenario}, the labels of the two scenarios in the tables of the
#'   reports (Go first, then NoGo).
#' @noRd
false_decision_labels <- function() {
  list(legend = c("False Go (Go criterion met when Go would be wrong)",
                  "False NoGo (NoGo criterion met when NoGo would be wrong)"),
       scenario = c("Go would be wrong", "NoGo would be wrong"))
}
