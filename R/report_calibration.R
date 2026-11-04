#' Calibration section of a report
#'
#' @param x A result of \code{qdm_calibrate()} in which both thresholds were
#'   found.
#' @param digits Number of decimal places of the probabilities.
#'
#' @return A list of blocks.
#' @noRd
report_calibration <- function(x, digits) {
  design <- x$design
  num <- report_num
  prob <- function(p) formatC(p, format = "f", digits = digits)
  cols <- scenario_columns(design)$names
  scenario <- function(s) {
    paste(sprintf("$%s = %s$", report_symbol(cols),
                  vapply(cols, function(nm) num(s[[nm]]), character(1))),
          collapse = ", ")
  }
  grid <- x$curve$gamma
  chosen <- sprintf(paste("The thresholds $\\gamma_{go} = %s$ and",
                          "$\\gamma_{nogo} = %s$ were chosen from %d",
                          "candidate values from %s to %s (%s)."),
                    num(x$gamma_go), num(x$gamma_nogo), length(grid),
                    num(min(grid)), num(max(grid)), x$method)
  rule <- sprintf(paste(
    "$\\gamma_{go}$ is the smallest candidate for which the probability that",
    "the Go criterion is met is at most %s under the scenario %s (the",
    "probability is %s), and $\\gamma_{nogo}$ is the smallest candidate for",
    "which the probability that the NoGo criterion is met is at most %s",
    "under the scenario %s (the probability is %s)."
  ), num(x$target_go), scenario(x$scenario_go), prob(x$prob_go),
  num(x$target_nogo), scenario(x$scenario_nogo), prob(x$prob_nogo))
  oc <- x$oc
  oc$scenario <- ifelse(oc$scenario == "scenario_go",
                        "Scenario of $\\gamma_{go}$",
                        "Scenario of $\\gamma_{nogo}$")
  names(oc)[names(oc) == "scenario"] <- "Scenario"
  list(report_block("heading", text = "Calibration of the thresholds"),
       report_block("par", text = paste(chosen, rule)),
       report_block("par", text = paste("The probabilities of the decisions",
                                        "with these thresholds are as",
                                        "follows.")),
       report_scenarios(oc, digits),
       report_block("par", text = report_legend(design)))
}
