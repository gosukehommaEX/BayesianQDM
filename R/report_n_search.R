#' Section of a report on the search for the sample size
#'
#' @param x An object of class \code{"qdm_n_search"} with a chosen sample
#'   size.
#' @param digits Number of decimal places of the probabilities.
#'
#' @return A list of blocks: how the sample size was chosen, the
#'   probabilities for the candidates and the decision probabilities under
#'   the two scenarios.
#' @noRd
report_n_search <- function(x, digits) {
  design <- x$design
  num <- report_num
  prob <- function(p) formatC(p, format = "f", digits = digits)
  cols <- scenario_columns(design)$names
  scenario <- function(s) {
    paste(sprintf("$%s = %s$", report_symbol(cols),
                  vapply(cols, function(nm) num(s[[nm]]), character(1))),
          collapse = ", ")
  }
  curve <- x$curve
  chosen <- sprintf(paste("The sample sizes $n_t = %d$ and $n_c = %d$ were",
                          "chosen from %d candidates, with $n_t$ from %d to",
                          "%d (%s)."),
                    x$n_t, x$n_c, nrow(curve), min(curve$n_t),
                    max(curve$n_t), x$method)
  rule <- sprintf(paste(
    "They are the smallest candidate for which the probability that the Go",
    "criterion is met is at most %s under the scenario %s (the probability",
    "is %s) and the probability that the NoGo criterion is met is at most %s",
    "under the scenario %s (the probability is %s)."
  ), num(x$target_go), scenario(x$scenario_go), prob(x$prob_go),
  num(x$target_nogo), scenario(x$scenario_nogo), prob(x$prob_nogo))
  if (is.na(x$n_t_stable)) {
    rule <- paste(rule, "The largest candidate does not meet the targets.")
  } else if (x$n_t_stable != x$n_t) {
    rule <- paste(rule, sprintf(paste(
      "The probabilities are not monotone in the sample size: all candidates",
      "from $n_t = %d$ and $n_c = %d$ on meet the targets."
    ), x$n_t_stable, x$n_c_stable))
  }
  rows <- cbind(as.character(curve$n_t), as.character(curve$n_c),
                prob(curve$prob_go), prob(curve$prob_nogo),
                ifelse(curve$met, "yes", "no"))
  oc <- x$oc
  scen_labels <- false_decision_labels()$scenario
  oc$scenario <- ifelse(oc$scenario == "scenario_go", scen_labels[1L],
                        scen_labels[2L])
  names(oc)[names(oc) == "scenario"] <- "Scenario"
  list(report_block("heading", text = "Sample size"),
       report_block("par", text = paste(chosen, rule)),
       report_block("table",
                    header = c("$n_t$", "$n_c$", "P(Go criterion met)",
                               "P(NoGo criterion met)", "Targets met"),
                    rows = unname(rows)),
       report_block("par", text = paste("The probabilities of the decisions",
                                        "with the chosen sample sizes are as",
                                        "follows.")),
       report_scenarios(oc, digits),
       report_block("par", text = report_legend(design)))
}
