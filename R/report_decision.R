#' Decision section of a report
#'
#' @param x A result of \code{qdm_decide()}.
#' @param digits Number of decimal places of the probabilities.
#'
#' @return A list of blocks: the observed data, the probabilities of the
#'   criteria compared with the thresholds, the decision, the probabilities
#'   of the regions (two endpoints) and the posterior distributions.
#' @noRd
report_decision <- function(x, digits) {
  design <- x$design
  crit <- design$criterion
  prob <- function(p) formatC(p, format = "f", digits = digits)
  ev <- report_events(crit)
  join <- if (crit$n_endpoints == 2L && crit$rule == "or") " or " else
    " and "
  compare <- function(e, p, gamma) {
    gamma <- rep_len(gamma, length(e))
    op <- ifelse(meets_threshold(p, gamma), "\\ge", "<")
    paste(sprintf("$%s = %s %s %s$", e, prob(p), op,
                  vapply(gamma, report_num, character(1))), collapse = join)
  }
  met <- function(m) if (m) "met" else "not met"
  data <- if (is.null(x$data_c)) {
    sprintf("The observed data of the treatment group are %s.",
            report_data(x$data_t))
  } else {
    sprintf(paste("The observed data are %s in the treatment group and %s in",
                  "the control group."), report_data(x$data_t),
            report_data(x$data_c))
  }
  out <- list(
    report_block("heading", text = "Decision for the observed data"),
    report_block("par", text = data),
    report_block("list", items = c(
      sprintf("Go criterion (%s): %s", met(x$met_go),
              compare(ev$go, x$prob_go, crit$gamma_go)),
      sprintf("NoGo criterion (%s): %s", met(x$met_nogo),
              compare(ev$nogo, x$prob_nogo, crit$gamma_nogo))
    )),
    report_block("par", text = sprintf("The decision is %s.", x$decision))
  )
  if (!is.null(x$prob_regions)) {
    pr <- x$prob_regions
    out <- c(out, list(
      report_block("par", text = "The probabilities of the regions are:"),
      report_block("table", header = c("Region", "Probability"),
                   rows = cbind(sprintf("$R_{%s}$", sub("^R", "", names(pr))),
                                prob(unname(pr))))
    ))
  }
  post_t <- report_posterior(design, x$posterior_t, "t")
  # A fixed control value is given in the section on the statistical model
  posterior <- if (startsWith(names(x$posterior_c)[1L], "value")) {
    sprintf("The posterior distribution of %s is %s.", post_t[1L], post_t[2L])
  } else {
    post_c <- report_posterior(design, x$posterior_c, "c")
    sprintf("The posterior distribution of %s is %s, and that of %s is %s.",
            post_t[1L], post_t[2L], post_c[1L], post_c[2L])
  }
  c(out, list(report_block("par", text = posterior)))
}
