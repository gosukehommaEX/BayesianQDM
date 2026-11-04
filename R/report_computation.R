#' Computation section of a report
#'
#' @param design A design.
#'
#' @return A list of blocks: the method of computation of the probabilities
#'   and, for a time-to-event endpoint, the conduct of the study and of the
#'   future study.
#' @noRd
report_computation <- function(design) {
  crit <- design$criterion
  fixed <- design$control == "fixed"
  posterior <- crit$type == "posterior"
  how <- switch(
    design$endpoint,
    binary = if (fixed) {
      "exact computation (fixed control)"
    } else if (posterior) {
      "numerical integration"
    } else {
      "exact summation"
    },
    continuous = method_label(design$method, design$nMC, fixed),
    continuous2 = method_label(design$method, design$nMC, fixed),
    binary2 = if (posterior) method_label("MC", design$nMC) else
      "exact computation",
    tte = tte_method_label(design)
  )
  out <- list(report_block("heading",
                           text = "Computation of the probabilities"),
              report_block("par", text = sprintf(
                "The probabilities of the decision criteria are computed by %s.",
                how
              )))
  if (design$endpoint == "tte" && !is.null(design$conduct)) {
    out <- c(out, list(
      report_block("par", text = paste("The operating characteristics are",
                                       "computed by simulating the conduct",
                                       "of the study:")),
      report_block("list", items = report_conduct(design$conduct))
    ))
  }
  if (design$endpoint == "tte" && !posterior) {
    out <- c(out, list(
      report_block("par", text = "The future studies are simulated with:"),
      report_block("list", items = report_conduct(design$conduct_future))
    ))
  }
  out
}
