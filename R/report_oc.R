#' Operating characteristics section of a report
#'
#' @param x A result of \code{qdm_oc()}.
#' @param digits Number of decimal places of the probabilities.
#'
#' @return A list of blocks.
#' @noRd
report_oc <- function(x, digits) {
  method <- attr(x, "method")
  intro <- paste("The probabilities of Go, NoGo, Gray and Miss decisions",
                 "under scenarios for the true parameters")
  if (!is.null(method)) {
    intro <- sprintf("%s (%s)", intro, method)
  }
  list(report_block("heading", text = "Operating characteristics"),
       report_block("par", text = paste(intro, "are as follows.")),
       report_scenarios(x, digits),
       report_block("par", text = report_legend(attr(x, "design"))))
}
