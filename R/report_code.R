#' Section of a report with the R code of the design
#'
#' @param x The object of the report (see \code{show_report()}).
#'
#' @return A list of blocks: a heading, a sentence and the code of
#'   \code{show_code()}.
#' @noRd
report_code <- function(x) {
  if (inherits(x, "qdm_spec")) {
    x <- qdm_build(x)
  }
  list(report_block("heading", text = "R code"),
       report_block("par", text = paste("The following R code creates the",
                                        "design with BayesianQDM.")),
       report_block("pre", lines = unclass(show_code(x))))
}
