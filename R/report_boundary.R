#' Decision boundary section of a report
#'
#' @param x A result of \code{qdm_boundary()}.
#'
#' @return A list of blocks, with the tables of \code{format()} as
#'   preformatted text (all values of the conditioning variable), or the
#'   printed data frame when the attributes of \code{x} were lost by
#'   subsetting.
#' @noRd
report_boundary <- function(x) {
  lines <- if (is.null(attr(x, "axes"))) {
    utils::capture.output(print(as.data.frame(x)))
  } else {
    format(x, n = nrow(x))
  }
  list(report_block("heading", text = "Decision boundaries"),
       report_block("par", text = paste("The observed data that lead to",
                                        "each decision are as follows.")),
       report_block("pre", lines = lines))
}
