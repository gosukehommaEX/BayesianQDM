#' Table of decision probabilities under scenarios in a report
#'
#' @param x A data frame with the columns of the scenarios (and possibly
#'   other columns) and the columns \code{Go}, \code{NoGo}, \code{Gray} and
#'   \code{Miss}, such as a result of \code{qdm_oc()}.
#' @param digits Number of decimal places of the probabilities.
#'
#' @return A block of type \code{"table"}. The columns of the scenarios have
#'   their symbols as headers (see \code{report_symbol()}), and missing
#'   probabilities (scenarios that cannot be evaluated) are shown as
#'   \code{"NA"}.
#' @noRd
report_scenarios <- function(x, digits) {
  probs <- c("Go", "NoGo", "Gray", "Miss")
  cols <- setdiff(names(x), probs)
  sym <- report_symbol(cols)
  header <- c(ifelse(is.na(sym), cols, paste0("$", sym, "$")), probs)
  n <- nrow(x)
  scen <- vapply(cols, function(nm) {
    v <- x[[nm]]
    if (is.numeric(v)) vapply(v, report_num, character(1)) else
      as.character(v)
  }, character(n))
  pr <- vapply(probs, function(nm) {
    p <- x[[nm]]
    ifelse(is.na(p), "NA", formatC(p, format = "f", digits = digits))
  }, character(n))
  rows <- cbind(matrix(scen, nrow = n), matrix(pr, nrow = n))
  report_block("table", header = unname(header), rows = unname(rows))
}
