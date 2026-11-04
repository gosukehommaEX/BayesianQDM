#' Numbers of a report
#'
#' Formats each number with up to four significant digits, separately (so
#' that the numbers of a vector do not share the number of decimals), and
#' separates them by commas.
#'
#' @param x Numeric vector.
#'
#' @return A character string.
#' @noRd
report_num <- function(x) {
  paste(vapply(x, format_num, character(1), digits = 4L), collapse = ", ")
}
