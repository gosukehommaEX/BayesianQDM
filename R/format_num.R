#' Format numbers as a comma-separated string
#'
#' @param x A numeric vector.
#' @param digits Number of significant digits.
#'
#' @return A single character string.
#' @noRd
format_num <- function(x, digits = 4L) {
  paste(format(x, digits = digits, trim = TRUE), collapse = ", ")
}
