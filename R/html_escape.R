#' Escape text for HTML
#'
#' Replaces \code{&}, \code{<} and \code{>} by their HTML entities. The
#' formulas in LaTeX are escaped as well: MathJax reads the text after the
#' entities are decoded.
#'
#' @param x Character vector.
#'
#' @return \code{x} with the characters replaced.
#' @noRd
html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}
