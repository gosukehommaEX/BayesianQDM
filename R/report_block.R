#' Block of a report
#'
#' A report is a list of blocks, which \code{report_text()} and
#' \code{report_html()} write as text or HTML. The text of the blocks may
#' contain formulas in LaTeX between dollar signs.
#'
#' @param type Type of the block: \code{"heading"} (element \code{text}),
#'   \code{"par"} (a paragraph, element \code{text}), \code{"list"} (element
#'   \code{items}, a character vector), \code{"table"} (elements
#'   \code{header}, a character vector, and \code{rows}, a character matrix
#'   with one column per element of \code{header}) or \code{"pre"}
#'   (preformatted text, element \code{lines}).
#' @param ... The elements of the block.
#'
#' @return A list with the element \code{type} and the elements in
#'   \code{...}.
#' @noRd
report_block <- function(type, ...) {
  c(list(type = type), list(...))
}
