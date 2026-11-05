#' Call as R code on one line
#'
#' @param node A character string of R code or a call from
#'   \code{code_node()}.
#'
#' @return A character string.
#' @noRd
code_flat <- function(node) {
  if (is.character(node)) {
    return(node)
  }
  nms <- names(node$args)
  if (is.null(nms)) {
    nms <- rep("", length(node$args))
  }
  parts <- vapply(seq_along(node$args), function(k) {
    v <- code_flat(node$args[[k]])
    if (nzchar(nms[k])) paste(nms[k], "=", v) else v
  }, character(1))
  paste0(node$fun, "(", paste(parts, collapse = ", "), ")")
}
