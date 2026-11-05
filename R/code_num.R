#' Numbers as R code
#'
#' Writes each number with the fewest significant digits that give it back
#' exactly when the code is parsed.
#'
#' @param x A numeric vector.
#' @param keep_type Whether to write integers with the suffix \code{L}, so
#'   that the code gives an integer vector. Arguments that the functions of
#'   the package convert to integers are written without it.
#'
#' @return A character string: the number, or a call of \code{c()} for
#'   several numbers.
#' @noRd
code_num <- function(x, keep_type = TRUE) {
  one <- function(v) {
    if (is.integer(v)) {
      return(if (keep_type) paste0(v, "L") else as.character(v))
    }
    exact <- function(t) identical(as.numeric(t), v)
    for (d in 1:16) {
      # formatC() can start the number with a blank (" 1" for 1)
      fixed <- trimws(formatC(v, digits = d, format = "fg"))
      sci <- trimws(formatC(v, digits = d, format = "g"))
      # Fixed notation unless the scientific one is much shorter
      if (exact(fixed) && (!exact(sci) || nchar(fixed) <= nchar(sci) + 3L)) {
        return(fixed)
      }
      if (exact(sci)) {
        return(sci)
      }
    }
    # 17 significant digits give every number exactly
    trimws(formatC(v, digits = 17L, format = "g"))
  }
  s <- vapply(as.list(unname(x)), one, character(1))
  if (length(s) == 1L) s else paste0("c(", paste(s, collapse = ", "), ")")
}
