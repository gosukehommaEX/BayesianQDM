#' Call as lines of R code
#'
#' Writes a call on one line if it fits in \code{width} characters from the
#' column \code{indent}. Otherwise the arguments are placed on lines aligned
#' after the opening parenthesis, as many on a line as fit, and an argument
#' that does not fit on a line is itself written over several lines. With
#' \code{hanging = TRUE}, the arguments start on the next line, indented by
#' two spaces, and the closing parenthesis is on a line of its own, as in the
#' examples of the package.
#'
#' @param node A character string of R code or a call from
#'   \code{code_node()}.
#' @param indent Column at which the call starts.
#' @param width Largest width of a line.
#' @param hanging Whether to indent the arguments by two spaces from the
#'   start of the line (for a call that starts a line).
#'
#' @return A character vector of lines. The first line has no indentation
#'   (the caller writes what comes before it); the other lines start with
#'   the spaces of their columns.
#' @noRd
code_format <- function(node, indent = 0L, width = 80L, hanging = FALSE) {
  flat <- code_flat(node)
  if (is.character(node) || indent + nchar(flat) <= width) {
    return(flat)
  }
  open <- paste0(node$fun, "(")
  args <- node$args
  n <- length(args)
  nms <- names(args)
  if (is.null(nms)) {
    nms <- rep("", n)
  }
  if (hanging) {
    inner <- 2L
    out <- open
    line <- strrep(" ", inner)
  } else {
    inner <- indent + nchar(open)
    out <- character(0)
    line <- open
  }
  col <- inner
  empty <- TRUE
  for (k in seq_len(n)) {
    lead <- if (nzchar(nms[k])) paste(nms[k], "= ") else ""
    tail <- if (k < n) "," else if (hanging) "" else ")"
    piece <- paste0(lead, code_flat(args[[k]]), tail)
    sep <- if (empty) "" else " "
    if (col + nchar(sep) + nchar(piece) <= width) {
      line <- paste0(line, sep, piece)
      col <- col + nchar(sep) + nchar(piece)
    } else if (!empty && inner + nchar(piece) <= width) {
      out <- c(out, line)
      line <- paste0(strrep(" ", inner), piece)
      col <- inner + nchar(piece)
    } else {
      if (!empty) {
        out <- c(out, line)
        line <- strrep(" ", inner)
        col <- inner
      }
      sub <- code_format(args[[k]], col + nchar(lead), width)
      sub[1L] <- paste0(lead, sub[1L])
      sub[length(sub)] <- paste0(sub[length(sub)], tail)
      line <- paste0(line, sub[1L])
      if (length(sub) > 1L) {
        out <- c(out, line, sub[-c(1L, length(sub))])
        line <- sub[length(sub)]
        col <- nchar(line)
      } else {
        col <- col + nchar(sub[1L])
      }
    }
    empty <- FALSE
  }
  out <- c(out, line)
  if (hanging) {
    out <- c(out, ")")
  }
  out
}
