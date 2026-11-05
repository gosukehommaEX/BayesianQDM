#' Call as lines of R code
#'
#' Writes a call on one line if it fits in \code{width} characters from the
#' column \code{indent}. Otherwise the arguments are placed on lines aligned
#' after the opening parenthesis, as many on a line as fit, and an argument
#' that does not fit on a line is itself written over several lines. If the
#' aligned arguments do not fit (for example because the call starts far to
#' the right), or with \code{hanging = TRUE}, the arguments start on the next
#' line, indented by two spaces from the line on which the call starts, and
#' the closing parenthesis is on a line of its own, as in the examples of the
#' package.
#'
#' @param node A character string of R code or a call from
#'   \code{code_node()}.
#' @param indent Column at which the call starts.
#' @param width Largest width of a line.
#' @param hanging Whether to start the arguments on the next line (for a
#'   call that starts a line).
#' @param base Indentation of the line on which the call starts.
#' @param tail Number of characters that the caller writes after the call on
#'   its last line (such as a comma or a closing parenthesis).
#' @param fallback Whether the arguments may start on the next line when the
#'   aligned arguments do not fit. The arguments of a call whose own
#'   arguments are aligned are aligned too, so that the line breaks are made
#'   in the outermost call that needs them.
#'
#' @return A character vector of lines. The first line has no indentation
#'   (the caller writes what comes before it); the other lines start with
#'   the spaces of their columns.
#' @noRd
code_format <- function(node, indent = 0L, width = 80L, hanging = FALSE,
                        base = 0L, tail = 0L, fallback = TRUE) {
  flat <- code_flat(node)
  if (is.character(node) || indent + nchar(flat) + tail <= width) {
    return(flat)
  }
  args <- node$args
  n <- length(args)
  nms <- names(args)
  if (is.null(nms)) {
    nms <- rep("", n)
  }
  place_args <- function(hang, sub_fallback) {
    open <- paste0(node$fun, "(")
    if (hang) {
      inner <- base + 2L
      out <- open
      line <- strrep(" ", inner)
      line_base <- inner
    } else {
      inner <- indent + nchar(open)
      out <- character(0)
      line <- open
      line_base <- base
    }
    col <- inner
    empty <- TRUE
    for (k in seq_len(n)) {
      lead <- if (nzchar(nms[k])) paste(nms[k], "= ") else ""
      own <- if (k < n) "," else if (hang) "" else ")"
      # The characters of the caller follow the last aligned argument
      extra <- if (k == n && !hang) tail else 0L
      piece <- paste0(lead, code_flat(args[[k]]), own)
      sep <- if (empty) "" else " "
      if (col + nchar(sep) + nchar(piece) + extra <= width) {
        line <- paste0(line, sep, piece)
        col <- col + nchar(sep) + nchar(piece)
      } else if (!empty && inner + nchar(piece) + extra <= width) {
        out <- c(out, line)
        line <- paste0(strrep(" ", inner), piece)
        col <- inner + nchar(piece)
        line_base <- inner
      } else {
        if (!empty) {
          out <- c(out, line)
          line <- strrep(" ", inner)
          col <- inner
          line_base <- inner
        }
        sub <- code_format(args[[k]], col + nchar(lead), width,
                           base = line_base, tail = nchar(own) + extra,
                           fallback = sub_fallback)
        sub[1L] <- paste0(lead, sub[1L])
        sub[length(sub)] <- paste0(sub[length(sub)], own)
        line <- paste0(line, sub[1L])
        if (length(sub) > 1L) {
          out <- c(out, line, sub[-c(1L, length(sub))])
          line <- sub[length(sub)]
          col <- nchar(line)
          line_base <- nchar(line) - nchar(trimws(line, which = "left"))
        } else {
          col <- col + nchar(sub[1L])
        }
      }
      empty <- FALSE
    }
    out <- c(out, line)
    if (hang) c(out, paste0(strrep(" ", base), ")")) else out
  }
  if (!hanging) {
    out <- place_args(FALSE, FALSE)
    widths <- c(indent + nchar(out[1L]), nchar(out[-1L]))
    widths[length(widths)] <- widths[length(widths)] + tail
    if (!fallback || max(widths) <= width) {
      return(out)
    }
  }
  place_args(TRUE, TRUE)
}
