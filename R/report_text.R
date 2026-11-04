#' Report as plain text
#'
#' Writes the blocks of a report as plain text in Markdown format: headings
#' with \code{#}, paragraphs, lists with \code{-}, tables with \code{|} and
#' preformatted text indented by four spaces, with the formulas in LaTeX.
#'
#' @param title Title of the report.
#' @param blocks List of blocks (see \code{report_block()}).
#' @param note Last line of the report.
#'
#' @return A character vector, one element per line.
#' @noRd
report_text <- function(title, blocks, note) {
  table <- function(header, rows) {
    # rbind() names the first row "header" and apply() would pass the row
    # names on to the lines, so the names are removed
    cells <- unname(rbind(header, rows))
    lines <- apply(cells, 1L, function(r) {
      paste0("| ", paste(r, collapse = " | "), " |")
    })
    c(lines[1L],
      paste0("|", paste(rep("---", ncol(cells)), collapse = "|"), "|"),
      lines[-1L])
  }
  out <- c(paste("#", title), "")
  for (b in blocks) {
    out <- c(out, switch(b$type,
                         heading = paste("##", b$text),
                         par = b$text,
                         list = paste("-", b$items),
                         table = table(b$header, b$rows),
                         pre = paste0("    ", b$lines)), "")
  }
  unname(c(out, note))
}
