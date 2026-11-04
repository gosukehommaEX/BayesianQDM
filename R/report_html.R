#' Report as HTML
#'
#' Writes the blocks of a report as a standalone HTML page. The formulas,
#' in LaTeX between dollar signs, are displayed by MathJax, loaded from a
#' content delivery network.
#'
#' @param title Title of the report.
#' @param blocks List of blocks (see \code{report_block()}).
#' @param note Last line of the report.
#'
#' @return A character vector, one element per line.
#' @noRd
report_html <- function(title, blocks, note) {
  esc <- html_escape
  row <- function(cells, tag) {
    paste0("<tr>", paste0("<", tag, ">", esc(cells), "</", tag, ">",
                          collapse = ""), "</tr>")
  }
  body <- character(0)
  for (b in blocks) {
    body <- c(body, switch(
      b$type,
      heading = paste0("<h2>", esc(b$text), "</h2>"),
      par = paste0("<p>", esc(b$text), "</p>"),
      list = c("<ul>", paste0("<li>", esc(b$items), "</li>"), "</ul>"),
      table = c("<table>", row(b$header, "th"),
                apply(b$rows, 1L, row, tag = "td"), "</table>"),
      pre = c("<pre>", esc(b$lines), "</pre>")
    ))
  }
  # unname(): the lines would keep the names of the text of the blocks (as
  # html_escape() keeps them) and of the rows of the tables
  unname(c("<!DOCTYPE html>",
    "<html lang=\"en\">",
    "<head>",
    "<meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
    paste0("<title>", esc(title), "</title>"),
    "<script>",
    "window.MathJax = {tex: {inlineMath: [['$', '$'], ['\\\\(', '\\\\)']]}};",
    "</script>",
    paste0("<script id=\"MathJax-script\" async src=\"https://cdn.jsdelivr",
           ".net/npm/mathjax@3.2.2/es5/tex-chtml.js\"></script>"),
    "<style>",
    paste("body {font-family: Arial, Helvetica, sans-serif; color: #222222;",
          "background: #ffffff; max-width: 52em; margin: 2em auto;",
          "padding: 0 1em; line-height: 1.5;}"),
    "table {border-collapse: collapse; margin: 1em 0;}",
    "th, td {border: 1px solid #999999; padding: 0.25em 0.6em;}",
    "th {background: #eeeeee;}",
    "pre {background: #f5f5f5; padding: 0.75em; overflow-x: auto;}",
    ".note {color: #666666; font-size: 0.9em;}",
    "</style>",
    "</head>",
    "<body>",
    paste0("<h1>", esc(title), "</h1>"),
    body,
    paste0("<p class=\"note\">", esc(note), "</p>"),
    "</body>",
    "</html>"))
}
