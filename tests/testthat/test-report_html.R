html_blocks <- function() {
  list(report_block("heading", text = "H"),
       report_block("par", text = "P $x < 1$ & more"),
       report_block("list", items = c("a", "b")),
       report_block("table", header = c("", "B"),
                    rows = matrix(c("r1", "1", "r2", "2"), nrow = 2L,
                                  byrow = TRUE)),
       report_block("pre", lines = c("x", "  y")))
}

test_that("report_html() writes the blocks as HTML", {
  html <- report_html("T & <1>", html_blocks(), "N")
  expect_identical(html[1:2], c("<!DOCTYPE html>", "<html lang=\"en\">"))
  expect_identical(html[length(html) - 0:2],
                   c("</html>", "</body>", "<p class=\"note\">N</p>"))
  expect_true("<title>T &amp; &lt;1&gt;</title>" %in% html)
  expect_true("<h1>T &amp; &lt;1&gt;</h1>" %in% html)
  body <- html[seq(match("<body>", html) + 2L, length(html) - 3L)]
  expect_identical(
    body,
    c("<h2>H</h2>", "<p>P $x &lt; 1$ &amp; more</p>", "<ul>", "<li>a</li>",
      "<li>b</li>", "</ul>", "<table>", "<tr><th></th><th>B</th></tr>",
      "<tr><td>r1</td><td>1</td></tr>", "<tr><td>r2</td><td>2</td></tr>",
      "</table>", "<pre>", "x", "  y", "</pre>")
  )
})

test_that("report_html() loads MathJax with dollar signs for inline formulas", {
  html <- report_html("T", html_blocks(), "N")
  expect_true(
    "window.MathJax = {tex: {inlineMath: [['$', '$'], ['\\\\(', '\\\\)']]}};" %in%
      html
  )
  expect_true(any(grepl("npm/mathjax@3.2.2/es5/tex-chtml.js", html,
                        fixed = TRUE)))
})

test_that("report_html() closes every element that it opens", {
  html <- paste(report_html("T", html_blocks(), "N"), collapse = "\n")
  count <- function(pattern) {
    lengths(regmatches(html, gregexpr(pattern, html)))
  }
  for (tag in c("html", "head", "title", "script", "style", "body", "h1",
                "h2", "p", "ul", "li", "table", "tr", "th", "td", "pre")) {
    expect_identical(count(sprintf("<%s[ >]", tag)),
                     count(sprintf("</%s>", tag)), info = tag)
  }
  # A missing closing tag is detected
  broken <- sub("</ul>", "", html, fixed = TRUE)
  expect_false(identical(
    lengths(regmatches(broken, gregexpr("<ul[ >]", broken))),
    lengths(regmatches(broken, gregexpr("</ul>", broken)))
  ))
})

test_that("report_html() gives lines without names", {
  rows <- matrix(c("r1", "1"), nrow = 1L, dimnames = list("row", c("a", "b")))
  blocks <- list(report_block("table", header = c("a", "b"), rows = rows),
                 report_block("pre", lines = c(first = "x")),
                 report_block("list", items = c(first = "y")))
  expect_null(names(report_html("T", blocks, "N")))
})
