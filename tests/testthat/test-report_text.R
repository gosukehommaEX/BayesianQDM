report_test_blocks <- function() {
  list(report_block("heading", text = "H"),
       report_block("par", text = "P $x < 1$ & more"),
       report_block("list", items = c("a", "b")),
       report_block("table", header = c("", "B"),
                    rows = matrix(c("r1", "1", "r2", "2"), nrow = 2L,
                                  byrow = TRUE)),
       report_block("pre", lines = c("x", "  y")))
}

test_that("report_text() writes the blocks as Markdown text", {
  expect_identical(
    report_text("T", report_test_blocks(), "N"),
    c("# T", "", "## H", "", "P $x < 1$ & more", "", "- a", "- b", "",
      "|  | B |", "|---|---|", "| r1 | 1 |", "| r2 | 2 |", "", "    x",
      "      y", "", "N")
  )
})

test_that("report_text() gives lines without names", {
  rows <- matrix(c("r1", "1"), nrow = 1L, dimnames = list("row", c("a", "b")))
  blocks <- list(report_block("table", header = c("a", "b"), rows = rows),
                 report_block("pre", lines = c(first = "x")))
  expect_null(names(report_text("T", blocks, "N")))
  expect_null(names(report_text("T", report_test_blocks(), "N")))
})
