test_that("report_block() stores the type and the elements of a block", {
  expect_identical(report_block("par", text = "A"),
                   list(type = "par", text = "A"))
  rows <- matrix(c("x", "y"), nrow = 1L)
  expect_identical(report_block("table", header = c("a", "b"), rows = rows),
                   list(type = "table", header = c("a", "b"), rows = rows))
})
