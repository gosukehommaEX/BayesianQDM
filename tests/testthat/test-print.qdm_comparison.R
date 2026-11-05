test_that("print.qdm_comparison() prints the method and the table", {
  cmp <- qdm_compare(list(a = ref_design("D1"), b = ref_design("D2")), scen_c)
  out <- capture.output(res <- withVisible(print(cmp)))
  expect_false(res$visible)
  expect_identical(res$value, cmp)
  expect_identical(out[1:2], c("Operating characteristics of 2 designs",
                               "Method: exact enumeration of the outcomes"))
  y <- cmp
  attr(y, "designs") <- NULL
  attr(y, "method") <- NULL
  class(y) <- "data.frame"
  table <- capture.output(print(y, digits = 4L, row.names = FALSE))
  expect_identical(out[-(1:2)], table)
  # One design, and designs with different methods
  one <- qdm_compare(list(ref_design("D1")), scen_c[1, ])
  expect_identical(capture.output(print(one))[1L],
                   "Operating characteristics of 1 design")
  mixed <- cmp
  attr(mixed, "method") <- c(a = "m1", b = "m2")
  expect_identical(capture.output(print(mixed))[2:3],
                   c("  a: m1", "  b: m2"))
  none <- cmp
  attr(none, "method") <- NULL
  expect_identical(capture.output(print(none))[-1L], table)
})
