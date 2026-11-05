test_that("print.qdm_code() prints the lines of code", {
  code <- show_code(ref_design("D1"))
  out <- capture.output(v <- withVisible(print(code)))
  expect_identical(out, unclass(code))
  expect_false(v$visible)
  expect_identical(v$value, code)
})
