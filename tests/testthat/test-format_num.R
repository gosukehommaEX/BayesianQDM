test_that("format_num() formats numbers as a comma-separated string", {
  expect_identical(format_num(c(0.5, 12)), "0.5, 12.0")
  expect_identical(format_num(1 / 3, digits = 3), "0.333")
})
