test_that("report_num() formats each number separately", {
  expect_identical(report_num(0.2), "0.2")
  expect_identical(report_num(30L), "30")
  expect_identical(report_num(c(-10, 1.23456)), "-10, 1.235")
  expect_identical(report_num(c(0.1, 0.15)), "0.1, 0.15")
  # format_num() gives the numbers of a vector the same number of decimals
  expect_identical(format_num(c(0.1, 0.15), 4L), "0.10, 0.15")
})
