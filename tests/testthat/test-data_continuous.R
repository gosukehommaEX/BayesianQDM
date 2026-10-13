test_that("data_continuous() stores the summary statistics", {
  d <- data_continuous(n = 12, mean = -1.5, sd = 2)
  expect_s3_class(d, c("qdm_data_continuous", "qdm_data"), exact = TRUE)
  expect_identical(unclass(d), list(n = 12L, mean = -1.5, sd = 2))
})

test_that("data_continuous() rejects invalid summaries", {
  expect_error(data_continuous(1, 0, 1), "'n' must be at least 2.")
  expect_error(data_continuous(10, 0, 0), "'sd' must be greater than 0.")
  expect_error(data_continuous(10, Inf, 1), "'mean' must be finite.")
})

test_that("data_continuous() rejects non-numeric and vector input", {
  expect_error(data_continuous(10, "1", 1), "'mean' must be numeric.")
  expect_error(data_continuous(10, c(1, 2), 1), "'mean' must have length 1.")
  expect_error(data_continuous(10.5, 1, 1), "whole numbers")
})
