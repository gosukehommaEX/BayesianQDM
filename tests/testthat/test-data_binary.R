test_that("data_binary() stores integer counts", {
  d <- data_binary(n = 12, y = 7)
  expect_s3_class(d, c("qdm_data_binary", "qdm_data"), exact = TRUE)
  expect_identical(d$n, 12L)
  expect_identical(d$y, 7L)
  expect_identical(data_binary(5, 0)$y, 0L)
  expect_identical(data_binary(5, 5)$y, 5L)
})

test_that("data_binary() rejects invalid counts", {
  expect_error(data_binary(0, 0), "'n' must be at least 1.")
  expect_error(data_binary(10, 11), "'y' must not exceed 'n'.")
  expect_error(data_binary(10, 2.5), "whole numbers")
  expect_error(data_binary(10, -1), "'y' must be at least 0.")
})

test_that("data_binary() accepts whole numbers stored as doubles only", {
  expect_identical(data_binary(12.0, 3.0)$n, 12L)
  expect_error(data_binary(c(10, 12), 3), "'n' must have length 1.")
  expect_error(data_binary(10.5, 3), "whole numbers")
})
