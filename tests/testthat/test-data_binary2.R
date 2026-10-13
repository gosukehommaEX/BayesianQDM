test_that("data_binary2() stores named integer counts and their total", {
  d <- data_binary2(c(4, 2, 3, 6))
  expect_s3_class(d, c("qdm_data_binary2", "qdm_data"), exact = TRUE)
  expect_identical(d$counts, c("00" = 4L, "01" = 2L, "10" = 3L, "11" = 6L))
  expect_identical(d$n, 15L)
})

test_that("data_binary2() rejects invalid counts", {
  expect_error(data_binary2(c(1, 2, 3)), "'counts' must have length 4.")
  expect_error(data_binary2(c(0, 0, 0, 0)), "at least one subject")
  expect_error(data_binary2(c(1, -1, 0, 0)), "'counts' must be at least 0.")
})

test_that("data_binary2() ignores input names and checks whole numbers", {
  d <- data_binary2(c(a = 1, b = 0, c = 0, d = 0))
  expect_identical(names(d$counts), c("00", "01", "10", "11"))
  expect_identical(d$n, 1L)
  expect_error(data_binary2(c(1, 0.5, 0, 0)), "whole numbers")
})
