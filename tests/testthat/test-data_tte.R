test_that("data_tte() stores the events and the total observed time", {
  d <- data_tte(events = 5, exposure = 40.5)
  expect_s3_class(d, c("qdm_data_tte", "qdm_data"), exact = TRUE)
  expect_identical(d$events, 5L)
  expect_identical(d$exposure, 40.5)
  expect_identical(data_tte(0, 1)$events, 0L)
})

test_that("data_tte() rejects invalid data", {
  expect_error(data_tte(-1, 10), "'events' must be at least 0.")
  expect_error(data_tte(2.5, 10), "whole numbers")
  expect_error(data_tte(2, 0), "'exposure' must be greater than 0.")
})
