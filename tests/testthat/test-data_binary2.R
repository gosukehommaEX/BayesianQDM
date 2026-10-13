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
