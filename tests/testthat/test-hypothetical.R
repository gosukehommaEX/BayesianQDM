test_that("hypothetical() stores data for binary endpoints", {
  d <- data_binary(15, 5)
  h <- hypothetical(data = d)
  expect_s3_class(h, "qdm_hypothetical", exact = TRUE)
  expect_identical(h$type, "data")
  expect_identical(h$data, d)
  expect_null(h$mean)
})

test_that("hypothetical() stores location and scale for continuous endpoints", {
  h1 <- hypothetical(mean = 1, r = 2)
  expect_identical(h1$type, "location")
  expect_identical(h1$mean, 1)
  expect_identical(h1$r, 2)
  expect_identical(hypothetical(mean = c(1, 2), r = 1)$mean, c(1, 2))
})

test_that("hypothetical() rejects incomplete or conflicting input", {
  expect_error(hypothetical(), "Specify either 'data'")
  expect_error(hypothetical(data = data_binary(15, 5), mean = 1),
               "Specify either 'data'")
  expect_error(hypothetical(mean = 1), "Both 'mean' and 'r' must be")
  expect_error(hypothetical(data = data_continuous(10, 0, 1)),
               "created by data_binary\\(\\) or data_binary2\\(\\)")
  expect_error(hypothetical(mean = 1:3, r = 1), "'mean' must have length 1")
  expect_error(hypothetical(mean = 1, r = 0), "'r' must be greater than 0.")
})
