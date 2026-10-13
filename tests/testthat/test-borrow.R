test_that("borrow() stores the external data and the weight", {
  d <- data_binary(20, 8)
  b <- borrow(d, a0 = 0.5)
  expect_s3_class(b, "qdm_borrow", exact = TRUE)
  expect_identical(b$data, d)
  expect_identical(b$a0, 0.5)
  expect_identical(borrow(d, 0)$a0, 0)
  expect_identical(borrow(d, 1)$a0, 1)
})

test_that("borrow() rejects invalid input", {
  expect_error(borrow(list(n = 20, y = 8), 0.5), "created by a data_")
  expect_error(borrow(data_binary(20, 8), 1.2), "'a0' must be at least 0")
})
