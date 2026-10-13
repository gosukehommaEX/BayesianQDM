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

test_that("borrow() accepts every kind of data", {
  for (d in list(data_continuous(20, 1, 2), data_binary2(c(1, 2, 3, 4)),
                 data_continuous2(20, c(0, 1), diag(2)))) {
    expect_identical(borrow(d, 0.3)$data, d)
  }
  expect_error(borrow(data_binary(20, 8), -0.1), "'a0' must be at least 0")
  expect_error(borrow(data_binary(20, 8), c(0.1, 0.2)),
               "'a0' must have length 1.")
})
