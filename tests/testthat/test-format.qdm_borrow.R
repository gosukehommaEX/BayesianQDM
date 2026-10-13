test_that("format.qdm_borrow() describes the data and the weight", {
  expect_identical(format(borrow(data_binary(20, 8), 0.5)),
                   "binary, n = 20, y = 8, a0 = 0.5")
})

test_that("format.qdm_borrow() works for every kind of external data", {
  expect_identical(format(borrow(data_continuous(20, 1.5, 2), 1)),
                   "continuous, n = 20, mean = 1.5, sd = 2, a0 = 1")
  expect_identical(format(borrow(data_binary2(c(1, 2, 3, 4)), 0)),
                   paste("two binary endpoints, n = 10,",
                         "counts (00, 01, 10, 11) = (1, 2, 3, 4), a0 = 0"))
})
