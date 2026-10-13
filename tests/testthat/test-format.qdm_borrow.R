test_that("format.qdm_borrow() describes the data and the weight", {
  expect_identical(format(borrow(data_binary(20, 8), 0.5)),
                   "binary, n = 20, y = 8, a0 = 0.5")
})
