test_that("print.qdm_borrow() prints the description and returns invisibly", {
  b <- borrow(data_binary(20, 8), 0.5)
  expect_output(res <- withVisible(print(b)),
                "External data: binary, n = 20, y = 8, a0 = 0.5")
  expect_false(res$visible)
  expect_identical(res$value, b)
})
