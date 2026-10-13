test_that("format.qdm_hypothetical() describes both kinds of control", {
  expect_identical(format(hypothetical(data = data_binary(15, 5))),
                   "data (binary, n = 15, y = 5)")
  expect_identical(format(hypothetical(mean = c(1, 2.5), r = 2)),
                   "mean = (1.0, 2.5), r = 2")
})
