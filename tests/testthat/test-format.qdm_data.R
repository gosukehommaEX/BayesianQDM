test_that("format.qdm_data() describes each kind of data", {
  expect_identical(format(data_binary(12, 7)), "binary, n = 12, y = 7")
  expect_identical(format(data_continuous(12, 3.2, 2)),
                   "continuous, n = 12, mean = 3.2, sd = 2")
  expect_identical(format(data_binary2(c(4, 2, 3, 6))),
                   paste("two binary endpoints, n = 15,",
                         "counts (00, 01, 10, 11) = (4, 2, 3, 6)"))
  expect_identical(format(data_continuous2(12, c(1, 2), diag(2))),
                   paste("two continuous endpoints, n = 12, mean = (1, 2),",
                         "cov = (1, 0; 0, 1)"))
})
