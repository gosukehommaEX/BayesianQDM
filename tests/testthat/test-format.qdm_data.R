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

test_that("format.qdm_data() honours digits and rejects unknown classes", {
  expect_identical(format(data_continuous(10, 1 / 3, 2 / 3), digits = 2),
                   "continuous, n = 10, mean = 0.33, sd = 0.67")
  fake <- structure(list(), class = c("qdm_data_other", "qdm_data"))
  expect_error(format(fake), "Unknown data class.")
})

test_that("format.qdm_data() describes time-to-event data", {
  expect_identical(format(data_tte(5, 40.25)),
                   "time-to-event, events = 5, exposure = 40.25")
})
