test_that("report_data() describes the data without the type of endpoint", {
  expect_identical(report_data(data_binary(20, 8)), "n = 20, y = 8")
  expect_identical(report_data(data_continuous(20, 1.5, 2)),
                   "n = 20, mean = 1.5, sd = 2")
  expect_identical(report_data(data_binary2(c(1, 2, 3, 4))),
                   "n = 10, counts (00, 01, 10, 11) = (1, 2, 3, 4)")
  expect_identical(report_data(data_tte(4, 30)),
                   "events = 4, total observed time = 30")
})

test_that("report_data() formats each number separately", {
  d <- data_continuous2(12, c(18, 2.5), matrix(c(3600, 120, 120, 25), 2L))
  expect_identical(report_data(d),
                   "n = 12, mean = (18, 2.5), cov = (3600, 120; 120, 25)")
  # format() gives the means the same number of decimals
  expect_match(format(d), "mean = (18.0, 2.5)", fixed = TRUE)
})

test_that("report_data() describes the data as format() does", {
  for (d in list(data_binary(20, 8), data_continuous(20, 1.5, 2),
                 data_binary2(c(1, 2, 3, 4)),
                 data_continuous2(20, c(-0.5, -1.5),
                                  matrix(c(24, 2, 2, 10), 2L)))) {
    expect_identical(report_data(d), sub("^[^,]*, ", "", format(d)))
  }
  expect_identical(report_data(data_tte(4, 30.5)),
                   sub("exposure", "total observed time",
                       sub("^[^,]*, ", "", format(data_tte(4, 30.5))),
                       fixed = TRUE))
  expect_error(report_data(structure(list(), class = "qdm_data")),
               "Unknown data class.", fixed = TRUE)
})
