test_that("report_data() describes the data without the type of endpoint", {
  expect_identical(report_data(data_binary(20, 8)), "n = 20, y = 8")
  expect_identical(report_data(data_continuous(20, 1.5, 2)),
                   "n = 20, mean = 1.5, sd = 2")
  expect_identical(report_data(data_binary2(c(1, 2, 3, 4))),
                   "n = 10, counts (00, 01, 10, 11) = (1, 2, 3, 4)")
  expect_identical(report_data(data_tte(4, 30)),
                   "events = 4, total observed time = 30")
})
