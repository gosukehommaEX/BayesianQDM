test_that("report_data() describes the data as a noun phrase", {
  expect_identical(report_data(data_binary(12, 9)),
                   "9 responders among 12 subjects")
  expect_identical(report_data(data_continuous(20, 1.5, 2)),
                   "20 subjects with mean 1.5 and standard deviation 2")
  expect_identical(report_data(data_binary2(c(1, 2, 3, 4))),
                   paste("10 subjects with the counts (1, 2, 3, 4) of the",
                         "responses 00, 01, 10 and 11"))
  expect_identical(report_data(data_tte(4, 30)),
                   "4 events in a total observed time of 30")
  expect_error(report_data(structure(list(), class = "qdm_data")),
               "Unknown data class.", fixed = TRUE)
})

test_that("report_data() formats each number separately", {
  d <- data_continuous2(12, c(18, 2.5), matrix(c(3600, 120, 120, 25), 2L))
  expect_identical(report_data(d),
                   paste("12 subjects with means (18, 2.5) and covariance",
                         "matrix $\\begin{pmatrix} 3600 & 120 \\\\ 120 & 25",
                         "\\end{pmatrix}$"))
  # format() gives the means the same number of decimals
  expect_match(format(d), "mean = (18.0, 2.5)", fixed = TRUE)
  expect_identical(report_data(data_tte(4, 30.5)),
                   "4 events in a total observed time of 30.5")
})

test_that("report_data() writes the singular after a count of one", {
  expect_identical(report_data(data_binary(1, 1)),
                   "1 responder among 1 subject")
  expect_identical(report_data(data_binary(5, 0)),
                   "0 responders among 5 subjects")
  expect_identical(report_data(data_binary2(c(0, 0, 0, 1))),
                   paste("1 subject with the counts (0, 0, 0, 1) of the",
                         "responses 00, 01, 10 and 11"))
  expect_identical(report_data(data_tte(1, 12)),
                   "1 event in a total observed time of 12")
  expect_identical(report_data(data_tte(0, 12)),
                   "0 events in a total observed time of 12")
})
