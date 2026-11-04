test_that("report_conduct() formats each number of the conduct separately", {
  expect_identical(
    report_conduct(conduct_example),
    c(paste("Enrollment: piecewise uniform on (0, 5, 10, 15) with rates",
            "(0.4, 2, 4.5)"),
      "Dropout: hazard 0.00226", "Analysis: at time 15.44")
  )
  # format() gives the rates the same number of decimals
  expect_identical(format(conduct_example)[1L],
                   paste("Enrollment: piecewise uniform on (0, 5, 10, 15)",
                         "with rates (0.4, 2.0, 4.5)"))
})

test_that("report_conduct() describes the conduct as format() does", {
  x <- conduct_tte(accrual_time = c(0, 12), cut_events = 20)
  expect_identical(report_conduct(x),
                   c("Enrollment: uniform from 0 to 12", "Dropout: none",
                     "Analysis: at event 20"))
  expect_identical(report_conduct(x), format(x))
  y <- conduct_tte(accrual_time = c(0, 3), dropout_time = c(0, 6),
                   dropout_hazard = c(0.01, 0.05), cut_time = 9)
  expect_identical(report_conduct(y)[2L],
                   "Dropout: hazards (0.01, 0.05) from times (0, 6)")
  expect_identical(report_conduct(y), format(y))
  expect_identical(report_conduct(conduct_complete), format(conduct_complete))
})
