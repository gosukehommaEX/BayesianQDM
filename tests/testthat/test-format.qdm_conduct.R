test_that("format.qdm_conduct() describes the conduct of a study", {
  expect_identical(format(conduct_tte(accrual_time = c(0, 12),
                                      cut_events = 20)),
                   c("Enrollment: uniform from 0 to 12", "Dropout: none",
                     "Analysis: at event 20"))
  expect_identical(
    format(conduct_tte(accrual_time = c(0, 5, 10, 15),
                       accrual_rate = c(0.4, 2, 4.5),
                       dropout_hazard = 0.00226, cut_time = 15.44)),
    c(paste("Enrollment: piecewise uniform on (0, 5, 10, 15) with rates",
            "(0.4, 2.0, 4.5)"),
      "Dropout: hazard 0.00226", "Analysis: at time 15.44")
  )
  expect_identical(
    format(conduct_tte(accrual_time = c(0, 3), dropout_time = c(0, 6),
                       dropout_hazard = c(0.01, 0.05), cut_time = 9))[2],
    "Dropout: hazards (0.01, 0.05) from times (0, 6)"
  )
})
