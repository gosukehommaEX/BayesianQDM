test_that("conduct_tte() stores the conduct of a study", {
  co <- conduct_tte(accrual_time = c(0, 5, 10, 15),
                    accrual_rate = c(0.4, 2, 4.5), dropout_hazard = 0.00226,
                    cut_time = 15.44)
  expect_s3_class(co, "qdm_conduct", exact = TRUE)
  expect_identical(co$accrual_time, c(0, 5, 10, 15))
  # Interval probabilities proportional to rate times length (2, 10, 22.5)
  expect_equal(co$accrual_cum, cumsum(c(2, 10, 22.5)) / 34.5)
  expect_identical(co$dropout_time, 0)
  expect_identical(co$dropout_hazard, 0.00226)
  expect_identical(co$cut_time, 15.44)
  expect_null(co$cut_events)
  co2 <- conduct_tte(accrual_time = c(0, 12), cut_events = 20)
  expect_identical(co2$accrual_rate, 1)
  expect_identical(co2$accrual_cum, 1)
  expect_identical(co2$cut_events, 20L)
  expect_null(co2$cut_time)
  # Intervals without enrollment
  co3 <- conduct_tte(accrual_time = c(0, 2, 4), accrual_rate = c(1, 0),
                     cut_time = 10)
  expect_identical(co3$accrual_cum, c(1, 1))
})

test_that("conduct_tte() validates its arguments", {
  expect_error(conduct_tte(accrual_time = 5, cut_time = 1),
               "'accrual_time' must be an increasing vector of at least two")
  expect_error(conduct_tte(accrual_time = c(1, 5), cut_time = 1),
               "starts at 0")
  expect_error(conduct_tte(accrual_time = c(0, 5, 5), cut_time = 1),
               "'accrual_time' must be an increasing vector")
  expect_error(conduct_tte(accrual_time = c(0, 5), accrual_rate = c(1, 2),
                           cut_time = 1),
               "'accrual_rate' must have length 1.")
  expect_error(conduct_tte(accrual_time = c(0, 5), accrual_rate = 0,
                           cut_time = 1),
               "'accrual_rate' must have a positive sum.")
  expect_error(conduct_tte(accrual_time = c(0, 5), dropout_time = c(1, 2),
                           dropout_hazard = c(0, 0), cut_time = 1),
               "'dropout_time' must be an increasing vector that starts at 0.")
  expect_error(conduct_tte(accrual_time = c(0, 5), dropout_time = c(0, 2),
                           dropout_hazard = 0.1, cut_time = 1),
               "'dropout_hazard' must have length 2.")
  expect_error(conduct_tte(accrual_time = c(0, 5)),
               "Specify exactly one of 'cut_time' and 'cut_events'.")
  expect_error(conduct_tte(accrual_time = c(0, 5), cut_time = 1,
                           cut_events = 3),
               "Specify exactly one of")
  expect_error(conduct_tte(accrual_time = c(0, 5), cut_time = 0),
               "'cut_time' must be greater than 0.")
  expect_error(conduct_tte(accrual_time = c(0, 5), cut_events = 0.5),
               "'cut_events' must contain whole numbers only.")
})
