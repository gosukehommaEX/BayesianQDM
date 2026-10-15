test_that("conduct_args() gives the arguments of the C++ functions", {
  co <- conduct_tte(accrual_time = c(0, 5, 10), accrual_rate = c(1, 3),
                    dropout_time = c(0, 2), dropout_hazard = c(0, 0.1),
                    cut_time = 12)
  expect_identical(conduct_args(co),
                   list(acc_time = c(0, 5, 10), acc_cum = c(0.25, 1),
                        drop_start = c(0, 2), drop_hazard = c(0, 0.1),
                        cut_time = 12, cut_events = 0L))
  co2 <- conduct_tte(accrual_time = c(0, 12), cut_events = 20)
  a2 <- conduct_args(co2)
  expect_identical(a2$cut_time, Inf)
  expect_identical(a2$cut_events, 20L)
})
