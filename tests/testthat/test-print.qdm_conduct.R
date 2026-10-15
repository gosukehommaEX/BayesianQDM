test_that("print.qdm_conduct() prints the description", {
  co <- conduct_tte(accrual_time = c(0, 12), cut_events = 20)
  expect_output(res <- print(co), "Enrollment: uniform from 0 to 12")
  expect_identical(res, co)
})
