test_that("check_conduct() accepts the conduct of a study", {
  co <- conduct_tte(accrual_time = c(0, 12), cut_events = 20)
  expect_identical(check_conduct(co, "conduct", 20L), co)
  expect_null(check_conduct(NULL, "conduct", 20L, allow_null = TRUE))
})

test_that("check_conduct() rejects other objects and too many events", {
  co <- conduct_tte(accrual_time = c(0, 12), cut_events = 20)
  expect_error(check_conduct(NULL, "conduct_future", 20L),
               "'conduct_future' must be created by conduct_tte\\(\\).")
  expect_error(check_conduct(list(), "conduct", 20L),
               "'conduct' must be created by conduct_tte\\(\\).")
  expect_error(check_conduct(co, "conduct", 19L),
               "The analysis of 'conduct' is at event 20, but the study has only 19 subjects.")
})
