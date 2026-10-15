test_that("check_flag() accepts TRUE and FALSE only", {
  expect_silent(check_flag(TRUE, "f"))
  expect_silent(check_flag(FALSE, "f"))
  expect_error(check_flag(NA, "f"), "'f' must be TRUE or FALSE.")
  expect_error(check_flag(c(TRUE, FALSE), "f"), "'f' must be TRUE or FALSE.")
  expect_error(check_flag(1, "f"), "'f' must be TRUE or FALSE.")
})
