test_that("check_numeric() accepts valid input and returns it invisibly", {
  expect_invisible(check_numeric(1.5, "x"))
  expect_identical(check_numeric(c(1, 2), "x", len = 2L), c(1, 2))
  expect_silent(check_numeric(Inf, "x", allow_inf = TRUE))
  expect_silent(check_numeric(3, "x", integer = TRUE, lower = 1))
})

test_that("check_numeric() rejects invalid input with informative messages", {
  expect_error(check_numeric("a", "x"), "'x' must be numeric.")
  expect_error(check_numeric(numeric(0), "x"), "'x' must be numeric.")
  expect_error(check_numeric(c(1, 2), "x"), "'x' must have length 1.")
  expect_error(check_numeric(NA_real_, "x"), "'x' must not contain missing")
  expect_error(check_numeric(Inf, "x"), "'x' must be finite.")
  expect_error(check_numeric(1.5, "x", integer = TRUE), "whole numbers")
  expect_error(check_numeric(0, "x", lower = 0, lower_open = TRUE),
               "'x' must be greater than 0.")
  expect_error(check_numeric(-1, "x", lower = 0), "'x' must be at least 0.")
  expect_error(check_numeric(1, "x", lower = 0, upper = 1, upper_open = TRUE),
               "'x' must be at least 0 and less than 1.")
})
