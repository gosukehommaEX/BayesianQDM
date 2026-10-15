test_that("check_fixed_c() returns NULL or the values", {
  expect_null(check_fixed_c(NULL, 1L, probability = TRUE))
  expect_identical(check_fixed_c(0.2, 1L, probability = TRUE), 0.2)
  expect_identical(check_fixed_c(c(0L, 1L), 2L, probability = TRUE), c(0, 1))
  expect_identical(check_fixed_c(-10, 1L, probability = FALSE), -10)
})

test_that("check_fixed_c() checks the length and the range", {
  expect_error(check_fixed_c(c(0.2, 0.3), 1L, probability = TRUE),
               "'fixed_c' must have length 1.")
  expect_error(check_fixed_c(0.2, 2L, probability = FALSE),
               "'fixed_c' must have length 2.")
  expect_error(check_fixed_c(1.2, 1L, probability = TRUE),
               "'fixed_c' must be at least 0 and at most 1.")
  expect_error(check_fixed_c(Inf, 1L, probability = FALSE),
               "'fixed_c' must be finite.")
  expect_error(check_fixed_c("a", 1L, probability = FALSE),
               "'fixed_c' must be numeric.")
})
