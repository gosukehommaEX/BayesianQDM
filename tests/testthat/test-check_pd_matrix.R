test_that("check_pd_matrix() accepts symmetric positive definite matrices", {
  expect_silent(check_pd_matrix(diag(2), "m"))
  expect_silent(check_pd_matrix(matrix(c(4, 1.2, 1.2, 2), 2), "m"))
})

test_that("check_pd_matrix() rejects other matrices", {
  expect_error(check_pd_matrix(diag(3), "m"), "'m' must be a numeric 2 x 2")
  expect_error(check_pd_matrix(c(1, 0, 0, 1), "m"), "'m' must be a numeric")
  expect_error(check_pd_matrix(matrix(c(1, 0.5, 0, 1), 2), "m"),
               "'m' must be symmetric.")
  expect_error(check_pd_matrix(matrix(c(1, 2, 2, 1), 2), "m"),
               "'m' must be positive definite.")
  expect_error(check_pd_matrix(matrix(c(1, NA, NA, 1), 2), "m"),
               "'m' must contain finite values only.")
})
