test_that("prior_niw() stores the hyperparameters", {
  L <- matrix(c(2, 0.5, 0.5, 1), 2)
  p <- prior_niw(mu0 = c(1, 2), kappa0 = 1, nu0 = 4, Lambda0 = L)
  expect_s3_class(p, c("qdm_prior_niw", "qdm_prior"), exact = TRUE)
  expect_identical(p$mu0, c(1, 2))
  expect_identical(p$Lambda0, L)
})

test_that("prior_niw() rejects invalid hyperparameters", {
  expect_error(prior_niw(1, 1, 4, diag(2)), "'mu0' must have length 2.")
  expect_error(prior_niw(c(0, 0), 1, 1, diag(2)),
               "'nu0' must be greater than 1.")
  expect_error(prior_niw(c(0, 0), 1, 4, matrix(c(1, 2, 2, 1), 2)),
               "'Lambda0' must be positive definite.")
})
