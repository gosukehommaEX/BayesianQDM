test_that("prior_beta() stores the shape parameters", {
  p <- prior_beta(0.5, 2)
  expect_s3_class(p, c("qdm_prior_beta", "qdm_prior"), exact = TRUE)
  expect_identical(p$a, 0.5)
  expect_identical(p$b, 2)
})

test_that("prior_beta() rejects non-positive shapes", {
  expect_error(prior_beta(0, 1), "'a' must be greater than 0.")
  expect_error(prior_beta(1, -1), "'b' must be greater than 0.")
  expect_error(prior_beta(c(1, 2), 1), "'a' must have length 1.")
})

test_that("prior_beta() rejects missing and non-numeric shapes", {
  expect_error(prior_beta(NA_real_, 1), "'a' must not contain missing")
  expect_error(prior_beta(1, "2"), "'b' must be numeric.")
  expect_error(prior_beta(Inf, 1), "'a' must be finite.")
})
