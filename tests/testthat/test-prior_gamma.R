test_that("prior_gamma() stores the shape and rate", {
  p <- prior_gamma(0.1, 2)
  expect_s3_class(p, c("qdm_prior_gamma", "qdm_prior"), exact = TRUE)
  expect_identical(p$shape, 0.1)
  expect_identical(p$rate, 2)
})

test_that("prior_gamma() rejects invalid parameters", {
  expect_error(prior_gamma(0, 1), "'shape' must be greater than 0.")
  expect_error(prior_gamma(1, -1), "'rate' must be greater than 0.")
  expect_error(prior_gamma(c(1, 2), 1), "'shape' must have length 1.")
})
