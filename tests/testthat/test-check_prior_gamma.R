test_that("check_prior_gamma() accepts only Gamma priors", {
  p <- prior_gamma(1, 1)
  expect_identical(check_prior_gamma(p, "prior_t"), p)
  expect_error(check_prior_gamma(prior_beta(1, 1), "prior_c"),
               "'prior_c' must be created by prior_gamma\\(\\).")
})
