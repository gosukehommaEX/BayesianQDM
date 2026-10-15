test_that("check_prior_beta() accepts Beta priors only", {
  pr <- prior_beta(1, 2)
  expect_identical(check_prior_beta(pr, "p"), pr)
  expect_invisible(check_prior_beta(pr, "p"))
  expect_error(check_prior_beta(prior_vague(), "prior_x"),
               "'prior_x' must be created by prior_beta().")
  expect_error(check_prior_beta(NULL, "prior_x"), "'prior_x' must be")
})
