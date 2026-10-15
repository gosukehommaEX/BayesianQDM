test_that("check_prior_control() requires NULL when the prior is not used", {
  pr <- prior_beta(1, 1)
  expect_null(check_prior_control(NULL, "fixed", check_prior_beta))
  expect_error(check_prior_control(pr, "fixed", check_prior_beta),
               "'prior_c' must be NULL for a fixed control, whose value is")
  expect_error(check_prior_control(prior_vague(), "hypothetical",
                                   check_prior_continuous,
                                   unused = c("hypothetical", "fixed")),
               "'prior_c' must be NULL for a hypothetical control, whose")
})

test_that("check_prior_control() checks the prior otherwise", {
  pr <- prior_beta(1, 1)
  expect_identical(check_prior_control(pr, "concurrent", check_prior_beta),
                   pr)
  expect_error(check_prior_control(NULL, "external", check_prior_beta),
               "'prior_c' must be created by prior_beta().")
  expect_error(check_prior_control(prior_vague(), "hypothetical",
                                   check_prior_beta),
               "'prior_c' must be created by prior_beta().")
})
