test_that("check_prior_continuous2() accepts NIW and vague priors", {
  expect_invisible(check_prior_continuous2(prior_vague(), "p"))
  pr <- prior_niw(c(0, 0), 1, 3, diag(2))
  expect_identical(check_prior_continuous2(pr, "p"), pr)
  expect_error(check_prior_continuous2(prior_nix(0, 1, 2, 3), "prior_c"),
               "'prior_c' must be created by prior_niw\\(\\) or prior_vague")
  expect_error(check_prior_continuous2(NULL, "prior_t"), "'prior_t' must be")
})
