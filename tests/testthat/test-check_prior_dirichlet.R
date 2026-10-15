test_that("check_prior_dirichlet() accepts only Dirichlet priors", {
  p <- prior_dirichlet(rep(1, 4))
  expect_identical(check_prior_dirichlet(p, "prior_t"), p)
  expect_error(check_prior_dirichlet(prior_vague(), "prior_c"),
               "'prior_c' must be created by prior_dirichlet\\(\\).")
})
