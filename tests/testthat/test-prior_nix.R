test_that("prior_nix() stores the hyperparameters", {
  p <- prior_nix(mu0 = 1, kappa0 = 2, nu0 = 3, sigma0 = 4)
  expect_s3_class(p, c("qdm_prior_nix", "qdm_prior"), exact = TRUE)
  expect_identical(unclass(p),
                   list(mu0 = 1, kappa0 = 2, nu0 = 3, sigma0 = 4))
})

test_that("prior_nix() rejects invalid hyperparameters", {
  expect_error(prior_nix(NA, 1, 1, 1), "'mu0' must not contain missing")
  expect_error(prior_nix(0, 0, 1, 1), "'kappa0' must be greater than 0.")
  expect_error(prior_nix(0, 1, 0, 1), "'nu0' must be greater than 0.")
  expect_error(prior_nix(0, 1, 1, -2), "'sigma0' must be greater than 0.")
})
