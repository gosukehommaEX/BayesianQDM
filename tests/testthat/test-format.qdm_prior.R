test_that("format.qdm_prior() describes each prior", {
  expect_identical(format(prior_beta(0.5, 0.5)), "Beta(a = 0.5, b = 0.5)")
  expect_identical(format(prior_dirichlet(rep(0.25, 4))),
                   "Dirichlet(alpha = (0.25, 0.25, 0.25, 0.25))")
  expect_identical(format(prior_nix(0, 1, 3, 2)),
                   "N-Inv-Chisq(mu0 = 0, kappa0 = 1, nu0 = 3, sigma0 = 2)")
  expect_identical(format(prior_niw(c(0, 1), 1, 4, diag(2))),
                   paste0("N-Inv-Wishart(mu0 = (0, 1), kappa0 = 1, nu0 = 4, ",
                          "Lambda0 = (1, 0; 0, 1))"))
  expect_identical(format(prior_vague()), "Vague")
})

test_that("format.qdm_prior() honours digits and rejects unknown classes", {
  expect_identical(format(prior_nix(1 / 3, 1, 1, 1), digits = 3),
                   "N-Inv-Chisq(mu0 = 0.333, kappa0 = 1, nu0 = 1, sigma0 = 1)")
  fake <- structure(list(), class = c("qdm_prior_other", "qdm_prior"))
  expect_error(format(fake), "Unknown prior class.")
})

test_that("format.qdm_prior() describes the Gamma prior", {
  expect_identical(format(prior_gamma(0.1, 1)), "Gamma(shape = 0.1, rate = 1)")
})
