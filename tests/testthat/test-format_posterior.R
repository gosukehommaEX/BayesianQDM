test_that("format_posterior() formats Beta parameters", {
  expect_identical(format_posterior(ref_design("D1"), c(a = 7.5, b = 5.5)),
                   "Beta(a = 7.5, b = 5.5)")
  expect_identical(format_posterior(ref_design("D1"), c(a = 1 / 3, b = 2),
                                    digits = 2),
                   "Beta(a = 0.33, b = 2)")
})

test_that("format_posterior() gives the t distribution of a continuous mean", {
  des <- ref_design_cont("C1")
  expect_identical(format_posterior(des, c(mu = 18, kappa = 12, nu = 11,
                                           sigma = 25)),
                   "mean ~ t(df = 11, location = 18, scale = 7.217)")
  expect_identical(format_posterior(des, c(mu = -1 / 3, kappa = 4, nu = 2.5,
                                           sigma = 2), digits = 2),
                   "mean ~ t(df = 2.5, location = -0.33, scale = 1)")
})

test_that("format_posterior() gives the Dirichlet parameters", {
  expect_identical(format_posterior(ref_design_bin2("B1"),
                                    c("00" = 2.25, "01" = 1.25, "10" = 3.25,
                                      "11" = 4.25)),
                   "Dirichlet(alpha = (2.25, 1.25, 3.25, 4.25))")
})

test_that("format_posterior() gives the t distribution of two means", {
  # kappa = 10, nu = 11: the scale matrix is Lambda / 100
  post <- c(mu1 = 1, mu2 = -2, kappa = 10, nu = 11, L11 = 200, L12 = 50,
            L22 = 100)
  expect_identical(format_posterior(ref_design_cont2("K1"), post),
                   paste("mean ~ t(df = 10, location = (1, -2),",
                         "scale = (2.0, 0.5; 0.5, 1.0))"))
})
